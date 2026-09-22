# ===============================================================
# Delta pooling with steepness-aware w_eff + precision log-pooling
#
# SINGLE-COUNTRY EXAMPLE VERSION
#
# Derived from the production script delta_log_pool_weff.R by removing
# the parts that only apply to a full estimation round. The algorithm
# itself (sections 1-7 below) is unchanged.
#
# Removed relative to the production script:
#   - the loop over IPV types and periods (one outcome, one period here)
#   - the WHO flattening ISO3 lists and replace_pooled_with_flat(),
#     together with the post-hoc "merge WHO flattening" step
#   - the iso_keep construction (previous_estimates_lifetime.xlsx and the
#     filtered_data_*.rds files), which is used only by the alternative
#     population-weighted parent; production uses the model-derived
#     global parent, which is what is used here
#   - the trailing block that writes the production SUMMARY grids
#   - K_last, which is declared in the production config but never used
#
# Inputs (all simulated - run make_example_data.R first to create them):
#   sim_country_X_draws.rds    - CTRY_X, ages 15-49, 2000-2023, 500 draws,
#   sim_country_Y_draws.rds      last_year 2012 / 2016 / 2016,
#   sim_country_Z_draws.rds      allow_extrapolation = 0
#   sim_global_draws.rds       - the global (parent) posterior the three
#                                countries share, ages 15-49, 2000-2023,
#                                500 draws
# ===============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

options(dplyr.summarise.inform = FALSE)

# ---------- Config ----------
ipv <- "PSIPV"
per <- "lifetime"

# Base weight and w_eff controls
w0            <- 0.50
beta_steep    <- 1.0
w_min         <- 0.05
w_max         <- 0.95
s0            <- 1e-6

# Floors
eps_var       <- 1e-8

set.seed(190575)

# ---------- Example inputs ----------
# The three example countries, and the global (parent) series they share.
# All are simulated: run make_example_data.R first to create them.
COUNTRY_FILES <- c(
  X = "sim_country_X_draws.rds",   # rise to a plateau, then a downward turn
  Y = "sim_country_Y_draws.rds",   # broadly flat, then an upward drift
  Z = "sim_country_Z_draws.rds"    # a sustained decline that keeps going
)
file_parent_global <- "sim_global_draws.rds"

# ---------- Helpers ----------
inv_logit  <- function(x) 1/(1+exp(-x))
clip01     <- function(p, eps=1e-6) pmin(pmax(p, eps), 1-eps)


process_one <- function(ipv, per){
  message("\n========== ", ipv, " · ", per, " ==========")

  if (!file.exists(file_in)) {
    message("✋ Input not found: ", file_in, " — skipping.")
    return(invisible(NULL))
  }

  # ---- Load country grid ----
  grid_ctry <- readRDS(file_in)
  draw_cols <- grep("^sample_", names(grid_ctry), value = TRUE)
  if (!length(draw_cols)) {
    message("✋ No draw columns found in ", file_in, " — skipping.")
    return(invisible(NULL))
  }

  # ===============================================================
  # 1) Parent draws: model-derived global posterior
  # ===============================================================
  message("→ Using MODEL-DERIVED global parent: ", file_parent_global)
  parent_draws <- readRDS(file_parent_global) %>%
    mutate(age_group = as.character(age_group)) %>%
    select(age_group, time, all_of(draw_cols)) %>%
    pivot_longer(all_of(draw_cols), names_to = "draw", values_to = "p_parent")

  # ---- one global parent delta per draw (age-specific, time-invariant) ----
  parent_d_global <- parent_draws %>%
    arrange(age_group, draw, time) %>%
    group_by(age_group, draw) %>%
    mutate(eta_p = qlogis(clip01(p_parent)),
           d_p   = eta_p - lag(eta_p)) %>%
    summarise(d_parent_global = mean(d_p, na.rm = TRUE), .groups = "drop")

  # Per age-group summary
  parent_dg_sum <- parent_d_global %>%
    group_by(age_group) %>%
    summarise(
      n_draws = sum(!is.na(d_parent_global)),
      mean_d  = ifelse(n_draws > 0, mean(d_parent_global, na.rm = TRUE), NA_real_),
      sd_d    = ifelse(n_draws > 1, sd(d_parent_global, na.rm = TRUE), NA_real_),
      q05     = ifelse(n_draws > 0, quantile(d_parent_global, 0.05, na.rm = TRUE), NA_real_),
      q95     = ifelse(n_draws > 0, quantile(d_parent_global, 0.95, na.rm = TRUE), NA_real_)
    ) %>%
    arrange(age_group)

  print(parent_dg_sum, n = Inf)

  # ===============================================================
  # 2) Long + logits + per-draw deltas for country & parent
  # ===============================================================
  grid_long <- grid_ctry %>%
    select(country, age_group, time, last_year, allow_extrapolation, all_of(draw_cols)) %>%
    pivot_longer(all_of(draw_cols), names_to = "draw", values_to = "p_country") %>%
    left_join(parent_draws, by = c("age_group","time","draw")) %>%
    mutate(
      eta_country = qlogis(clip01(p_country)),
      eta_parent  = qlogis(clip01(p_parent))
    ) %>%
    arrange(country, age_group, draw, time) %>%
    group_by(country, age_group, draw) %>%
    mutate(
      d_country = eta_country - dplyr::lag(eta_country),
      d_parent  = eta_parent  - dplyr::lag(eta_parent)
    ) %>%
    ungroup()

  grid_long <- grid_long %>%
    left_join(parent_d_global, by = c("age_group","draw")) %>%
    mutate(d_parent_use = d_parent_global)

  rm(parent_draws); gc()

  # ===============================================================
  # 3) Time-specific steepness -> w_eff_t per (country, age_group, time)
  #    Compare |delta|country vs |delta|parent at EACH time t
  # ===============================================================

  # (we already have d_country and d_parent computed per draw, per time)

  abs_delta_med_ctry <- grid_long %>%
    filter(!is.na(d_country)) %>%
    group_by(country, age_group, time) %>%
    summarise(abs_d_c_med = median(abs(d_country), na.rm = TRUE),
              .groups = "drop")

  abs_delta_med_par <- grid_long %>%
    filter(!is.na(d_parent)) %>%
    group_by(age_group, time) %>%
    summarise(abs_d_p_med = median(abs(d_parent), na.rm = TRUE),
              .groups = "drop")

  # Time-specific ratio and effective base weight
  w_eff_tbl_t <- abs_delta_med_ctry %>%
    left_join(abs_delta_med_par, by = c("age_group","time")) %>%
    mutate(
      ratio   = (abs_d_p_med + s0) / (abs_d_c_med + s0),
      w_eff_t = pmin(pmax(w0 * (ratio ^ beta_steep), w_min), w_max)
    ) %>%
    select(country, age_group, time, w_eff_t)

  # Fill any missing time points with w0
  w_eff_tbl_t <- grid_ctry %>%
    distinct(country, age_group, time) %>%
    left_join(w_eff_tbl_t, by = c("country","age_group","time")) %>%
    mutate(
      w_eff_t = ifelse(is.na(w_eff_t), w0, w_eff_t)
    ) %>%
    select(country, age_group, time, w_eff_t)

  # ===============================================================
  # 4) Precision-aware log-pooling weights (per time)
  # ===============================================================
  var_by_time <- grid_long %>%
    group_by(country, age_group, time) %>%
    summarise(
      var_c = var(d_country, na.rm = TRUE),
      var_p = var(d_parent,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      var_c = ifelse(is.na(var_c) | var_c < eps_var, eps_var, var_c),
      var_p = ifelse(is.na(var_p) | var_p < eps_var, eps_var, var_p)
    )


  grid_long <- grid_long %>%
    left_join(var_by_time, by = c("country","age_group","time")) %>%
    left_join(w_eff_tbl_t, by = c("country","age_group","time")) %>%
    mutate(
      # parent precision share under log-pooling (now using w_eff_t)
      alpha_parent = ((1 - w_eff_t) / var_p) / (((w_eff_t) / var_c) + ((1 - w_eff_t) / var_p)),
      alpha_parent = pmin(pmax(alpha_parent, 0), 1),
      d_pool_logpool_weff = (1 - alpha_parent) * d_country + alpha_parent * d_parent_use
    )

  rm(var_by_time, w_eff_tbl_t); gc()

  # ===============================================================
  # 5) Rebuild eta recursively after last_year (only when blocked)
  # ===============================================================
  grid_long <- grid_long %>%
    arrange(country, age_group, draw, time) %>%
    group_by(country, age_group, draw) %>%
    mutate(
      eta_pooled_logpool_weff = {
        out <- eta_country
        ly  <- unique(last_year)[1]
        pool_this <- isTRUE(unique(allow_extrapolation)[1] == 0)
        if (pool_this) {
          idx <- which(time > ly)
          if (length(idx) > 0) {
            for (k in seq_along(idx)) {
              t_i <- idx[k]
              out[t_i] <- out[t_i - 1] + d_pool_logpool_weff[t_i]
            }
          }
        }
        out
      }
    ) %>%
    ungroup() %>%
    mutate(
      p_country             = inv_logit(eta_country),
      p_pooled_logpool_weff = inv_logit(eta_pooled_logpool_weff)
    )

  # ===============================================================
  # 6) Reshape & Save pooled draws (ONLY country & logpool_weff)
  # ===============================================================
  wide_pooled <- grid_long %>%
    select(country, age_group, time, draw,
           p_country, p_pooled_logpool_weff) %>%
    pivot_wider(names_from = draw,
                values_from = c(p_country, p_pooled_logpool_weff))

  grid_out <- grid_ctry %>%
    select(-all_of(draw_cols)) %>%
    left_join(wide_pooled, by = c("country","age_group","time"))

  saveRDS(grid_out, file_out)
  message("✅ Saved: ", file_out)

  # ===============================================================
  # 7) Plot: Country vs log-pool-w_eff (15-49)
  # ===============================================================
  grid_1549 <- grid_out %>% filter(as.character(age_group) == "15-49")

  draw_ctry <- grep("^p_country_",               names(grid_1549), value = TRUE)
  draw_log  <- grep("^p_pooled_logpool_weff_",   names(grid_1549), value = TRUE)

  if (!length(draw_ctry) || !length(draw_log)) {
    message("✋ Missing draw columns for plotting in ", ipv, " · ", per, " — skipping PDF.")
    rm(grid_ctry, grid_long, wide_pooled, grid_out, grid_1549); gc()
    return(invisible(NULL))
  }

  summarise_draws <- function(df, draw_cols, label){
    df %>%
      select(country, time, last_year, all_of(draw_cols)) %>%
      pivot_longer(all_of(draw_cols), names_to = "draw", values_to = "p") %>%
      group_by(country, time, last_year) %>%
      summarise(
        median = median(p, na.rm = TRUE),
        lci    = quantile(p, 0.025, na.rm = TRUE),
        uci    = quantile(p, 0.975, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(type = label)
  }

  df_ctry <- summarise_draws(grid_1549, draw_ctry, "Country")
  df_log  <- summarise_draws(grid_1549, draw_log,  "Pooled (log-pool, w_eff)")
  plot_df <- bind_rows(df_ctry, df_log)

  plot_one <- function(df_country){
    iso  <- unique(df_country$country)
    last <- unique(df_country$last_year)
    ggplot(df_country, aes(x = time, y = median, colour = type, linetype = type)) +
      geom_ribbon(aes(ymin = lci, ymax = uci, fill = type),
                  alpha = .20, colour = NA, show.legend = FALSE) +
      geom_line(size = .9) +
      scale_colour_manual(values = c(
        "Country"                 = "orange",
        "Pooled (log-pool, w_eff)"= "darkgreen"
      )) +
      scale_linetype_manual(values = c(
        "Country"                 = "solid",
        "Pooled (log-pool, w_eff)"= "dashed"
      )) +
      geom_vline(xintercept = last, colour = "grey40") +
      labs(title = paste(iso, "–", ipv, "·", per),
           x = "Year", y = "Prevalence 15–49") +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.title = element_blank())
  }

  countries <- sort(unique(plot_df$country))
  pdf(pdf_out, width = 7, height = 5)
  walk(countries, ~ print(plot_one(filter(plot_df, country == .x))))
  dev.off()
  message("✅ PDF saved: ", pdf_out)

  invisible(plot_df)
}

# ---------- Run the example for each country ----------
all_res <- list()
for (k in names(COUNTRY_FILES)) {
  file_in  <- COUNTRY_FILES[[k]]
  file_out <- sprintf("example_output_pooled_%s.rds", k)
  pdf_out  <- sprintf("example_plot_country_%s.pdf", k)
  all_res[[k]] <- process_one(ipv, per)
}

summary_tbl <- bind_rows(all_res) %>%
  filter(time %in% range(time)) %>%
  mutate(across(c(median, lci, uci), ~ round(100 * .x, 1)),
         width = round(uci - lci, 1)) %>%
  arrange(country, time, type) %>%
  select(country, time, last_year, type, median, lci, uci, width)

cat("\n=== Worked example: first and final year, each country ===\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
