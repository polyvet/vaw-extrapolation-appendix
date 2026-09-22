# ===============================================================
# Delta pooling with steepness-aware w_eff + precision log-pooling
#
# SINGLE-COUNTRY EXAMPLE — IMPLEMENTED EXACTLY AS DESCRIBED IN THE PAPER
#
# This is a diagnostic companion to extrapolation_example.R.
# Everything is identical except the final pooling equation.
#
#   As implemented in production (delta_log_pool_weff.R):
#       d*_t = (1 - alpha_p,t) * d_c,t + alpha_p,t * d_bar_p
#     where d_bar_p is the parent delta AVERAGED OVER ALL YEARS, i.e. one
#     time-invariant constant per age group per posterior draw:
#       d_bar_p = (1/(T-1)) * sum_t d_p,t = (eta_p,T - eta_p,1) / (T - 1)
#
#   As written in the paper's Methods:
#       d*_t = (1 - alpha_p,t) * d_c,t + alpha_p,t * d_p,t
#     using the YEAR-SPECIFIC parent delta.
#
# This script implements the second form, so the difference between the
# two can be quantified for the worked example. The steepness ratio and
# the precision weights are unchanged; they already use the year-specific
# parent delta in both versions.
#
# Inputs: sim_country_{X,Y,Z}_draws.rds and sim_global_draws.rds, all
# produced by make_example_data.R. Run extrapolation_example.R first, so
# the as-implemented outputs are on disk to compare against.
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
COUNTRY_FILES <- c(
  X = "sim_country_X_draws.rds",
  Y = "sim_country_Y_draws.rds",
  Z = "sim_country_Z_draws.rds"
)
file_parent_global <- "sim_global_draws.rds"

# ---------- Helpers ----------
inv_logit  <- function(x) 1/(1+exp(-x))
clip01     <- function(p, eps=1e-6) pmin(pmax(p, eps), 1-eps)


process_one <- function(ipv, per){
  message("\n========== ", ipv, " · ", per, "  (paper equation) ==========")

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

  # NOTE: the production script computes here a single time-averaged parent
  # delta per (age_group, draw) and uses it in section 4. The paper's
  # equation uses the year-specific d_p,t instead, so that step is absent.

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

  rm(parent_draws); gc()

  # ===============================================================
  # 3) Time-specific steepness -> w_eff_t per (country, age_group, time)
  # ===============================================================
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

  w_eff_tbl_t <- abs_delta_med_ctry %>%
    left_join(abs_delta_med_par, by = c("age_group","time")) %>%
    mutate(
      ratio   = (abs_d_p_med + s0) / (abs_d_c_med + s0),
      w_eff_t = pmin(pmax(w0 * (ratio ^ beta_steep), w_min), w_max)
    ) %>%
    select(country, age_group, time, w_eff_t)

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
      alpha_parent = ((1 - w_eff_t) / var_p) / (((w_eff_t) / var_c) + ((1 - w_eff_t) / var_p)),
      alpha_parent = pmin(pmax(alpha_parent, 0), 1),
      # >>> THE ONE LINE THAT DIFFERS FROM THE PRODUCTION SCRIPT <<<
      # production: alpha_parent * d_parent_use   (time-averaged parent delta)
      # paper:      alpha_parent * d_parent       (year-specific parent delta)
      d_pool_logpool_weff = (1 - alpha_parent) * d_country + alpha_parent * d_parent
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
  # 6) Reshape & Save pooled draws
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

  invisible(grid_out)
}


# ---------------------------------------------------------------
# Run the paper-equation version for each country, and compare
# ---------------------------------------------------------------
summarise_draws <- function(df, prefix, label){
  cols <- grep(prefix, names(df), value = TRUE)
  df %>%
    filter(as.character(age_group) == "15-49") %>%
    select(country, time, last_year, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = "draw", values_to = "p") %>%
    group_by(country, time, last_year) %>%
    summarise(
      median = median(p, na.rm = TRUE),
      lci    = quantile(p, 0.025, na.rm = TRUE),
      uci    = quantile(p, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(type = label)
}

plot_all <- list()
for (k in names(COUNTRY_FILES)) {
  file_in  <- COUNTRY_FILES[[k]]
  file_out <- sprintf("example_output_pooled_paper_version_%s.rds", k)
  grid_paper <- process_one(ipv, per)

  file_as_coded <- sprintf("example_output_pooled_%s.rds", k)
  if (is.null(grid_paper) || !file.exists(file_as_coded)) {
    message("Run extrapolation_example.R first to produce ", file_as_coded)
    next
  }
  grid_coded <- readRDS(file_as_coded)
  plot_all[[k]] <- bind_rows(
    summarise_draws(grid_coded, "^p_country_",             "Country (unconstrained)"),
    summarise_draws(grid_coded, "^p_pooled_logpool_weff_", "Pooled - as implemented"),
    summarise_draws(grid_paper, "^p_pooled_logpool_weff_", "Pooled - paper equation")
  )
}

plot_df <- bind_rows(plot_all)

if (nrow(plot_df)) {
  vl <- plot_df %>% distinct(country, last_year)

  p <- ggplot(plot_df, aes(x = time, y = 100 * median, colour = type, linetype = type)) +
    geom_ribbon(aes(ymin = 100 * lci, ymax = 100 * uci, fill = type),
                alpha = .14, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = .8) +
    geom_vline(data = vl, aes(xintercept = last_year), colour = "grey40", linewidth = .4) +
    facet_wrap(~ country, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = c(
      "Country (unconstrained)" = "orange",
      "Pooled - as implemented" = "darkgreen",
      "Pooled - paper equation" = "purple")) +
    scale_fill_manual(values = c(
      "Country (unconstrained)" = "orange",
      "Pooled - as implemented" = "darkgreen",
      "Pooled - paper equation" = "purple")) +
    scale_linetype_manual(values = c(
      "Country (unconstrained)" = "solid",
      "Pooled - as implemented" = "dashed",
      "Pooled - paper equation" = "dotdash")) +
    labs(title = "Pooling as implemented vs as written in the Methods",
         x = "Year", y = "Prevalence 15\u201349 (%)") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))

  pdf("example_plot_code_vs_paper.pdf", width = 11, height = 4.6)
  print(p)
  dev.off()
  message("\u2705 PDF saved: example_plot_code_vs_paper.pdf")

  cat("\n=== Final year (2023): as implemented vs paper equation ===\n")
  print(
    plot_df %>% group_by(country) %>% filter(time == max(time)) %>% ungroup() %>%
      mutate(across(c(median, lci, uci), ~ round(100 * .x, 1)),
             width = round(uci - lci, 1)) %>%
      select(country, type, median, lci, uci, width) %>% as.data.frame(),
    row.names = FALSE
  )
}
