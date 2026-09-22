# ===============================================================
# Generate the simulated example data used by the appendix example.
#
# NOTHING IN THIS SCRIPT IS DERIVED FROM WHO DATA. The global parent and
# all example countries are simulated from the parametric model below,
# so the example can be distributed without releasing any prevalence
# estimate.
#
# The simulation reproduces the structural features the pooling
# algorithm depends on, and nothing else:
#
#   1. Each posterior draw is a SMOOTH curve across years. The algorithm
#      works on year-on-year changes on the logit scale, so draws that
#      were independent noise across years would make those changes
#      meaningless.
#   2. Each country posterior is narrow inside its observed window and
#      fans out after the last observed year, as a fitted spline does
#      once it leaves the data.
#   3. Each country trajectory is steeper than the global one after the
#      last observation, so the steepness weight has something to act on.
#   4. The global series has a broad, roughly flat credible band across
#      the window, but the across-draw variance of its YEAR-ON-YEAR
#      CHANGE grows towards the end of the window. That second property,
#      not the width of the band itself, is what the pooling algorithm
#      consumes, and it is why a window-averaged parent delta is more
#      stable than a year-specific one.
#
# Three example countries are generated, covering the situations in
# which the algorithm behaves differently:
#
#   CTRY_X  rises to a plateau, then the unconstrained extrapolation
#           turns downward. Long horizon (11 years).
#   CTRY_Y  broadly flat through the observed window, then the
#           unconstrained extrapolation drifts upward. Short horizon.
#   CTRY_Z  a sustained decline that the unconstrained extrapolation
#           carries on, at a gently easing pace. Also has a poorly
#           determined start to the window.
#
# In both Y and Z the post-horizon drift eases in from zero slope, so the
# mean path leaves the observed window without a corner.
#
# Running this script writes:
#   sim_global_draws.rds
#   sim_country_X_draws.rds, sim_country_Y_draws.rds, sim_country_Z_draws.rds
#
# set.seed() makes the output reproducible, so the files can be
# regenerated identically rather than distributed.
# ===============================================================

set.seed(20260922)

YEARS   <- 2000:2023
N_DRAWS <- 500
AGE     <- "15-49"
T_MID   <- mean(YEARS)
HALF    <- (max(YEARS) - min(YEARS)) / 2
NY      <- length(YEARS)

logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# time bases. tc is centred (-1..1): its year-on-year increment is constant,
# so it sets a floor on the variance of the global change. u4 is a quartic in
# forward time: its increment is near zero early and largest at the end of the
# window, which is what makes the variance of the global change grow with time
# without inflating the width of the band itself.
tc     <- (YEARS - T_MID) / HALF                              # -1 .. 1
u      <- (YEARS - min(YEARS)) / (max(YEARS) - min(YEARS))    # 0 .. 1
b_lin  <- tc
b_u4   <- u^4 - mean(u^4)

# ---------------------------------------------------------------
# Global (parent) series
#   mean path: 32% in 2000 declining to 27% in 2023, on the logit scale.
#   The level spread is dominated by a draw-level intercept, so the band
#   is wide and roughly flat; the change-variance grows through b_u4.
# ---------------------------------------------------------------
g_start <- logit(0.32)
g_end   <- logit(0.27)
g_mean  <- seq(g_start, g_end, length.out = NY)

g_a <- rnorm(N_DRAWS, 0, 0.190)   # level          - sets the width of the band
g_b <- rnorm(N_DRAWS, 0, 0.074)   # linear         - floor on the change-variance
g_c <- rnorm(N_DRAWS, 0, 0.140)   # late curvature - growth in the change-variance

eta_g <- outer(rep(1, NY), g_a) + outer(b_lin, g_b) + outer(b_u4, g_c)
eta_g <- sweep(eta_g, 1, g_mean, "+")
p_g   <- inv_logit(eta_g)

# ---------------------------------------------------------------
# Country generator
#   mean_path : vector of length NY, on the logit scale
#   last_yr   : last observed year
#   sd_in     : posterior spread inside the observed window
#   a_out,b_out,p_out : growth of the spread beyond the last observation,
#                       sd_out = a_out * h + b_out * h^p_out, h = years since
#   early     : extra spread at the very start of the window (0 = none)
# ---------------------------------------------------------------
make_country <- function(mean_path, last_yr, sd_in, a_out, b_out, p_out, early = 0){
  h      <- pmax(0, YEARS - last_yr)
  sd_out <- a_out * h + b_out * h^p_out
  sd_beg <- early * pmax(0, 1 - (YEARS - min(YEARS)) / 6)^2   # decays over ~6 years

  z1 <- rnorm(N_DRAWS, 0, 1); z2 <- rnorm(N_DRAWS, 0, 1)
  z3 <- rnorm(N_DRAWS, 0, 1); z4 <- rnorm(N_DRAWS, 0, 1)

  eta <- outer(rep(1, NY), z1) * sd_in +
         outer(b_lin, z2) * sd_in * 0.6 +
         outer(sd_out, z3) +
         outer(sd_beg, z4)
  eta <- sweep(eta, 1, mean_path, "+")
  inv_logit(eta)
}

smooth_ramp <- function(from_year, to_year){
  r <- (YEARS - from_year) / (to_year - from_year)
  r <- pmin(pmax(r, 0), 1)
  3 * r^2 - 2 * r^3
}

# Post-horizon drift that eases in. h is years since the last observation.
# The value and the slope are both zero at h = 0, so the mean path leaves the
# observed window without a corner; the slope approaches `rate` only slowly.
drift <- function(h, rate, tau) rate * h^2 / (h + tau)

# --- CTRY_X : rise to a plateau, then a downward turn, long horizon ------
LY_X   <- 2012
h_X    <- pmax(0, YEARS - LY_X)
mean_X <- logit(0.08) + (logit(0.24) - logit(0.08)) * smooth_ramp(2000, LY_X) -
          0.050 * h_X + 0.0011 * h_X^2
p_X <- make_country(mean_X, LY_X, sd_in = 0.055, a_out = 0.040, b_out = 0.0160, p_out = 1.6,
                    early = 0.22)

# --- CTRY_Y : broadly flat, then an upward drift, short horizon ----------
LY_Y   <- 2016
h_Y    <- pmax(0, YEARS - LY_Y)
mean_Y <- logit(0.34) + (logit(0.29) - logit(0.34)) * smooth_ramp(2000, LY_Y) +
          drift(h_Y, rate = 0.022, tau = 3)
p_Y <- make_country(mean_Y, LY_Y, sd_in = 0.048, a_out = 0.020, b_out = 0.0025, p_out = 2.0)

# --- CTRY_Z : sustained decline that keeps going, weak start -------------
LY_Z   <- 2016
h_Z    <- pmax(0, YEARS - LY_Z)
mean_Z <- logit(0.22) + (logit(0.23) - logit(0.22)) * smooth_ramp(2000, 2005) +
          (logit(0.16) - logit(0.23)) * smooth_ramp(2005, LY_Z) -
          drift(h_Z, rate = 0.075, tau = 3)
p_Z <- make_country(mean_Z, LY_Z, sd_in = 0.105, a_out = 0.055, b_out = 0.0065,
                    p_out = 2.0, early = 0.30)

# ---------------------------------------------------------------
# Assemble in the layout the pooling code expects
# ---------------------------------------------------------------
draw_names <- paste0("sample_", seq_len(N_DRAWS))

write_country <- function(p, iso, last_yr, file){
  df <- data.frame(country = iso, age_group = AGE, time = as.numeric(YEARS),
                   last_year = last_yr, allow_extrapolation = 0,
                   stringsAsFactors = FALSE)
  df <- cbind(df, as.data.frame(p)); names(df)[6:ncol(df)] <- draw_names
  saveRDS(df, file); df
}

global <- data.frame(country = "GLOBAL", age_group = AGE, time = as.numeric(YEARS),
                     stringsAsFactors = FALSE)
global <- cbind(global, as.data.frame(p_g)); names(global)[4:ncol(global)] <- draw_names
saveRDS(global, "sim_global_draws.rds")

write_country(p_X, "CTRY_X", LY_X, "sim_country_X_draws.rds")
write_country(p_Y, "CTRY_Y", LY_Y, "sim_country_Y_draws.rds")
write_country(p_Z, "CTRY_Z", LY_Z, "sim_country_Z_draws.rds")

cat("Wrote sim_global_draws.rds and sim_country_{X,Y,Z}_draws.rds\n\n")
summ <- function(m, lab){
  q <- apply(m, 1, quantile, c(.025, .5, .975))
  data.frame(series = lab, year = YEARS,
             median = round(100 * q[2, ], 1), lower = round(100 * q[1, ], 1),
             upper  = round(100 * q[3, ], 1), width = round(100 * (q[3, ] - q[1, ]), 1))
}
out <- rbind(summ(p_g, "global"), summ(p_X, "CTRY_X"), summ(p_Y, "CTRY_Y"), summ(p_Z, "CTRY_Z"))
print(out[out$year %in% c(2000, 2010, 2016, 2020, 2023), ], row.names = FALSE)
