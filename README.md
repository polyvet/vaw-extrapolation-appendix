# Worked example — steepness- and precision-weighted pooling for extrapolation

Companion code for the methods paper describing the extrapolation step used in
the global estimates of violence against women: a steepness- and
precision-weighted log-pooling of a country's year-on-year change towards the
change implied by the global (parent) series, applied beyond a country's last
observed survey year.

Everything here is self-contained and runs on **simulated** data. No prevalence
estimate and no posterior draw from the fitted model is distributed in this
repository.

The repository holds **source only**. The data files and every output are
produced by running the scripts below, in order — the generator is seeded, so
the simulated data come back bit-identical on any machine.

---

## Quick start

    Rscript make_example_data.R                     # writes the simulated data
    Rscript extrapolation_example.R                 # runs the pooling algorithm
    Rscript extrapolation_example_paper_version.R   # compares against the
                                                    # equation as first written up

Requires R with dplyr, tidyr, purrr and ggplot2.

---

## Scripts

**make_example_data.R** — generates the simulated global (parent) series and the
three example countries. set.seed() makes it reproducible, so the .rds files can
be regenerated rather than trusted. Writes sim_global_draws.rds and
sim_country_X_draws.rds, sim_country_Y_draws.rds, sim_country_Z_draws.rds.

**extrapolation_example.R** — the worked example. Applies the pooling algorithm
to each of the three example countries and writes the post-extrapolation fit and
a plot for each: example_output_pooled_X.rds and example_plot_country_X.pdf, and
likewise for Y and Z.

**extrapolation_example_paper_version.R** — diagnostic companion. Identical
except for the final pooling equation, implemented as the Methods section
originally described it (year-specific global change) rather than as production
implements it (global change averaged over the observation window). The only
other difference is that the two lines constructing the time-averaged parent
delta are dropped, since that version does not use them. Writes
example_output_pooled_paper_version_X/Y/Z.rds and the three-panel comparison
figure example_plot_code_vs_paper.pdf.

---

## The simulated data

The four .rds files are produced by make_example_data.R from the parametric
model written out in that script. They reproduce the structural features the
pooling algorithm consumes, and nothing else:

1. each posterior draw is a **smooth curve across years** — the algorithm works
   on year-on-year changes on the logit scale, so draws that were independent
   noise across years would make those changes meaningless;
2. each country posterior is **narrow inside its observed window and fans out**
   after the last observed year, as a fitted spline does once it leaves the data;
3. each country trajectory is **steeper than the global one** after the last
   observation, so the steepness weight has something to act on;
4. the global series has a broad, roughly flat credible band, but the
   across-draw **variance of its year-on-year change grows** towards the end of
   the window. That second property, not the width of the band, is what the
   algorithm consumes, and it is why a window-averaged parent delta is more
   stable than a year-specific one.

Each series covers ages 15–49, 2000–2023, with 500 posterior draws. The three
countries cover the situations in which the algorithm behaves differently:

- **CTRY_X** (last observed year 2012) — rises to a plateau, then the
  unconstrained extrapolation turns downward. Long horizon, and a poorly
  determined start to the window.
- **CTRY_Y** (2016) — broadly flat through the observed window, then the
  unconstrained extrapolation drifts gently upward. Short horizon.
- **CTRY_Z** (2016) — a sustained decline that the unconstrained extrapolation
  carries on, at an easing pace. Wide posterior throughout.

In both Y and Z the post-horizon drift eases in from zero slope, so the mean
path leaves the observed window without a corner.

---

## What the example shows

Up to the last observed survey year the two trajectories are identical: the
algorithm does not touch the fitted period. Beyond it, the unconstrained
country-specific extrapolation drifts and its 95% uncertainty interval widens
steeply, while the pooled trajectory stays anchored.

**2023**

| | | median | 95% interval | width (pp) |
|---|---|---:|:---:|---:|
| CTRY_X | country-specific | 16.2% | 1.9 – 62.3 | 60.4 |
| | pooled | 22.0% | 18.2 – 25.8 | 7.6 |
| CTRY_Y | country-specific | 31.0% | 20.9 – 43.1 | 22.2 |
| | pooled | 27.7% | 24.9 – 31.0 | 6.1 |
| CTRY_Z | country-specific | 11.0% | 3.1 – 33.9 | 30.8 |
| | pooled | 15.0% | 12.2 – 18.4 | 6.2 |

In 2000 the two are identical in every country (CTRY_X 8.0%, width 6.7; CTRY_Y
34.0%, width 4.7; CTRY_Z 22.1%, width 22.5).

### As implemented vs as originally written up

The paper-version script quantifies the difference between the two forms of the
final pooling equation. The medians agree to within 0.1 pp; the year-specific
parent delta inflates the interval, because the across-draw variance of the
global year-on-year change grows towards the end of the window.

**2023**

| | | median | 95% interval | width (pp) |
|---|---|---:|:---:|---:|
| CTRY_X | as implemented | 22.0% | 18.2 – 25.8 | 7.6 |
| | paper equation | 21.9% | 16.9 – 27.2 | 10.3 |
| CTRY_Y | as implemented | 27.7% | 24.9 – 31.0 | 6.1 |
| | paper equation | 27.7% | 23.5 – 32.4 | 8.9 |
| CTRY_Z | as implemented | 15.0% | 12.2 – 18.4 | 6.2 |
| | paper equation | 14.9% | 11.7 – 19.5 | 7.8 |

The Methods section has since been corrected to the implemented form.

---

## Provenance of the code

The algorithm in extrapolation_example.R is reproduced **verbatim** from the
production script delta_log_pool_weff.R, including the configuration constants
(w0 = 0.50, beta_steep = 1.0, w_min = 0.05, w_max = 0.95, s0 = 1e-6,
eps_var = 1e-8) and set.seed(190575). Lines 199–333 of the production script,
with comments stripped, surrounding whitespace trimmed and blank lines dropped,
give a 98-line block with MD5 checksum 19ef95d0b0f4155cf0b920f2e8481098; that
block appears verbatim, in the same order, inside extrapolation_example.R.

What was removed is only the surrounding production plumbing that is not needed
for a single country: the loop over violence types and periods, the iso_keep
construction (used only by the alternative population-weighted parent branch;
production uses the model-derived global parent), the post-hoc flattening lists
and merge, the trailing block that writes the production summary grids, and
K_last, which is declared in the production configuration but never used.

The example reproduces production exactly. Run against the real production
posteriors for three countries that were not flattened, it returns the same 2023
pooled estimate as the production run to within 2.8e-15 percentage points.
