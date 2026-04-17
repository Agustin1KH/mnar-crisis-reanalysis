# scripts/01_replicate_table3.R
# Replicate Table 3 of Metrick & Schmelzing (2024): GDP gap regressions.
#
# Target numbers from the paper (p. 32, Table 3):
#   Combined sample: INCOME (log) coef = -0.2287, SE = 0.0370, N = 334
#     (when DEBT/EXP excluded — larger sample)
#   Combined sample with DEBT/EXP: INCOME coef = -0.2672, SE = 0.0426, N = 261
#   Canonical-only: INCOME coef = -0.2306 (no fiscal), -0.2815 (with fiscal)
#   Candidate-only: INCOME coef = -0.3621 (no fiscal), -0.2783 (with fiscal)
#
# Note: paper uses OLS with unadjusted SEs and acknowledges cross-sectional
# and time-series dependence. We report both OLS-naive and HC3-robust.

source(here::here("R", "utils.R"))
load_libs()

library(fixest)
library(sandwich)
library(lmtest)

df <- read.csv(here::here("data", "processed", "analysis_data.csv")) |>
  dplyr::filter(year >= 1800)  # paper: "Data begins in 1800"

# --- Spec A: no fiscal (larger sample) ---
# Regressors: candidate, INCOME (log), Polity, FX regime, year
spec_a <- gap_sum ~ candidate + log_lagged_gdppc + polity + currency_regime + year

fit_a_combined  <- lm(spec_a, data = df)
fit_a_canonical <- lm(spec_a, data = df |> dplyr::filter(candidate == 0))
fit_a_candidate <- lm(spec_a, data = df |> dplyr::filter(candidate == 1))

# --- Spec B: with fiscal variables (smaller, Table 3 col 2/4/6) ---
spec_b <- gap_sum ~ candidate + log_lagged_gdppc + polity + currency_regime +
  debt + exp + year

fit_b_combined  <- lm(spec_b, data = df)
fit_b_canonical <- lm(spec_b, data = df |> dplyr::filter(candidate == 0))
fit_b_candidate <- lm(spec_b, data = df |> dplyr::filter(candidate == 1))

# --- Pretty table ---
models <- list(
  "Combined, no fiscal" = fit_a_combined,
  "Combined, + fiscal" = fit_b_combined,
  "Canonical, no fiscal" = fit_a_canonical,
  "Canonical, + fiscal" = fit_b_canonical,
  "Candidate, no fiscal" = fit_a_candidate,
  "Candidate, + fiscal" = fit_b_candidate
)

tab <- modelsummary::modelsummary(
  models,
  stars = TRUE,
  gof_omit = "AIC|BIC|Log.Lik|F|RMSE",
  coef_rename = c(
    "log_lagged_gdppc" = "INCOME (log)",
    "currency_regime" = "FX regime",
    "polity" = "Polity",
    "candidate" = "Candidate",
    "debt" = "DEBT/GDP",
    "exp" = "EXP/GDP",
    "year" = "Year"
  ),
  output = "markdown"
)

# bugfix (output only, not a spec change): under modelsummary >= 2.x the
# markdown return is an S4 tinytable object, so cat() fails with
# "argument 1 (type 'S4') cannot be handled by 'cat'". Coerce to character.
writeLines(as.character(tab),
           here::here("output", "tables", "01_table3_replication.md"))

# --- Print side-by-side to console for quick eyeballing ---
cat("\n=== Table 3 replication ===\n")
for (nm in names(models)) {
  fit <- models[[nm]]
  beta <- coef(fit)["log_lagged_gdppc"]
  se <- summary(fit)$coefficients["log_lagged_gdppc", "Std. Error"]
  cat(sprintf("%-26s  INCOME = %6.4f  (SE %5.4f)  N = %d\n",
              nm, beta, se, nobs(fit)))
}

# --- Save main estimates for downstream scripts to compare against ---
saveRDS(
  list(models = models, spec_a = spec_a, spec_b = spec_b),
  here::here("data", "processed", "table3_ols_fits.rds")
)
