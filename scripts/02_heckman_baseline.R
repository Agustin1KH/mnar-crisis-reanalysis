# scripts/02_heckman_baseline.R
# Bayesian-flavoured baseline: Heckman Type II Tobit via `sampleSelection`.
# Frequentist MLE, but strictly more general than the paper's OLS: models
# the selection into the complete-case sample jointly with the GDP gap outcome.
#
# The argument:
#   * Table 3 OLS uses 260 of 910 crises. The 650 dropped crises are
#     systematically older, smaller-country, less-documented.
#   * Under MNAR (the Heckman ρ ≠ 0 case), OLS on the complete-case sample
#     is biased for the population parameter.
#   * `selection()` with method = "ml" jointly estimates a probit for R = 1
#     (completeness) and the OLS for gap_sum, with bivariate normal errors
#     linked by correlation ρ.
#
# Exclusion restriction (variable in selection but not outcome):
#   n_chron_clean — number of canonical chronologies (RR, ST, LV, BVX) that
#   list the crisis. Affects the probability that a crisis has full regressor
#   data (more-covered crises have more archival follow-up) but, given the
#   economic covariates, plausibly does not enter the GDP gap directly.

source(here::here("R", "utils.R"))
load_libs()

library(sampleSelection)

df <- read.csv(here::here("data", "processed", "analysis_data.csv")) |>
  dplyr::filter(year >= 1800)

# --- Selection equation: predicts R = 1 (complete case) ---
# Include the outcome-equation regressors plus the exclusion restriction.
# Selection equation needs variables observed for ALL crises — so no Exp/Debt
# here (they are what's being selected on). Use year, candidate, maddison_priority,
# and the chronology coverage as the exclusion restriction.
#
# Key identification: n_chron_clean must affect R but not (conditional on the
# rest) the gap_sum outcome. This is the testable-but-not-from-data piece.
sel_eq <- r ~ candidate + year + maddison_priority + n_chron_clean

# --- Outcome equation: the Table 3 Spec A (no fiscal, larger sample) ---
out_eq_a <- gap_sum ~ candidate + log_lagged_gdppc + polity + currency_regime + year

# --- Fit by MLE ---
fit_heckman_a <- selection(
  selection = sel_eq,
  outcome   = out_eq_a,
  data      = df,
  method    = "ml"
)

summary(fit_heckman_a)

# --- Two-step as a sanity check (classical Heckit) ---
fit_heckman_a_2step <- selection(
  selection = sel_eq,
  outcome   = out_eq_a,
  data      = df,
  method    = "2step"
)

# --- Full-spec outcome (with fiscal variables) ---
out_eq_b <- gap_sum ~ candidate + log_lagged_gdppc + polity + currency_regime +
  debt + exp + year

fit_heckman_b <- selection(
  selection = sel_eq,
  outcome   = out_eq_b,
  data      = df,
  method    = "ml"
)

# --- Pull out the key coefficient: INCOME on gap_sum ---
extract_income <- function(fit, label) {
  co <- summary(fit)$estimate
  row <- grep("log_lagged_gdppc", rownames(co))
  if (length(row) == 0) return(NULL)
  co_row <- co[row[1], ]
  cat(sprintf("%-22s  INCOME = %7.4f  (SE %5.4f)  rho = %6.3f\n",
              label,
              co_row["Estimate"],
              co_row["Std. error"],
              if ("rho" %in% rownames(co)) co["rho", "Estimate"] else NA))
}

cat("\n=== Heckman vs OLS on INCOME (log_lagged_gdppc → gap_sum) ===\n")
ols_ref <- readRDS(here::here("data", "processed", "table3_ols_fits.rds"))
cat(sprintf("%-22s  INCOME = %7.4f  (SE %5.4f)\n",
            "OLS complete-case A",
            coef(ols_ref$models[["Combined, no fiscal"]])["log_lagged_gdppc"],
            summary(ols_ref$models[["Combined, no fiscal"]])$coefficients["log_lagged_gdppc", "Std. Error"]))
extract_income(fit_heckman_a, "Heckman MLE A")
extract_income(fit_heckman_a_2step, "Heckman 2-step A")
cat(sprintf("%-22s  INCOME = %7.4f  (SE %5.4f)\n",
            "OLS complete-case B",
            coef(ols_ref$models[["Combined, + fiscal"]])["log_lagged_gdppc"],
            summary(ols_ref$models[["Combined, + fiscal"]])$coefficients["log_lagged_gdppc", "Std. Error"]))
extract_income(fit_heckman_b, "Heckman MLE B")

# --- Likelihood-ratio test: ρ = 0 (MAR) vs ρ ≠ 0 (MNAR) ---
# If we cannot reject ρ = 0, the MNAR correction is not needed for this outcome.
rho_hat <- fit_heckman_a$estimate["rho"]
rho_se  <- sqrt(diag(vcov(fit_heckman_a)))["rho"]
cat(sprintf("\nMAR-vs-MNAR test (spec A): rho = %.3f (SE %.3f), |t| = %.2f\n",
            rho_hat, rho_se, abs(rho_hat / rho_se)))

# --- Save ---
saveRDS(
  list(a_ml = fit_heckman_a, a_2step = fit_heckman_a_2step, b_ml = fit_heckman_b),
  here::here("data", "processed", "heckman_fits.rds")
)
