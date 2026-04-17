# scripts/01b_replicate_table2.R
# Replicate Table 2 of Metrick & Schmelzing (2024): logit regressions of
# intervention-type dummies on country/crisis covariates.
#
# Paper sample (Table 2, p. 30): N = 273 crises with complete data on all
# regressors (INCOME log, POLITY, FX regime, DEBT/GDP, EXP/GDP). We use the
# same complete-case filter via r == 1 in our analysis CSV. Our N = 272
# (paper N = 273) is the documented +/-1 pipeline-provenance delta.
#
# Outcomes (in paper's column order, p. 30):
#   1. guarantees_d
#   2. lending_d
#   3. capital_injections_d
#   4. restructuring_d
#   5. asset_management_d
#   6. rules_d
#   7. other_d
#
# Specification (per paper text, p. 29-30):
#   logit(P(y=1)) = candidate + log_lagged_gdppc + polity + currency_regime
#                 + debt + exp + year
#
# Output:
#   data/processed/table2_logit_fits.rds  (list of 7 glm fits + metadata)
#
# Convergence note: rare-outcome columns (especially asset_management_d,
# rules_d, other_d on the 272-row complete-case sample) may produce glm
# convergence warnings. We surface (do NOT muffle) these warnings, capture
# them per spec, and persist them alongside the fits so 01c can tag affected
# rows in the comparison tibble.

source(here::here("R", "utils.R"))
load_libs()

df <- read.csv(here::here("data", "processed", "analysis_data.csv"))

df_complete <- df |> dplyr::filter(r == 1)

stopifnot(
  nrow(df_complete) >= 260L, nrow(df_complete) <= 280L,
  all(df_complete$year >= 1800)
)

outcomes <- c(
  "guarantees_d",
  "lending_d",
  "capital_injections_d",
  "restructuring_d",
  "asset_management_d",
  "rules_d",
  "other_d"
)

rhs <- "candidate + log_lagged_gdppc + polity + currency_regime + debt + exp + year"

fit_one <- function(outcome) {
  f <- as.formula(paste(outcome, "~", rhs))
  warns <- character()
  fit <- withCallingHandlers(
    glm(f, family = binomial("logit"), data = df_complete),
    warning = function(w) {
      warns[length(warns) + 1L] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  for (msg in warns) {
    warning(sprintf("[%s] %s", outcome, msg), call. = FALSE)
  }
  list(fit = fit, warnings = warns, converged = isTRUE(fit$converged))
}

fits_with_meta <- lapply(outcomes, fit_one)
names(fits_with_meta) <- outcomes

models <- lapply(fits_with_meta, `[[`, "fit")
warnings_per_spec <- lapply(fits_with_meta, `[[`, "warnings")
converged_per_spec <- vapply(fits_with_meta, `[[`, logical(1), "converged")

cat("\n=== Table 2 replication (logits, r == 1 sample) ===\n")
cat(sprintf("Complete-case N: %d\n", nrow(df_complete)))
for (i in seq_along(outcomes)) {
  o <- outcomes[i]
  fit <- models[[o]]
  beta <- coef(fit)["log_lagged_gdppc"]
  se <- summary(fit)$coefficients["log_lagged_gdppc", "Std. Error"]
  flag <- if (length(warnings_per_spec[[o]]) > 0L) " [WARN]" else ""
  cat(sprintf("%-22s INCOME = %+7.4f  (SE %5.4f)  N = %d  conv=%s%s\n",
              o, beta, se, nobs(fit),
              ifelse(converged_per_spec[o], "TRUE", "FALSE"),
              flag))
  if (length(warnings_per_spec[[o]]) > 0L) {
    for (msg in warnings_per_spec[[o]]) cat("    !", msg, "\n")
  }
}

dir.create(here::here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    models = models,
    warnings_per_spec = warnings_per_spec,
    converged_per_spec = converged_per_spec,
    rhs = rhs,
    sample_n = nrow(df_complete)
  ),
  here::here("data", "processed", "table2_logit_fits.rds")
)

cat(sprintf("\nSaved: %s\n",
            here::here("data", "processed", "table2_logit_fits.rds")))
