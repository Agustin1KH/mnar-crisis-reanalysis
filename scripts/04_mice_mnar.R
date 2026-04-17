# scripts/04_mice_mnar.R
# Multiple imputation under MNAR using a Heckman selection model (Galimard
# et al. 2016, 2018). Produces m imputed datasets, each runs the Table 3
# regression, and we pool the m estimates via Rubin's rules.
#
# Why this complements 02_heckman_baseline.R:
#   * The selection() baseline corrects the outcome equation directly but
#     doesn't give you imputed values for the dropped 650 crises.
#   * MICE-MNAR produces actual imputed gap_sum, log_lagged_gdppc, etc. for
#     those crises — so you can then run any downstream analysis on the
#     "complete" dataset and compare.
#   * Useful for also imputing the currency_regime variable (60% missing),
#     which is the binding constraint on the complete-case sample.
#
# TODO:
#   1. library(miceMNAR) + library(mice)
#   2. Set up mice() with method = "hecknorm" for continuous outcomes with
#      MNAR, method = "heckprob" for the binary intervention dummies.
#   3. Run m = 20 imputations with the same exclusion restriction (n_chron_clean).
#   4. Refit Table 3 on each imputed dataset, pool via mice::pool().
#   5. Compare pooled INCOME estimate to OLS complete-case and to Heckman.

source(here::here("R", "utils.R"))
load_libs()
# library(mice)
# library(miceMNAR)

stop("TODO: implement after 03_gjrm_copula.R validates the selection structure.")
