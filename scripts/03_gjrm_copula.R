# scripts/03_gjrm_copula.R
# Stress-test the Heckman baseline by relaxing joint normality of errors.
# GJRM (Marra & Radice) fits the same selection structure but with a
# copula between outcome and selection errors. If the sign / magnitude of
# INCOME is stable across N, Clayton, Joe, Gumbel, Frank, and t copulas,
# the baseline result isn't an artefact of the normality assumption.
#
# TODO once the Heckman baseline results look sensible:
#   1. Fit GJRM Model = "BSS" (Bivariate with Sample Selection) for each copula.
#   2. Compare INCOME estimate across copulas.
#   3. Report AIC / BIC of copulas to pick the data-preferred one.
#   4. Compare its INCOME estimate to the paper's OLS.

source(here::here("R", "utils.R"))
load_libs()
# library(GJRM)

stop("TODO: implement after 02_heckman_baseline.R produces sensible results.")
