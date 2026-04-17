# scripts/05_sensitivity.R
# Sensitivity analysis: the MNAR result is only as credible as the
# identifying assumptions. Report how INCOME moves across two dimensions:
#
#   1. rho grid — fix the selection/outcome error correlation at
#      {-0.7, -0.4, -0.2, 0, 0.2, 0.4, 0.7} instead of estimating it,
#      refit, and report how INCOME shifts. If the qualitative conclusion
#      holds across reasonable rho, the result is robust.
#
#   2. Exclusion-restriction swap — rerun 02_heckman_baseline.R with:
#      (a) n_chron_clean only
#      (b) maddison_priority only
#      (c) both
#      (d) neither (identification via normality only — fragile, cited as
#          evidence of what Heckit loses without an instrument)
#
#   3. Manski-Horowitz worst-case bounds on the population INCOME effect
#      under unrestricted MNAR — as a floor. If parametric posterior sits
#      well inside these bounds, parametric assumptions are doing
#      substantial work.
#
# TODO after 02 runs.

source(here::here("R", "utils.R"))
load_libs()

stop("TODO: implement after 02_heckman_baseline.R produces stable estimates.")
