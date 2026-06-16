# scripts/phase2_run_full_pipeline.R
# Phase-2 SLOW stage: runs the heavy-compute scripts (03 -> 03a -> 04 -> 05)
# under each of the three shadow values, after the FAST stage's checkpoint
# has passed (see scripts/phase2_run_shadow_grid.R).
#
# Assumes:
#   - scripts/00_etl.R has been run (analysis_data.csv has the new columns)
#   - scripts/phase2_run_shadow_grid.R has been run (heckman_fits.rds is
#     in the namespaced cache for each of the three shadows)
#   - output/tables/phase2_checkpoint_validation.md says CHECKPOINT PASSED
#     (this driver enforces the precondition by reading that file and
#     halting if the verdict line is absent or contradictory)
#
# Wall-clock budget: ~60-100 min total; MICE-MNAR is the slow step
# (~20 min per shadow at m=20). GJRM bootstrap (~10 min per shadow);
# 03a (~30 s); 05 (~5 min per shadow).
#
# Cache reuse: if a previous run completed for a particular shadow x
# script combination, the cached RDS in
# data/processed/phase2_shadow_<sv>/ will be reused (e.g. mice_mnar.rds
# is the largest single cost; once it exists it's not re-fitted). Force
# a recomputation by deleting the cached RDS for that shadow.

source(here::here("R", "utils.R"))
load_libs()

suppressPackageStartupMessages({
  library(cli)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

SHADOWS <- c("n_chron_clean", "n_chron_3chron", "n_chron_3chron_raw")
SCRIPTS <- c(
  "scripts/03_gjrm_copula.R",
  "scripts/03a_exclusion_diagnostic.R",
  "scripts/04_mice_mnar.R",
  "scripts/05_sensitivity.R"
)

# Precondition: checkpoint must have passed.
cp_path <- here::here("output", "tables",
                      "phase2_checkpoint_validation.md")
if (!file.exists(cp_path)) {
  stop("Pre-condition unmet: checkpoint file missing. Run ",
       "scripts/phase2_run_shadow_grid.R first.", call. = FALSE)
}
cp_lines <- readLines(cp_path)
if (!any(grepl("CHECKPOINT PASSED", cp_lines, fixed = TRUE))) {
  stop("Pre-condition unmet: checkpoint file does not say CHECKPOINT PASSED. ",
       "Investigate ", cp_path, " before running the slow stage.",
       call. = FALSE)
}
cli::cli_alert_success("Checkpoint precondition satisfied.")

# Run loop.
run_script <- function(shadow, script) {
  cli::cli_h2("[{shadow}] {script}")
  Sys.setenv(SHADOW_VAR = shadow)
  on.exit(Sys.unsetenv("SHADOW_VAR"), add = TRUE)
  t0 <- Sys.time()
  res <- system2(
    command = "Rscript",
    args    = shQuote(here::here(script)),
    stdout  = TRUE, stderr = TRUE
  )
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- attr(res, "status") %||% 0L
  if (!isTRUE(status == 0L)) {
    cat(paste(res, collapse = "\n"), "\n")
    stop(sprintf("Script %s failed under SHADOW_VAR=%s (exit=%s, %.1f s).",
                 script, shadow, status, dt), call. = FALSE)
  }
  cli::cli_alert_success("  done in {round(dt, 1)} s.")
  cat(paste(tail(res, 5L), collapse = "\n"), "\n")
}

cli::cli_h1("Phase-2 SLOW stage: 4 scripts x 3 shadows")
total_t0 <- Sys.time()

for (shadow in SHADOWS) {
  cli::cli_h1("=== SHADOW_VAR = {shadow} ===")
  for (script in SCRIPTS) {
    run_script(shadow, script)
  }
}

total_dt <- as.numeric(difftime(Sys.time(), total_t0, units = "secs"))
cli::cli_h1("Phase-2 SLOW stage complete in {round(total_dt / 60, 1)} min.")
