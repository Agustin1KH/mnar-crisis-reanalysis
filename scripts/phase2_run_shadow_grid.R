# scripts/phase2_run_shadow_grid.R
# Phase-2 driver: run the FAST stage (scripts 02 + 03b) for each of the
# three shadow variables, then stitch the per-shadow CSVs into the two
# headline cross-shadow tables and validate against the Python Phase-1
# predictions.
#
# Outputs (always written, regardless of validation outcome):
#   output/tables/phase2_exclusion_across_shadows.md
#       8 outcomes x 3 shadows x [coef, SE, |t|, p, verdict]; the
#       headline grid asked for in cursor-phase2-task-spec.md Task 1.
#   output/tables/phase2_gap_sum_heckman_across_shadows.md
#       gap_sum Spec A / Spec B Heckman ML / 2-step INCOME coefficients
#       across the three shadows, plus complete-case OLS reference.
#   output/tables/phase2_checkpoint_validation.md
#       Pass/fail report against the Python predictions, with explicit
#       per-cell tolerances (Heckman 2-step coefs +/- 0.01; 03b |t|
#       stats +/- 0.1).
#
# Side effects: invokes scripts/02_heckman_baseline.R and
# scripts/03b_exclusion_insample_test.R via Rscript with SHADOW_VAR set
# in the environment, three times each.
#
# This is the FAST path (~2-3 min wall-clock total). The slow stage
# (scripts/03_gjrm_copula.R, 03a_exclusion_diagnostic.R, 04_mice_mnar.R,
# 05_sensitivity.R) is gated by the checkpoint validation outcome and is
# launched separately by scripts/phase2_run_full_pipeline.R.

source(here::here("R", "utils.R"))
load_libs()

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(cli)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

SHADOWS <- c("n_chron_clean", "n_chron_3chron", "n_chron_3chron_raw")

# Python Phase-1 predictions (from files (1)/phase1_insample_test.csv and
# the Heckman 2-step run in files (1)/03_heckman_2step.py). Tolerances
# per the user's checkpoint instruction:
#   Heckman 2-step INCOME coefs: |delta| <= 0.01
#   03b |t| stats:               |delta| <= 0.10
PYTHON_HECKMAN_2STEP_GAPA_INCOME <- c(
  n_chron_clean      = -0.2364,
  n_chron_3chron     = -0.2521,
  n_chron_3chron_raw = -0.2638
)
TOL_HECKMAN_INCOME <- 0.01

# 03b |t| stats from files (1)/phase1_insample_test.csv. Read on demand
# below so we don't have to copy them by hand.

# =============================================================================
# Step 1: invoke 02 + 03b for each shadow in clean R sessions
# =============================================================================
run_shadow <- function(shadow) {
  cli::cli_h1("=== SHADOW_VAR = {shadow} ===")
  # Set SHADOW_VAR in this process; spawned Rscripts inherit the env. We
  # reset to empty after the loop so the driver leaves no stale env behind.
  Sys.setenv(SHADOW_VAR = shadow)
  on.exit(Sys.unsetenv("SHADOW_VAR"), add = TRUE)
  for (script in c("scripts/02_heckman_baseline.R",
                   "scripts/03b_exclusion_insample_test.R")) {
    cli::cli_alert_info("Invoking {.path {script}} (SHADOW_VAR={shadow}) ...")
    res <- system2(
      command = "Rscript",
      args    = shQuote(here::here(script)),
      stdout  = TRUE, stderr = TRUE
    )
    status <- attr(res, "status") %||% 0L
    if (!isTRUE(status == 0L)) {
      cat(paste(res, collapse = "\n"), "\n")
      stop(sprintf("Script %s exited non-zero (%s) under SHADOW_VAR=%s",
                   script, status, shadow), call. = FALSE)
    }
    # Print the last 6 lines (success markers) for visibility.
    tail_lines <- tail(res, 6L)
    cat(paste(tail_lines, collapse = "\n"), "\n")
  }
}

for (sv in SHADOWS) {
  run_shadow(sv)
}

# =============================================================================
# Step 2: collect the per-shadow 03b CSVs into one tidy table
# =============================================================================
collect_03b <- function(shadow) {
  sub <- if (shadow == "n_chron_clean") "" else paste0("phase2_shadow_", shadow)
  fp  <- phase2_output_path("03b_results.csv", sub, ensure_dir = FALSE)
  if (!file.exists(fp)) {
    stop(sprintf("Expected 03b CSV missing for shadow=%s: %s", shadow, fp),
         call. = FALSE)
  }
  read.csv(fp, stringsAsFactors = FALSE)
}

grid_long <- dplyr::bind_rows(lapply(SHADOWS, collect_03b))
# Restrict to the SHADOW_VAR rows (each per-shadow CSV also includes the
# placebo `maddison_priority`; we keep it but flag it separately).
grid_focal   <- grid_long |> dplyr::filter(focal == shadow)
grid_placebo <- grid_long |> dplyr::filter(focal == "maddison_priority")

# =============================================================================
# Step 3: assemble the headline cross-shadow exclusion-grid table
# =============================================================================
fmt_n <- function(x, d = 4) {
  if (length(x) == 0L) return("—")
  vapply(x, function(xi) {
    if (is.na(xi)) "—" else formatC(xi, format = "f", digits = d)
  }, character(1))
}
fmt_p <- function(x) {
  if (length(x) == 0L) return("—")
  vapply(x, function(xi) {
    if (is.na(xi))     "—"
    else if (xi < 1e-4) "<0.0001"
    else                formatC(xi, format = "f", digits = 4)
  }, character(1))
}

OUTCOME_ORDER <- c(
  "gap_sum",              "guarantees_d",       "lending_d",
  "capital_injections_d", "restructuring_d",    "asset_management_d",
  "rules_d",              "other_d"
)

# Build the wide cross-shadow table per outcome:
#   one row per outcome
#   columns grouped by shadow: coef | SE | |t| | p | verdict
grid_focal <- grid_focal |>
  dplyr::mutate(
    outcome = factor(outcome, levels = OUTCOME_ORDER),
    shadow  = factor(shadow,  levels = SHADOWS)
  ) |>
  dplyr::arrange(outcome, shadow)

shadow_label <- function(s) switch(s,
  n_chron_clean      = "n_chron_clean (4-chron, paper)",
  n_chron_3chron     = "n_chron_3chron (3-chron, fallback)",
  n_chron_3chron_raw = "n_chron_3chron_raw (3-chron, no fallback)",
  s)

# Header
md <- c(
  "# Phase-2 cross-shadow exclusion-restriction in-sample grid",
  "",
  "Generated by `scripts/phase2_run_shadow_grid.R`. Headline deliverable for ",
  "Task 1 in `files (1)/cursor-phase2-task-spec.md`. The full 8-outcome x ",
  "3-shadow grid for the in-sample partial-correlation test of the active ",
  "exclusion variable in the Heckman selection equation, against the ",
  "complete-case (year >= 1800, r == 1) sample.",
  "",
  "Per outcome and shadow, fits ",
  "",
  "```",
  "<outcome> ~ candidate + log_lagged_gdppc + polity + currency_regime +",
  "            debt + exp + year + <shadow>",
  "```",
  "",
  "with `lm()` for the continuous outcome (`gap_sum`) and ",
  "`glm(family = binomial(\"logit\"))` for the seven binary outcomes. ",
  "Verdict thresholds (per Prompt 3e spec):",
  "",
  "- `|t| < 1.5` AND `p > 0.10` -> **defensible**",
  "- `1.5 <= |t| < 2.0` -> **weak violation**",
  "- `|t| >= 2.0` -> **in-sample violation**",
  "",
  "**Caveat on `other_d`.** EPV ~ 1.7 in complete-case (Peduzzi et al. ",
  "1996); the standalone partial-correlation t-stat is itself noisy. ",
  "Included for grid completeness; do not weight in the synthesis.",
  "",
  "## Per-outcome x shadow grid",
  "",
  paste(c("| outcome",
          paste0(SHADOWS, "<br/>coef"),
          paste0(SHADOWS, "<br/>|t|"),
          paste0(SHADOWS, "<br/>verdict")),
        collapse = " | ") |> paste0(" |"),
  paste(c("| :---", rep("---:", 2 * length(SHADOWS)),
          rep(":---:", length(SHADOWS))),
        collapse = " | ") |> paste0(" |")
)

# Render: one row per outcome, with all 3 shadows side by side.
for (oc in OUTCOME_ORDER) {
  rows <- grid_focal |> dplyr::filter(outcome == oc)
  coef_cells <- vapply(SHADOWS, function(s) {
    r <- rows[rows$shadow == s, ]
    if (nrow(r) == 0L) "—" else fmt_n(r$coef[1], 4)
  }, character(1))
  t_cells <- vapply(SHADOWS, function(s) {
    r <- rows[rows$shadow == s, ]
    if (nrow(r) == 0L) "—" else fmt_n(r$abs_test[1], 2)
  }, character(1))
  verdict_cells <- vapply(SHADOWS, function(s) {
    r <- rows[rows$shadow == s, ]
    if (nrow(r) == 0L) "—" else r$verdict_short[1]
  }, character(1))
  md <- c(md, paste(c(sprintf("| `%s`", oc),
                      coef_cells, t_cells, verdict_cells),
                    collapse = " | ") |> paste0(" |"))
}

# Per-shadow violation summary
md <- c(md, "",
  "## Verdict summary by shadow",
  "",
  "| shadow | defensible | weak | violation |",
  "| :--- | ---: | ---: | ---: |"
)
for (s in SHADOWS) {
  rows <- grid_focal |> dplyr::filter(shadow == s)
  md <- c(md, sprintf(
    "| %s | %d | %d | %d |",
    shadow_label(s),
    sum(rows$verdict_short == "defensible"),
    sum(rows$verdict_short == "weak violation"),
    sum(rows$verdict_short == "in-sample violation")
  ))
}

# Placebo cross-shadow table
md <- c(md, "",
  "## Placebo (`maddison_priority`) cross-shadow grid",
  "",
  "Provides a calibration baseline: `maddison_priority` is in the ",
  "selection equation but is not claimed as the exclusion restriction. ",
  "The placebo's |t| stats should be invariant to the choice of ",
  "SHADOW_VAR (the partial-correlation regression has the same ",
  "covariates and the same complete-case sample across shadows; only the ",
  "placebo focal stays the same). Any non-trivial drift across the three ",
  "columns would indicate a sample-composition issue.",
  "",
  paste(c("| outcome",
          paste0(SHADOWS, "<br/>coef"),
          paste0(SHADOWS, "<br/>|t|")),
        collapse = " | ") |> paste0(" |"),
  paste(c("| :---", rep("---:", 2 * length(SHADOWS))),
        collapse = " | ") |> paste0(" |")
)
for (oc in OUTCOME_ORDER) {
  rows <- grid_placebo |> dplyr::filter(outcome == oc)
  coef_cells <- vapply(SHADOWS, function(s) {
    r <- rows[rows$shadow == s, ]
    if (nrow(r) == 0L) "—" else fmt_n(r$coef[1], 4)
  }, character(1))
  t_cells <- vapply(SHADOWS, function(s) {
    r <- rows[rows$shadow == s, ]
    if (nrow(r) == 0L) "—" else fmt_n(r$abs_test[1], 2)
  }, character(1))
  md <- c(md, paste(c(sprintf("| `%s`", oc), coef_cells, t_cells),
                    collapse = " | ") |> paste0(" |"))
}

writeLines(md, here::here("output", "tables",
                          "phase2_exclusion_across_shadows.md"))
cli::cli_alert_success("Wrote {.path output/tables/phase2_exclusion_across_shadows.md}")

# =============================================================================
# Step 4: assemble the gap_sum Heckman cross-shadow comparison
# =============================================================================
collect_gap <- function(shadow) {
  sub <- if (shadow == "n_chron_clean") "" else paste0("phase2_shadow_", shadow)
  fp  <- phase2_output_path("02_gap_sum_heckman_summary.csv", sub,
                             ensure_dir = FALSE)
  if (!file.exists(fp)) {
    stop(sprintf("Expected 02 gap-sum CSV missing for shadow=%s: %s",
                 shadow, fp), call. = FALSE)
  }
  read.csv(fp, stringsAsFactors = FALSE)
}

gap_long <- dplyr::bind_rows(lapply(SHADOWS, collect_gap))
gap_long$shadow <- factor(gap_long$shadow, levels = SHADOWS)
gap_long$fit_label <- factor(gap_long$fit_label,
                              levels = c("gap_a_2step", "gap_a_ml",
                                         "gap_b_2step", "gap_b_ml"))

gap_md <- c(
  "# Phase-2 cross-shadow Heckman comparison: `gap_sum` Spec A and B INCOME",
  "",
  "Generated by `scripts/phase2_run_shadow_grid.R`. Companion to ",
  "`phase2_exclusion_across_shadows.md`. Reports the Heckman ML and 2-step ",
  "INCOME (`log_lagged_gdppc`) coefficient on `gap_sum` (Spec A and B) ",
  "across the three Phase-2 shadow variables, plus the complete-case OLS ",
  "(r == 1) reference for each spec.",
  "",
  "Per fit: `Heckman ML` (sampleSelection joint MLE) vs `Heckman 2-step` ",
  "(Greene/Heckman 1979 two-step). The 2-step is the calibration target ",
  "for the Phase-1 Python predictions in ",
  "`files (1)/cursor-phase2-task-spec.md`.",
  "",
  "## Spec A INCOME comparison (gap_sum ~ candidate + INCOME + polity + ",
  "currency_regime + year)",
  "",
  "| shadow | OLS(cc) | OLS N | Heckman 2-step | Heckman ML | rho ML | |t_rho| ML |",
  "| :--- | ---: | ---: | ---: | ---: | ---: | ---: |"
)
for (s in SHADOWS) {
  ols_row  <- gap_long |> dplyr::filter(shadow == s, fit_label == "gap_a_ml")
  ml_row   <- ols_row
  twos_row <- gap_long |> dplyr::filter(shadow == s, fit_label == "gap_a_2step")
  gap_md <- c(gap_md, sprintf(
    "| %s | %s | %d | %s | %s | %s | %s |",
    shadow_label(s),
    fmt_n(ols_row$income_ols_cc[1]),
    ols_row$n_ols_cc[1],
    fmt_n(twos_row$income_heck[1]),
    fmt_n(ml_row$income_heck[1]),
    fmt_n(ml_row$rho[1]),
    fmt_n(abs(ml_row$rho_t[1] %||% NA), 2)
  ))
}

gap_md <- c(gap_md, "",
  "## Spec B INCOME comparison (Spec A + debt + exp)",
  "",
  "| shadow | OLS(cc) | OLS N | Heckman 2-step | Heckman ML | rho ML | |t_rho| ML |",
  "| :--- | ---: | ---: | ---: | ---: | ---: | ---: |"
)
for (s in SHADOWS) {
  ols_row  <- gap_long |> dplyr::filter(shadow == s, fit_label == "gap_b_ml")
  ml_row   <- ols_row
  twos_row <- gap_long |> dplyr::filter(shadow == s, fit_label == "gap_b_2step")
  gap_md <- c(gap_md, sprintf(
    "| %s | %s | %d | %s | %s | %s | %s |",
    shadow_label(s),
    fmt_n(ols_row$income_ols_cc[1]),
    ols_row$n_ols_cc[1],
    fmt_n(twos_row$income_heck[1]),
    fmt_n(ml_row$income_heck[1]),
    fmt_n(ml_row$rho[1]),
    fmt_n(abs(ml_row$rho_t[1] %||% NA), 2)
  ))
}

writeLines(gap_md, here::here("output", "tables",
                               "phase2_gap_sum_heckman_across_shadows.md"))
cli::cli_alert_success(
  "Wrote {.path output/tables/phase2_gap_sum_heckman_across_shadows.md}")

# =============================================================================
# Step 5: validate against Python predictions; emit checkpoint report
# =============================================================================
cli::cli_h1("Cross-shadow checkpoint validation")

# 5a: Heckman 2-step gap_sum Spec A INCOME vs Python predictions
hi_check <- vapply(SHADOWS, function(s) {
  twos_row <- gap_long |> dplyr::filter(shadow == s, fit_label == "gap_a_2step")
  r_val <- twos_row$income_heck[1]
  py_val <- PYTHON_HECKMAN_2STEP_GAPA_INCOME[[s]]
  abs(r_val - py_val)
}, numeric(1))
hi_pass <- hi_check <= TOL_HECKMAN_INCOME
hi_overall <- all(hi_pass)

# 5b: 03b |t| stats vs Python predictions (read directly from
# files (1)/phase1_insample_test.csv).
py_csv_path <- here::here("files (1)", "phase1_insample_test.csv")
py_03b_pass <- NA
py_03b_table <- NULL
if (file.exists(py_csv_path)) {
  py <- read.csv(py_csv_path, stringsAsFactors = FALSE,
                 check.names = FALSE)
  py <- py[py$focal %in% SHADOWS, c("outcome", "focal", "|t/z|", "verdict")]
  names(py)[3:4] <- c("py_abs_t", "py_verdict")
  r_focal <- grid_focal[, c("outcome", "shadow", "abs_test",
                            "verdict_short")]
  names(r_focal)[2:4] <- c("focal", "r_abs_t", "r_verdict")
  py_03b_table <- merge(py, r_focal, by = c("outcome", "focal"))
  py_03b_table$delta_abs_t <- py_03b_table$r_abs_t - py_03b_table$py_abs_t
  py_03b_table$within_tol <- abs(py_03b_table$delta_abs_t) <= 0.10
  py_03b_pass <- all(py_03b_table$within_tol, na.rm = TRUE)
}

# Emit the checkpoint markdown.
cp_md <- c(
  "# Phase-2 checkpoint validation",
  "",
  "Generated by `scripts/phase2_run_shadow_grid.R`. Validates the R ",
  "outputs against the Python Phase-1 predictions in `files (1)/`. ",
  "Two tolerances per the user's checkpoint instruction:",
  "",
  "- Heckman 2-step `gap_sum` Spec A INCOME coefficient: `|delta| <= 0.01`",
  "- 03b in-sample `|t|` stats (8 outcomes x 3 shadows): `|delta| <= 0.10`",
  "",
  "## Heckman 2-step gap_sum Spec A INCOME",
  "",
  "| shadow | R | Python (predicted) | |delta| | within tol? |",
  "| :--- | ---: | ---: | ---: | :---: |"
)
for (s in SHADOWS) {
  twos_row <- gap_long |> dplyr::filter(shadow == s, fit_label == "gap_a_2step")
  r_val <- twos_row$income_heck[1]
  py_val <- PYTHON_HECKMAN_2STEP_GAPA_INCOME[[s]]
  cp_md <- c(cp_md, sprintf(
    "| %s | %s | %s | %s | %s |",
    s, fmt_n(r_val, 4), fmt_n(py_val, 4),
    fmt_n(abs(r_val - py_val), 4),
    if (abs(r_val - py_val) <= TOL_HECKMAN_INCOME) "**yes**" else "**NO**"
  ))
}
cp_md <- c(cp_md, "",
  sprintf("Heckman 2-step gap_sum Spec A overall: **%s** (max |delta| = %s).",
          if (hi_overall) "PASS" else "FAIL",
          fmt_n(max(hi_check), 4)),
  ""
)

if (!is.null(py_03b_table)) {
  cp_md <- c(cp_md,
    "## 03b in-sample exclusion test |t| stats (full grid)",
    "",
    "| outcome | shadow | R |t| | Python |t| | |delta| | within tol? |",
    "| :--- | :--- | ---: | ---: | ---: | :---: |"
  )
  py_03b_sorted <- py_03b_table[order(
    factor(py_03b_table$outcome, levels = OUTCOME_ORDER),
    factor(py_03b_table$focal, levels = SHADOWS)), ]
  for (i in seq_len(nrow(py_03b_sorted))) {
    r <- py_03b_sorted[i, ]
    cp_md <- c(cp_md, sprintf(
      "| `%s` | %s | %s | %s | %s | %s |",
      r$outcome, r$focal,
      fmt_n(r$r_abs_t, 2), fmt_n(r$py_abs_t, 2),
      fmt_n(abs(r$delta_abs_t), 2),
      if (isTRUE(r$within_tol)) "**yes**" else "**NO**"
    ))
  }
  cp_md <- c(cp_md, "",
    sprintf("03b |t| overall: **%s** (max |delta| = %s).",
            if (isTRUE(py_03b_pass)) "PASS" else "FAIL",
            fmt_n(max(abs(py_03b_table$delta_abs_t), na.rm = TRUE), 4)),
    ""
  )
} else {
  cp_md <- c(cp_md,
    sprintf("**Python 03b CSV not found at %s; skipping 03b validation.**",
            py_csv_path),
    "")
}

overall_pass <- hi_overall && (is.na(py_03b_pass) || isTRUE(py_03b_pass))
cp_md <- c(cp_md,
  "## Overall verdict",
  "",
  if (overall_pass) {
    "**CHECKPOINT PASSED.** R reproduces the Python Phase-1 numbers within tolerance on both axes; safe to proceed to the slow stage (03/03a/04/05) for the two non-default shadows."
  } else {
    "**CHECKPOINT FAILED.** R does NOT reproduce the Python predictions within tolerance. STOP and investigate the per-row deltas above before continuing."
  },
  ""
)

writeLines(cp_md, here::here("output", "tables",
                              "phase2_checkpoint_validation.md"))
cli::cli_alert_success(
  "Wrote {.path output/tables/phase2_checkpoint_validation.md}")

cli::cli_h1("Phase-2 fast stage complete")
cli::cli_alert_info("Heckman 2-step Spec A: max |delta| = {fmt_n(max(hi_check), 4)} \\
                    ({if (hi_overall) 'PASS' else 'FAIL'}).")
if (!is.na(py_03b_pass)) {
  cli::cli_alert_info("03b |t| stats: max |delta| = \\
                      {fmt_n(max(abs(py_03b_table$delta_abs_t), na.rm = TRUE), 4)} \\
                      ({if (isTRUE(py_03b_pass)) 'PASS' else 'FAIL'}).")
}
if (overall_pass) {
  cli::cli_alert_success("OVERALL: CHECKPOINT PASSED.")
} else {
  cli::cli_alert_danger("OVERALL: CHECKPOINT FAILED.")
}
