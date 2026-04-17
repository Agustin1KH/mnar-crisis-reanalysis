# scripts/03b_exclusion_insample_test.R
# In-sample partial-correlation test of the Heckman exclusion restriction
# (Prompt 3e). Companion diagnostic to 03a_exclusion_diagnostic.R.
#
# Question: conditional on the full paper specification (candidate,
# log_lagged_gdppc, polity, currency_regime, debt, exp, year), does
# `n_chron_clean` -- the shadow / exclusion variable in the Heckman
# selection equation -- have any RESIDUAL partial correlation with the
# outcome in the complete-case (r == 1) sample?
#
# If `n_chron_clean` has no observable-channel power conditional on the
# regressors, that strengthens the exclusion claim: the "shadow variable"
# is not sneaking in via a measurable backdoor on the outcome equation.
# If it does, the exclusion is in-sample-falsified and the memo framing
# needs to acknowledge a specific weakness.
#
# Three headline outcomes (matching the original paper's Table 3 / Table 2
# headline interventions):
#   - gap_sum        : continuous, lm(...)
#   - guarantees_d   : binary, glm(... , family = binomial("logit"))
#   - lending_d      : binary, glm(... , family = binomial("logit"))
#
# Per outcome, fit
#     <outcome> ~ candidate + log_lagged_gdppc + polity + currency_regime +
#                 debt + exp + year + <variable of interest>
# with <variable of interest> = n_chron_clean (the actual claim) and then
# = maddison_priority (placebo / calibration).
#
# Placebo logic: maddison_priority is in the paper's selection equation
# but is NOT claimed as an exclusion restriction (it's also a covariate
# implicitly via the Schmelzing / Maddison data-availability story). If
# it shows similar in-sample partial correlation to n_chron_clean, that
# pattern is consistent with noise at this sample size; if it's clean
# and n_chron_clean is not, then n_chron_clean has a specific weakness.
#
# Verdict thresholds (per Prompt 3e spec):
#   |t| < 1.5 AND p > 0.1   --> "defensible in-sample"
#   1.5 <= |t| < 2.0        --> "weak violation -- flag for memo"
#   |t| >= 2.0              --> "in-sample violation -- exclusion
#                                restriction undermined"
#
# Output: output/tables/03b_insample_exclusion_test.md
# Constraints: no modifications to analysis_data.csv or to any
# transformation script; no refits of prior-prompt models; one page max.

source(here::here("R", "utils.R"))
load_libs()

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(cli)
})

# -- data --------------------------------------------------------------------
df <- read.csv(here::here("data", "processed", "analysis_data.csv"),
               stringsAsFactors = FALSE)

df_complete <- df |>
  dplyr::filter(year >= 1800, r == 1L)

cli::cli_h1("In-sample exclusion-restriction partial-correlation test")
cli::cli_alert_info(
  "Complete-case sample (year >= 1800, r == 1): N = {nrow(df_complete)}.")

# -- spec -------------------------------------------------------------------
# Three outcomes. lm for the continuous one, logit glm for the two binary
# ones. Two "variables of interest": the actual exclusion restriction
# (n_chron_clean) and the placebo / calibration covariate
# (maddison_priority).
outcomes <- list(
  list(name = "gap_sum",      kind = "continuous"),
  list(name = "guarantees_d", kind = "binary"),
  list(name = "lending_d",    kind = "binary")
)
focal_vars <- c("n_chron_clean", "maddison_priority")

base_rhs <- "candidate + log_lagged_gdppc + polity + currency_regime + debt + exp + year"

# -- fit + extract focal-variable row ---------------------------------------
fit_and_extract <- function(outcome_name, kind, focal_var, dat) {
  rhs <- paste(base_rhs, "+", focal_var)
  fml <- as.formula(paste(outcome_name, "~", rhs))
  fit <- if (kind == "continuous") {
    stats::lm(fml, data = dat)
  } else {
    stats::glm(fml, data = dat, family = stats::binomial("logit"))
  }
  cf <- summary(fit)$coefficients
  if (!focal_var %in% rownames(cf)) {
    stop(sprintf("Focal var %s not in coefficient matrix for %s.",
                 focal_var, outcome_name), call. = FALSE)
  }
  row <- cf[focal_var, , drop = TRUE]
  # lm gives "t value" (col 3) and "Pr(>|t|)" (col 4); glm gives
  # "z value" (col 3) and "Pr(>|z|)" (col 4). Either way the third entry
  # is the test statistic and the fourth is the two-sided p-value.
  list(
    outcome  = outcome_name,
    kind     = kind,
    focal    = focal_var,
    coef     = unname(row[[1]]),
    se       = unname(row[[2]]),
    test     = unname(row[[3]]),
    abs_test = abs(unname(row[[3]])),
    p_value  = unname(row[[4]]),
    n_used   = nobs(fit)
  )
}

verdict_for <- function(focal, abs_test, p_value, outcome) {
  if (is.na(abs_test) || is.na(p_value)) {
    return("undefined (model failed to identify the focal coefficient)")
  }
  if (abs_test < 1.5 && p_value > 0.10) {
    sprintf("exclusion restriction defensible in-sample on %s", outcome)
  } else if (abs_test < 2.0) {
    sprintf("weak violation -- %s has marginal in-sample partial correlation; flag for memo",
            focal)
  } else {
    sprintf("in-sample violation -- %s has a direct partial channel to %s; exclusion restriction is undermined",
            focal, outcome)
  }
}

verdict_short <- function(abs_test, p_value) {
  if (is.na(abs_test) || is.na(p_value)) return("undefined")
  if (abs_test < 1.5 && p_value > 0.10) return("defensible")
  if (abs_test < 2.0)                   return("weak violation")
  "in-sample violation"
}

# -- run all 6 fits ---------------------------------------------------------
rows <- list()
for (out in outcomes) {
  for (fv in focal_vars) {
    r <- fit_and_extract(out$name, out$kind, fv, df_complete)
    r$verdict       <- verdict_for(fv, r$abs_test, r$p_value, out$name)
    r$verdict_short <- verdict_short(r$abs_test, r$p_value)
    rows[[length(rows) + 1L]] <- r
    cli::cli_alert_info(
      "{out$name} (kind={out$kind}) ~ ... + {fv}: \\
       coef={sprintf('%.4f', r$coef)}, |t|={sprintf('%.2f', r$abs_test)}, \\
       p={sprintf('%.4f', r$p_value)}, n={r$n_used} -> {r$verdict_short}")
  }
}
res_tbl <- dplyr::bind_rows(lapply(rows, tibble::as_tibble))

# -- splits for the two interpretation paragraphs ---------------------------
res_nchron  <- res_tbl |> dplyr::filter(focal == "n_chron_clean")
res_placebo <- res_tbl |> dplyr::filter(focal == "maddison_priority")

# -- render markdown --------------------------------------------------------
fmt_n <- function(x, d = 4) {
  if (is.null(x)) return("—")
  vapply(x, function(xi) {
    if (is.na(xi)) "—" else formatC(xi, format = "f", digits = d)
  }, character(1))
}
fmt_p <- function(x) {
  if (is.null(x)) return("—")
  vapply(x, function(xi) {
    if (is.na(xi))     "—"
    else if (xi < 1e-4) "<0.0001"
    else                formatC(xi, format = "f", digits = 4)
  }, character(1))
}

lines <- c(
  "# In-sample exclusion-restriction partial-correlation test",
  "",
  "Generated by `scripts/03b_exclusion_insample_test.R`. Companion to ",
  "`output/tables/03a_exclusion_diagnostic.md` (Prompt 3d) and the GJRM ",
  "cross-package check in `output/tables/03_interpretation.md` (Phase 3a). ",
  "Question per the Perplexity review: within the complete-case sample ",
  sprintf("(year >= 1800, r == 1, N = %d), does `n_chron_clean` have any ",
          nrow(df_complete)),
  "residual partial correlation with each of the ",
  "three headline outcomes (`gap_sum`, `guarantees_d`, `lending_d`) once the ",
  "full paper specification is conditioned on? If not, the exclusion ",
  "restriction is defensible in observables; if so, it is in-sample-falsified.",
  "",
  "Specification (per outcome): ",
  "",
  "```",
  "<outcome> ~ candidate + log_lagged_gdppc + polity + currency_regime +",
  "            debt + exp + year + <focal variable>",
  "```",
  "",
  "Focal variables: `n_chron_clean` (the actual Heckman exclusion ",
  "restriction) and `maddison_priority` (placebo / calibration -- a ",
  "covariate in the selection equation but not claimed as an exclusion ",
  "restriction). Continuous outcome uses `lm()`; binary outcomes use ",
  "`glm(family = binomial(\"logit\"))`. Verdict thresholds: ",
  "`|t| < 1.5` AND `p > 0.10` -> defensible; `1.5 <= |t| < 2.0` -> weak ",
  "violation; `|t| >= 2.0` -> in-sample violation.",
  "",
  "## Summary table",
  "",
  "| outcome | variable | coef | SE | \\|t\\| | p-value | N | verdict |",
  "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | :--- |"
)
for (i in seq_len(nrow(res_tbl))) {
  r <- res_tbl[i, ]
  lines <- c(lines, sprintf(
    "| `%s` | `%s` | %s | %s | %s | %s | %d | %s |",
    r$outcome, r$focal,
    fmt_n(r$coef, 4), fmt_n(r$se, 4),
    fmt_n(r$abs_test, 2), fmt_p(r$p_value),
    r$n_used, r$verdict_short
  ))
}

# -- two-paragraph interpretation -------------------------------------------
all_pass_nchron <- all(res_nchron$verdict_short == "defensible")
fail_outcomes_nchron <- res_nchron$outcome[res_nchron$verdict_short != "defensible"]
weak_n  <- sum(res_nchron$verdict_short == "weak violation")
viol_n  <- sum(res_nchron$verdict_short == "in-sample violation")

placebo_any_nondef <- any(res_placebo$verdict_short != "defensible")
placebo_fails  <- res_placebo$outcome[res_placebo$verdict_short != "defensible"]
nchron_fails   <- res_nchron$outcome[res_nchron$verdict_short != "defensible"]
shared_fails   <- intersect(placebo_fails, nchron_fails)
nchron_unique  <- setdiff(nchron_fails, placebo_fails)

# Per-outcome summary string for the first paragraph
fmt_per_outcome <- function(r) {
  sprintf("`%s` (|t| = %s, p = %s)",
          r$outcome, fmt_n(r$abs_test, 2), fmt_p(r$p_value))
}
nchron_lines <- vapply(seq_len(nrow(res_nchron)),
                       function(i) fmt_per_outcome(res_nchron[i, ]),
                       character(1))

para1 <- if (all_pass_nchron) {
  sprintf(paste0(
    "**The exclusion restriction is defensible in observables across all ",
    "three headline outcomes.** Conditional on the full Table 3 / Table 2 ",
    "specification (`candidate + log_lagged_gdppc + polity + ",
    "currency_regime + debt + exp + year`), `n_chron_clean` shows no ",
    "significant partial correlation with any outcome: %s. All three ",
    "satisfy `|t| < 1.5` AND `p > 0.10`, the strict-reading threshold for ",
    "in-sample defensibility."),
    paste(nchron_lines, collapse = "; "))
} else if (viol_n > 0L) {
  sprintf(paste0(
    "**The exclusion restriction is in-sample-violated on at least one ",
    "outcome.** %s. Conditional on the full paper specification, ",
    "`n_chron_clean` has a partial t-statistic of magnitude >= 2.0 on the ",
    "named outcome(s), which falsifies the literal-reading exclusion ",
    "claim in observables. The remaining outcome diagnostics: %s."),
    paste(sprintf("`%s` exceeds the |t| = 2.0 threshold (|t| = %s, p = %s)",
                  res_nchron$outcome[res_nchron$verdict_short == "in-sample violation"],
                  fmt_n(res_nchron$abs_test[res_nchron$verdict_short == "in-sample violation"], 2),
                  fmt_p(res_nchron$p_value[res_nchron$verdict_short == "in-sample violation"])),
          collapse = "; "),
    paste(nchron_lines, collapse = "; "))
} else {
  sprintf(paste0(
    "**The exclusion restriction is weakly violated on at least one ",
    "outcome but not in-sample-falsified outright.** Per outcome: %s. ",
    "On %s, `n_chron_clean` falls in the `1.5 <= |t| < 2.0` band -- a ",
    "marginal in-sample partial correlation that is below conventional ",
    "rejection thresholds for the exclusion claim, but warrants explicit ",
    "memo-level acknowledgement rather than a clean-pass framing."),
    paste(nchron_lines, collapse = "; "),
    paste(sprintf("`%s`", fail_outcomes_nchron), collapse = ", "))
}

placebo_lines <- vapply(seq_len(nrow(res_placebo)),
                        function(i) fmt_per_outcome(res_placebo[i, ]),
                        character(1))

para2 <- if (all_pass_nchron && !placebo_any_nondef) {
  sprintf(paste0(
    "**Placebo comparison.** The placebo (`maddison_priority`, which is ",
    "in the selection equation as a covariate but is not the headline ",
    "exclusion restriction) is also clean across all three outcomes: %s. ",
    "Both the focal and the placebo pass the in-sample test, which is ",
    "the calibration we would expect: at this sample size (N = %d) the ",
    "complete-case partial-correlation test does not detect spurious ",
    "in-sample violations on either of two non-target variables, so the ",
    "`n_chron_clean` clean-pass result is not an artifact of low test ",
    "power."),
    paste(placebo_lines, collapse = "; "),
    nrow(df_complete))
} else if (all_pass_nchron && placebo_any_nondef) {
  sprintf(paste0(
    "**Placebo comparison.** The placebo (`maddison_priority`, which is ",
    "in the selection equation as a covariate but is not claimed as an ",
    "exclusion restriction) shows non-defensible in-sample partial ",
    "correlation on %s while passing the others (full per-outcome ",
    "placebo: %s). Because `n_chron_clean` itself passes on all three ",
    "outcomes, the placebo's own failure is a calibration point about ",
    "noise at this sample size (N = %d) rather than a finding about the ",
    "exclusion restriction; the headline conclusion is unchanged."),
    paste(sprintf("`%s`", placebo_fails), collapse = ", "),
    paste(placebo_lines, collapse = "; "),
    nrow(df_complete))
} else if (!placebo_any_nondef && !all_pass_nchron) {
  sprintf(paste0(
    "**Placebo comparison.** The placebo (`maddison_priority`) is clean ",
    "across all three outcomes: %s. By contrast, `n_chron_clean` shows ",
    "non-defensible in-sample partial correlation on %s. The placebo / ",
    "focal contrast indicates that the in-sample violation on ",
    "`n_chron_clean` is not a generic noise artifact at this sample size ",
    "(N = %d) but a specific weakness of the chosen exclusion variable; ",
    "this needs memo-level acknowledgment."),
    paste(placebo_lines, collapse = "; "),
    paste(sprintf("`%s`", nchron_fails), collapse = ", "),
    nrow(df_complete))
} else if (length(nchron_unique) == 0L) {
  sprintf(paste0(
    "**Placebo comparison.** The placebo (`maddison_priority`) shows ",
    "in-sample partial correlation of similar or larger magnitude on the ",
    "same outcome(s) that flagged for `n_chron_clean` (shared failures: ",
    "%s). Because the placebo is not claimed as an exclusion restriction, ",
    "this co-failure pattern is a calibration point: at this sample size ",
    "(N = %d) in-sample partial correlations on the chosen covariates ",
    "are noisy enough that the `n_chron_clean` result should be ",
    "interpreted as consistent with noise rather than as a structural ",
    "in-sample failure. Full per-outcome placebo: %s."),
    paste(sprintf("`%s`", shared_fails), collapse = ", "),
    nrow(df_complete),
    paste(placebo_lines, collapse = "; "))
} else {
  shared_text <- if (length(shared_fails) == 0L) {
    paste0("none (the placebo's own failure(s) on ",
           paste(sprintf("`%s`", placebo_fails), collapse = ", "),
           " do not coincide with `n_chron_clean`'s failure(s))")
  } else {
    paste(sprintf("`%s`", shared_fails), collapse = ", ")
  }
  sprintf(paste0(
    "**Placebo comparison.** The placebo (`maddison_priority`) is ",
    "non-defensible on %s; `n_chron_clean` is non-defensible on %s. ",
    "Shared failures between focal and placebo: %s. Where ",
    "`n_chron_clean` flags but the placebo passes (specifically: %s), ",
    "the exclusion variable carries a specific weakness that the placebo ",
    "does not, which is the substantive piece to flag in the memo: this ",
    "is the spec's literal-reading 'placebo clean, focal not clean' case ",
    "on those outcomes, even though the placebo has its own unrelated ",
    "failure on %s. Full per-outcome placebo: %s. Sample size N = %d."),
    paste(sprintf("`%s`", placebo_fails), collapse = ", "),
    paste(sprintf("`%s`", nchron_fails), collapse = ", "),
    shared_text,
    paste(sprintf("`%s`", nchron_unique), collapse = ", "),
    paste(sprintf("`%s`", setdiff(placebo_fails, nchron_fails)), collapse = ", "),
    paste(placebo_lines, collapse = "; "),
    nrow(df_complete))
}

note <- paste0(
  "**Note.** This is a within-sample observable-channel diagnostic only. ",
  "No partial-correlation test of this form can prove that ",
  "`n_chron_clean` is excludable on unobservable channels; the exclusion ",
  "restriction is, in the end, an identifying assumption, not a ",
  "testable hypothesis. The result here strengthens or weakens the ",
  "argument but does not settle it. Triangulation with the Prompt 3d ",
  "balance test (archival vs severity proxies dominating the OLS of ",
  "`n_chron_clean`) and the 3d falsification test (|Δ INCOME| = 0.0029 ",
  "when `n_chron_clean` is dropped from the Heckman selection equation) ",
  "remains the substantive defense; this diagnostic is one further ",
  "consistency check on top of those."
)

lines <- c(lines, "",
  "## Interpretation",
  "",
  para1,
  "",
  para2,
  "",
  note,
  ""
)

out_path <- here::here("output", "tables", "03b_insample_exclusion_test.md")
writeLines(lines, out_path)
cli::cli_alert_success("Wrote {.path output/tables/03b_insample_exclusion_test.md}")
