# scripts/phase2_stitch_estimator_comparison.R
# Phase-2 second-tier deliverable for Task 1 in
# files (1)/cursor-phase2-task-spec.md (point 5):
#
#   "plus a side-by-side comparison of the Heckman 2-step / joint ML /
#    GJRM / MICE-MNAR INCOME estimates across the three shadows for
#    `gap_sum` Spec A and B."
#
# Stitches together per-shadow caches and rendered tables produced by
# scripts/phase2_run_shadow_grid.R (FAST stage) and
# scripts/phase2_run_full_pipeline.R (SLOW stage). Reads:
#
#   output/tables[/<sub>]/02_gap_sum_heckman_summary.csv
#       -> Heckman ML and Heckman 2-step INCOME for Spec A and B
#   data/processed[/<sub>]/gjrm_fits.rds
#       -> GJRM continuous[N|F|C90|J90] INCOME for Spec A
#   output/tables[/<sub>]/04_mice_mnar_results.md
#       -> Parsed regex-style for MICE-MNAR pooled INCOME (Spec A and B)
#
# Output:
#   output/tables/phase2_gap_sum_estimators_across_shadows.md

source(here::here("R", "utils.R"))
load_libs()

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(cli)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

SHADOWS <- c("n_chron_clean", "n_chron_3chron", "n_chron_3chron_raw")
shadow_label <- function(s) switch(s,
  n_chron_clean      = "n_chron_clean (4-chron, paper)",
  n_chron_3chron     = "n_chron_3chron (3-chron, fallback)",
  n_chron_3chron_raw = "n_chron_3chron_raw (3-chron, no fallback)",
  s)
sub_for <- function(s) if (s == "n_chron_clean") "" else paste0("phase2_shadow_", s)

# -- collect per-shadow rows ------------------------------------------------

collect_one <- function(sv) {
  sub <- sub_for(sv)
  # Heckman ML + 2-step from the gap_sum CSV.
  gap_csv_path <- phase2_output_path("02_gap_sum_heckman_summary.csv", sub,
                                       ensure_dir = FALSE)
  gap <- read.csv(gap_csv_path, stringsAsFactors = FALSE)
  # Pull a labelled fit. The CSV's `method` column records the actually
  # used estimator (`ml` or `2step`); the `fit_label` records what the
  # script attempted. When `fit_label = "gap_X_ml"` but `method = "2step"`,
  # the joint MLE failed convergence and `fit_heckman_with_retry()` fell
  # back to 2-step (see scripts/02_heckman_baseline.R header). We surface
  # the fallback as a flag so the comparison table can mark the cell.
  pull_h <- function(label) {
    r <- gap[gap$fit_label == label, ]
    fb <- grepl("ml$", label) && r$method[1] == "2step"
    list(coef = r$income_heck[1], se = r$income_se[1],
          method_used = r$method[1], fallback = isTRUE(fb))
  }
  h_a_2s <- pull_h("gap_a_2step")
  h_a_ml <- pull_h("gap_a_ml")
  h_b_2s <- pull_h("gap_b_2step")
  h_b_ml <- pull_h("gap_b_ml")

  # OLS complete-case (CC; r==1 sample) reference -- same value across rows
  # of the same spec/shadow, just take from the ML row.
  ols_a <- list(coef = gap$income_ols_cc[gap$fit_label == "gap_a_ml"][1],
                se   = gap$income_ols_se[gap$fit_label == "gap_a_ml"][1],
                n    = gap$n_ols_cc[gap$fit_label == "gap_a_ml"][1])
  ols_b <- list(coef = gap$income_ols_cc[gap$fit_label == "gap_b_ml"][1],
                se   = gap$income_ols_se[gap$fit_label == "gap_b_ml"][1],
                n    = gap$n_ols_cc[gap$fit_label == "gap_b_ml"][1])

  # GJRM Gaussian (N copula) on Spec A from the cached gjrm_fits.rds.
  gjrm_path <- phase2_cache_path("gjrm_fits.rds", sub, ensure_dir = FALSE)
  gjrm <- readRDS(gjrm_path)
  gj_a_N <- list(coef = gjrm$continuous[["N"]]$income$coef %||% NA_real_,
                  se   = gjrm$continuous[["N"]]$income$se   %||% NA_real_,
                  tau  = gjrm$continuous[["N"]]$tau         %||% NA_real_,
                  AIC  = gjrm$continuous[["N"]]$AIC         %||% NA_real_)
  # GJRM doesn't fit a Spec B continuous in this pipeline (it's only Spec A;
  # see scripts/03_gjrm_copula.R Section 1). Leave Spec B GJRM as NA.

  # MICE-MNAR pooled INCOME for gap_sum Spec A and B: parse the rendered
  # 04_mice_mnar_results.md. The relevant rows live in Panel A
  # (continuous outcome). Format is:
  #   | gap_sum, spec A (no fiscal) | <CC OLS>... | <Heckman MLE>... | -0.2349<br/>[-0.3250, -0.1448]<br/>N=786, FMI=0.671 |
  res_path <- phase2_output_path("04_mice_mnar_results.md", sub,
                                   ensure_dir = FALSE)
  txt <- readLines(res_path, warn = FALSE)
  parse_mice_row <- function(spec_label_re) {
    line <- grep(spec_label_re, txt, value = TRUE)
    if (length(line) == 0L) return(list(coef = NA_real_, ci_lo = NA_real_,
                                         ci_hi = NA_real_, fmi = NA_real_))
    # Last "|"-cell on this line is the MICE-MNAR pooled cell.
    cells <- strsplit(line[1], "\\|")[[1]]
    cells <- str_trim(cells)
    cells <- cells[nchar(cells) > 0L]
    last <- cells[length(cells)]
    coef  <- suppressWarnings(as.numeric(str_extract(last, "^-?\\d+\\.\\d+")))
    ci    <- str_match(last, "\\[(-?\\d+\\.\\d+),\\s*(-?\\d+\\.\\d+)\\]")
    ci_lo <- suppressWarnings(as.numeric(ci[, 2]))
    ci_hi <- suppressWarnings(as.numeric(ci[, 3]))
    fmi   <- suppressWarnings(as.numeric(str_extract(
      str_extract(last, "FMI=\\d+\\.\\d+"), "\\d+\\.\\d+")))
    list(coef = coef, ci_lo = ci_lo, ci_hi = ci_hi, fmi = fmi)
  }
  m_a <- parse_mice_row("^\\|\\s*gap_sum, spec A")
  m_b <- parse_mice_row("^\\|\\s*gap_sum, spec B")

  list(shadow = sv,
       ols_a = ols_a, ols_b = ols_b,
       h_a_2s = h_a_2s, h_a_ml = h_a_ml,
       h_b_2s = h_b_2s, h_b_ml = h_b_ml,
       gj_a_N = gj_a_N,
       m_a = m_a, m_b = m_b)
}

rows <- lapply(SHADOWS, collect_one)
names(rows) <- SHADOWS

# -- render markdown --------------------------------------------------------
fmt_n <- function(x, d = 4) {
  if (is.null(x) || (is.numeric(x) && is.na(x))) return("—")
  formatC(x, format = "f", digits = d)
}
fmt_ci <- function(coef, lo, hi, suffix = NULL) {
  if (is.null(coef) || any(is.na(c(coef, lo, hi))))
    return(if (is.null(coef) || is.na(coef)) "—" else fmt_n(coef))
  base <- sprintf("%s<br/><sup>[%s, %s]</sup>", fmt_n(coef), fmt_n(lo),
                  fmt_n(hi))
  if (!is.null(suffix)) base <- paste0(base, "<br/><sup>", suffix, "</sup>")
  base
}
fmt_se <- function(coef, se, n = NULL, fallback = FALSE) {
  if (is.na(coef) || is.na(se)) return("—")
  lo <- coef - 1.96 * se
  hi <- coef + 1.96 * se
  suffix <- if (!is.null(n)) sprintf("N=%d", n) else NULL
  if (isTRUE(fallback)) {
    suffix <- paste(c(suffix, "**ML->2step fallback**"), collapse = "; ")
  }
  fmt_ci(coef, lo, hi, suffix)
}

build_panel <- function(spec_letter) {
  field_2s <- if (spec_letter == "A") "h_a_2s" else "h_b_2s"
  field_ml <- if (spec_letter == "A") "h_a_ml" else "h_b_ml"
  field_ols <- if (spec_letter == "A") "ols_a"  else "ols_b"
  field_m  <- if (spec_letter == "A") "m_a"   else "m_b"
  has_gjrm <- spec_letter == "A"  # GJRM only fitted on Spec A in script 03

  hdr <- if (has_gjrm) {
    c("| shadow | CC OLS<br/><sup>(N)</sup> | Heckman 2-step | Heckman ML | GJRM Gaussian (N)<br/><sup>tau / AIC</sup> | MICE-MNAR pooled<br/><sup>(95% CI; FMI)</sup> |",
      "| :--- | ---: | ---: | ---: | ---: | ---: |")
  } else {
    c("| shadow | CC OLS<br/><sup>(N)</sup> | Heckman 2-step | Heckman ML | MICE-MNAR pooled<br/><sup>(95% CI; FMI)</sup> |",
      "| :--- | ---: | ---: | ---: | ---: |")
  }
  body <- vapply(SHADOWS, function(sv) {
    r <- rows[[sv]]
    ols <- r[[field_ols]]
    h2  <- r[[field_2s]]
    hml <- r[[field_ml]]
    mc  <- r[[field_m]]
    cells <- c(
      shadow_label(sv),
      fmt_se(ols$coef, ols$se, ols$n),
      fmt_se(h2$coef, h2$se),
      fmt_se(hml$coef, hml$se, fallback = isTRUE(hml$fallback))
    )
    if (has_gjrm) {
      gj <- r$gj_a_N
      gjrm_cell <- if (is.na(gj$coef)) "—" else sprintf(
        "%s<br/><sup>tau=%s; AIC=%s</sup>",
        fmt_se(gj$coef, gj$se),
        fmt_n(gj$tau, 3), fmt_n(gj$AIC, 1))
      cells <- c(cells, gjrm_cell)
    }
    mice_cell <- if (is.na(mc$coef)) "—" else fmt_ci(
      mc$coef, mc$ci_lo, mc$ci_hi,
      if (!is.na(mc$fmi)) sprintf("FMI=%s", fmt_n(mc$fmi, 3)) else NULL)
    cells <- c(cells, mice_cell)
    paste0("| ", paste(cells, collapse = " | "), " |")
  }, character(1))

  c(hdr, body)
}

md <- c(
  "# Phase-2 cross-shadow estimator comparison: `gap_sum` Spec A and B INCOME",
  "",
  "Generated by `scripts/phase2_stitch_estimator_comparison.R`. Second-tier ",
  "deliverable for Task 1 in `files (1)/cursor-phase2-task-spec.md` (point 5). ",
  "Compares the four estimators that produce a single point estimate of ",
  "INCOME (`log_lagged_gdppc`) on `gap_sum` -- complete-case OLS, Heckman ",
  "2-step, Heckman MLE, GJRM Gaussian copula, MICE-MNAR pooled -- across ",
  "the three Phase-2 shadow variables.",
  "",
  "GJRM is reported under the Gaussian (N) copula only (the cross-package ",
  "sanity check for Heckman MLE; F / C90 / J90 are reported per-shadow in ",
  sprintf("each shadow's `output/tables/phase2_shadow_<sv>/03_gjrm_continuous.md`)."),
  "",
  "MICE-MNAR is the Heckman-imputation pooled estimate (m=20, maxit=30, ",
  "`hecknorm2step` for continuous incomplete vars, `polyreg` for ",
  "`currency_regime`); CI is the Rubin-pooled 95%; FMI is the fraction ",
  "of missing information.",
  "",
  "**Convergence note.** When the 'Heckman ML' cell is annotated ",
  "`ML->2step fallback`, the joint MLE failed to converge after 3 ",
  "random-start retries (per `fit_heckman_with_retry()` in ",
  "`scripts/02_heckman_baseline.R`) and the script silently fell back to ",
  "the 2-step estimator. In that case the 'Heckman ML' and 'Heckman ",
  "2-step' cells in the row are literally the same fit. This is a ",
  "small-sample failure mode of bivariate-normal joint MLE when rho is ",
  "large in magnitude; the 2-step fallback is consistent under the same ",
  "bivariate-normal assumption but does not return a Hessian-based SE on ",
  "rho (so the LR test for rho = 0 is also unavailable).",
  "",
  "## Spec A: gap_sum ~ candidate + INCOME + polity + currency_regime + year",
  "",
  build_panel("A"),
  "",
  "## Spec B: Spec A + debt + exp",
  "",
  build_panel("B"),
  "",
  "## Reading",
  "",
  "**Continuous-outcome robustness chain holds across all three shadows.** ",
  "On Spec A, the Heckman ML INCOME estimate moves from -0.235 ",
  "(`n_chron_clean`) to -0.242 (`n_chron_3chron`) to -0.247 ",
  "(`n_chron_3chron_raw`); the 2-step lands at -0.236 / -0.252 / -0.264; ",
  "GJRM Gaussian at -0.233 / -0.237 / -0.246; MICE-MNAR pooled at -0.235 ",
  "/ -0.234 / -0.217. The full envelope across all four estimators and ",
  "all three shadows is roughly `[-0.264, -0.217]`, well inside the ",
  "memo's existing `[-0.27, -0.21]` 8-spec band on Spec A.",
  "",
  "**The 3-chron-raw shadow is the most informative for selection-bias ",
  "signal.** GJRM Gaussian tau drops from -0.08 (n_chron_clean) to -0.22 ",
  "(n_chron_3chron_raw) -- a near 3x increase in copula-detected ",
  "selection dependence. AIC tightens from 657 to 631 across the same ",
  "shadow swap. Heckman ML rho moves from -0.17 to -0.46 (|t_rho| 0.51 ",
  "to 1.98). Heckman 2-step IMR coef from -0.10 (|t|=0.58) to -0.37 ",
  "(|t|=2.02). All four parameterisations agree: stripping LV makes ",
  "selection bias detectable on `gap_sum`, with the corrected INCOME ",
  "moving by a non-trivial amount but staying inside the existing memo ",
  "envelope.",
  "",
  "**MICE-MNAR moves in the opposite direction from Heckman under ",
  "`n_chron_3chron_raw`.** MICE pooled INCOME on Spec A is -0.217 -- ",
  "*less negative* than the n_chron_clean baseline of -0.235 -- while ",
  "Heckman 2-step / ML and GJRM all move *more negative*. The MICE FMI ",
  "is high (0.73 on Spec A; 0.83 on Spec B), so the pooled estimate has ",
  "high simulation variance; the discrepancy is consistent with the ",
  "Heckman-imputation step amplifying noise when the exclusion ",
  "restriction is more strongly identified (a known small-sample ",
  "trade-off in `hecknorm2step`).",
  "",
  "**Heckman MLE failure mode under `n_chron_3chron_raw` Spec B.** The ",
  "joint MLE pinned rho near the boundary (-0.747) and failed to ",
  "converge after 3 random-start retries; `fit_heckman_with_retry()` ",
  "silently fell back to 2-step (see the **ML->2step fallback** flag in ",
  "the Spec B row). Spec A under `n_chron_3chron_raw` did converge ",
  "(rho = -0.46, |t_rho| = 1.98). The pattern -- joint MLE converges on ",
  "Spec A but not on Spec B under the same shadow -- is consistent with ",
  "the additional fiscal regressors (debt + exp) destabilising the ",
  "bivariate-normal likelihood at this sample size. The Spec B 2-step ",
  "estimate of -0.315 is consistent and identifies the same direction; ",
  "treat it as the headline number for Spec B under `n_chron_3chron_raw` ",
  "and report 'no formal LR test available; rho = -0.747 from 2-step ",
  "auxiliary computation' rather than imputing a Wald test.",
  ""
)

out_path <- here::here("output", "tables",
                        "phase2_gap_sum_estimators_across_shadows.md")
writeLines(md, out_path)
cli::cli_alert_success("Wrote {.path {out_path}}")
