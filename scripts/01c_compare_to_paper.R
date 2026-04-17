# scripts/01c_compare_to_paper.R
# Build the side-by-side comparison of our replicated coefficients against
# the published Metrick & Schmelzing (2024) Tables 2 and 3, render the
# comparison tables, and emit the full six- and seven-column regression
# tables with iid + cluster-iso3 SEs side-by-side via modelsummary.
#
# Replication threshold (per project memo):
#   replicated == TRUE  iff  |delta_coef| <= 0.02 AND |delta_n| <= 5
#
# Sign-flip rule:
#   For any spec with a non-NA paper_coef, sign(ours_coef) MUST equal
#   sign(paper_coef). A flip stops the script (replication is broken).

source(here::here("R", "utils.R"))
load_libs()

library(fixest)
library(sandwich)
library(lmtest)
library(modelsummary)
library(kableExtra)
library(checkmate)
library(tibble)
library(tidyr)
library(purrr)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Load the previously fit models
# ---------------------------------------------------------------------------

t3_path <- here::here("data", "processed", "table3_ols_fits.rds")
t2_path <- here::here("data", "processed", "table2_logit_fits.rds")
checkmate::assert_file_exists(t3_path)
checkmate::assert_file_exists(t2_path)

t3 <- readRDS(t3_path)        # $models (length 6), $spec_a, $spec_b
t2 <- readRDS(t2_path)        # $models (length 7), $warnings_per_spec, ...

stopifnot(length(t3$models) == 6L, length(t2$models) == 7L)

# ---------------------------------------------------------------------------
# HARD-CODED PAPER VALUES
# ---------------------------------------------------------------------------
# Coefficients are the INCOME (log lagged real GDP per capita) point estimates
# and standard errors as printed in the paper. Values for which the paper
# table cell is not extractable from the public PDF text layer are recorded
# as NA with a `notes` flag so the tibble surfaces a "needs paper value" row
# rather than silently producing a misleading delta.

paper_table3 <- list(
  "Combined, no fiscal" = list(
    # Table 3, col 1: full sample, no fiscal regressors. Paper p. 32.
    coef = -0.2287, se = 0.0370, n = 334L,
    src  = "M&S (2024) Table 3 col 1, paper p. 32"
  ),
  "Combined, + fiscal" = list(
    # Table 3, col 2: full sample including DEBT/EXP. Paper p. 32.
    coef = -0.2672, se = 0.0426, n = 261L,
    src  = "M&S (2024) Table 3 col 2, paper p. 32"
  ),
  "Canonical, no fiscal" = list(
    # Table 3, col 3: canonical (Candidate=0), no fiscal. Paper p. 32.
    # Coefficient from script header; SE/N not in extractable PDF text.
    coef = -0.2306, se = NA_real_, n = NA_integer_,
    src  = "M&S (2024) Table 3 col 3, paper p. 32 (SE/N TODO from print copy)"
  ),
  "Canonical, + fiscal" = list(
    # Table 3, col 4: canonical, + fiscal. Paper p. 32.
    coef = -0.2815, se = NA_real_, n = NA_integer_,
    src  = "M&S (2024) Table 3 col 4, paper p. 32 (SE/N TODO from print copy)"
  ),
  "Candidate, no fiscal" = list(
    # Table 3, col 5: candidate (Candidate=1), no fiscal. Paper p. 32.
    coef = -0.3621, se = NA_real_, n = NA_integer_,
    src  = "M&S (2024) Table 3 col 5, paper p. 32 (SE/N TODO from print copy)"
  ),
  "Candidate, + fiscal" = list(
    # Table 3, col 6: candidate, + fiscal. Paper p. 32.
    coef = -0.2783, se = NA_real_, n = NA_integer_,
    src  = "M&S (2024) Table 3 col 6, paper p. 32 (SE/N TODO from print copy)"
  )
)

# Table 2 paper-value scaffolding. Per paper p. 30, N = 273 across all 7
# columns (complete-case sample). The paper text identifies INCOME as
# positive-and-significant for guarantees (col 1) and lending (col 2);
# the seven INCOME point estimates and SEs are not extractable from the
# PDF text layer, so we record NA + a TODO note. Sign for guarantees and
# lending is hard-coded so the sign-flip guard still bites where the
# paper is explicit.
paper_table2 <- list(
  "guarantees_d"        = list(coef = NA_real_, se = NA_real_, n = 273L,
                               sign = +1L,
                               src = "M&S (2024) Table 2 col 1, paper p. 30 (sign documented; coef TODO)"),
  "lending_d"           = list(coef = NA_real_, se = NA_real_, n = 273L,
                               sign = +1L,
                               src = "M&S (2024) Table 2 col 2, paper p. 30 (sign documented; coef TODO)"),
  "capital_injections_d"= list(coef = NA_real_, se = NA_real_, n = 273L,
                               sign = NA_integer_,
                               src = "M&S (2024) Table 2 col 3, paper p. 30 (coef/SE TODO)"),
  "restructuring_d"     = list(coef = NA_real_, se = NA_real_, n = 273L,
                               sign = NA_integer_,
                               src = "M&S (2024) Table 2 col 4, paper p. 30 (coef/SE TODO)"),
  "asset_management_d"  = list(coef = NA_real_, se = NA_real_, n = 273L,
                               sign = NA_integer_,
                               src = "M&S (2024) Table 2 col 5, paper p. 30 (coef/SE TODO)"),
  "rules_d"             = list(coef = NA_real_, se = NA_real_, n = 273L,
                               sign = NA_integer_,
                               src = "M&S (2024) Table 2 col 6, paper p. 30 (coef/SE TODO)"),
  "other_d"             = list(coef = NA_real_, se = NA_real_, n = 273L,
                               sign = NA_integer_,
                               src = "M&S (2024) Table 2 col 7, paper p. 30 (coef/SE TODO)")
)

# ---------------------------------------------------------------------------
# Helpers to extract our coef / se / n
# ---------------------------------------------------------------------------

extract_one <- function(fit, term = "log_lagged_gdppc") {
  cm <- summary(fit)$coefficients
  if (!term %in% rownames(cm)) {
    return(list(coef = NA_real_, se = NA_real_, n = NA_integer_))
  }
  list(
    coef = unname(cm[term, "Estimate"]),
    se   = unname(cm[term, "Std. Error"]),
    n    = as.integer(stats::nobs(fit))
  )
}

build_row <- function(spec_name, fit, paper_entry, table_name,
                      warns = character(), converged = TRUE) {
  ours <- extract_one(fit)
  pcoef <- paper_entry$coef
  pse   <- paper_entry$se
  pn    <- paper_entry$n

  delta_coef <- if (is.na(pcoef)) NA_real_ else ours$coef - pcoef
  delta_se   <- if (is.na(pse))   NA_real_ else ours$se   - pse
  delta_n    <- if (is.na(pn))    NA_integer_ else as.integer(ours$n - pn)

  replicated <- if (is.na(delta_coef) || is.na(delta_n)) {
    NA
  } else {
    abs(delta_coef) <= 0.02 && abs(delta_n) <= 5L
  }

  notes_v <- character()
  if (length(warns) > 0L)        notes_v <- c(notes_v, "glm-warning")
  if (!isTRUE(converged))        notes_v <- c(notes_v, "non-converged")
  if (is.na(pcoef))              notes_v <- c(notes_v, "paper-coef-NA")
  if (is.na(pn))                 notes_v <- c(notes_v, "paper-n-NA")
  notes <- if (length(notes_v) == 0L) "" else paste(notes_v, collapse = ";")

  tibble::tibble(
    table       = table_name,
    spec        = spec_name,
    target_coef = "log_lagged_gdppc",
    ours_coef   = ours$coef,
    ours_se     = ours$se,
    ours_n      = ours$n,
    paper_coef  = pcoef,
    paper_se    = pse,
    paper_n     = pn,
    delta_coef  = delta_coef,
    delta_se    = delta_se,
    delta_n     = delta_n,
    replicated  = replicated,
    notes       = notes
  )
}

# ---------------------------------------------------------------------------
# Comparison tibble: Table 3
# ---------------------------------------------------------------------------

cmp_t3 <- purrr::map2_dfr(
  names(t3$models),
  t3$models,
  function(nm, fit) build_row(nm, fit, paper_table3[[nm]], "Table 3")
)

# Comparison tibble: Table 2
cmp_t2 <- purrr::map2_dfr(
  names(t2$models),
  t2$models,
  function(nm, fit) {
    build_row(nm, fit, paper_table2[[nm]], "Table 2",
              warns = t2$warnings_per_spec[[nm]] %||% character(),
              converged = t2$converged_per_spec[[nm]])
  }
)

# ---------------------------------------------------------------------------
# SIGN-FLIP GUARD
# ---------------------------------------------------------------------------
# For Table 3 we have explicit paper signs (negative across all 6 specs).
# For Table 2 we encode known signs in `paper_table2[[*]]$sign` (positive
# for guarantees and lending; NA elsewhere). Compare only where known.

flips <- character()

for (i in seq_len(nrow(cmp_t3))) {
  pcoef <- cmp_t3$paper_coef[i]
  if (!is.na(pcoef) && sign(cmp_t3$ours_coef[i]) != sign(pcoef)) {
    flips <- c(flips, sprintf("[%s / %s] ours=%+.4f vs paper=%+.4f",
                              cmp_t3$table[i], cmp_t3$spec[i],
                              cmp_t3$ours_coef[i], pcoef))
  }
}

for (i in seq_len(nrow(cmp_t2))) {
  spec <- cmp_t2$spec[i]
  paper_sign <- paper_table2[[spec]]$sign
  if (!is.na(paper_sign) && sign(cmp_t2$ours_coef[i]) != sign(paper_sign)) {
    flips <- c(flips, sprintf("[%s / %s] ours=%+.4f vs paper sign=%+d",
                              cmp_t2$table[i], spec,
                              cmp_t2$ours_coef[i], paper_sign))
  }
}

if (length(flips) > 0L) {
  msg <- paste(c("INCOME sign flipped on:", paste0("  ", flips)),
               collapse = "\n")
  stop(msg, call. = FALSE)
} else {
  cat("INCOME sign-flip guard: PASSED (",
      sum(!is.na(cmp_t3$paper_coef)) +
        sum(!vapply(paper_table2, function(x) is.na(x$sign), logical(1))),
      " specs checked).\n", sep = "")
}

# ---------------------------------------------------------------------------
# Print to console
# ---------------------------------------------------------------------------

cat("\n=== Table 3 vs paper ===\n"); print(cmp_t3, n = nrow(cmp_t3))
cat("\n=== Table 2 vs paper ===\n"); print(cmp_t2, n = nrow(cmp_t2))

# ---------------------------------------------------------------------------
# Render comparison tibbles to .md and .html
# ---------------------------------------------------------------------------

out_dir <- here::here("output", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

format_num <- function(x, digits = 4) {
  ifelse(is.na(x), "—", formatC(x, format = "f", digits = digits))
}
format_int <- function(x) {
  ifelse(is.na(x), "—", formatC(x, format = "d"))
}

format_replicated <- function(x) {
  ifelse(is.na(x), "needs review",
         ifelse(x, "replicated", "drift"))
}

format_for_render <- function(df) {
  df |>
    dplyr::mutate(
      ours_coef  = format_num(ours_coef),
      ours_se    = format_num(ours_se),
      ours_n     = format_int(ours_n),
      paper_coef = format_num(paper_coef),
      paper_se   = format_num(paper_se),
      paper_n    = format_int(paper_n),
      delta_coef = format_num(delta_coef, digits = 4),
      delta_se   = format_num(delta_se,   digits = 4),
      delta_n    = format_int(delta_n),
      status     = format_replicated(replicated)
    ) |>
    dplyr::select(spec, target_coef, ours_coef, ours_se, ours_n,
                  paper_coef, paper_se, paper_n,
                  delta_coef, delta_se, delta_n, status, notes)
}

write_md_table <- function(df, path, title) {
  rendered <- format_for_render(df)
  md <- knitr::kable(rendered, format = "pipe",
                     caption = title, align = "l")
  writeLines(as.character(md), path)
}

write_html_table <- function(df, path, title) {
  rendered <- format_for_render(df)
  kt <- kableExtra::kbl(rendered, format = "html", caption = title) |>
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))
  for (i in seq_len(nrow(df))) {
    rep_i <- df$replicated[i]
    if (isFALSE(rep_i)) {
      kt <- kt |> kableExtra::row_spec(i, background = "#fde2e2")  # drift
    } else if (is.na(rep_i)) {
      kt <- kt |> kableExtra::row_spec(i, background = "#fff3cd")  # needs review
    }
  }
  writeLines(as.character(kt), path)
}

write_md_table(cmp_t3, file.path(out_dir, "01_table3_vs_paper.md"),
               "Table 3 vs paper (INCOME coefficient)")
write_html_table(cmp_t3, file.path(out_dir, "01_table3_vs_paper.html"),
                 "Table 3 vs paper (INCOME coefficient)")
write_md_table(cmp_t2, file.path(out_dir, "01_table2_vs_paper.md"),
               "Table 2 vs paper (INCOME coefficient)")
write_html_table(cmp_t2, file.path(out_dir, "01_table2_vs_paper.html"),
                 "Table 2 vs paper (INCOME coefficient)")

# ---------------------------------------------------------------------------
# Side-by-side regression tables: iid + cluster-iso3 SEs
# ---------------------------------------------------------------------------
#
# We must reload the analysis CSV here because the saved lm/glm fits do NOT
# carry iso3 in their stored model frame (iso3 is not in the regression
# formula), and `insight::get_data()` cannot reconstruct it from the call
# alone after fits are deserialized in a fresh session. We refit Table 3
# with fixest::feols on the same formulas (per project memo: use feols for
# clustered SEs on continuous outcomes; do not hand-code), and pre-compute
# clustered vcov matrices for the Table 2 logits via sandwich::vcovCL.

df_t3 <- read.csv(here::here("data", "processed", "analysis_data.csv")) |>
  dplyr::filter(year >= 1800)
df_t2 <- read.csv(here::here("data", "processed", "analysis_data.csv")) |>
  dplyr::filter(r == 1)

# ---- Table 3: refit with feols on each subsample so cluster ~ iso3 works
feols_subsamples <- list(
  "Combined, no fiscal"  = df_t3,
  "Combined, + fiscal"   = df_t3,
  "Canonical, no fiscal" = df_t3 |> dplyr::filter(candidate == 0),
  "Canonical, + fiscal"  = df_t3 |> dplyr::filter(candidate == 0),
  "Candidate, no fiscal" = df_t3 |> dplyr::filter(candidate == 1),
  "Candidate, + fiscal"  = df_t3 |> dplyr::filter(candidate == 1)
)
feols_specs <- list(
  "Combined, no fiscal"  = t3$spec_a,
  "Combined, + fiscal"   = t3$spec_b,
  "Canonical, no fiscal" = t3$spec_a,
  "Canonical, + fiscal"  = t3$spec_b,
  "Candidate, no fiscal" = t3$spec_a,
  "Candidate, + fiscal"  = t3$spec_b
)

t3_feols <- mapply(
  function(spec, data) fixest::feols(spec, data = data),
  feols_specs, feols_subsamples,
  SIMPLIFY = FALSE
)

# Spot-check that feols and lm agree on point estimates
for (nm in names(t3_feols)) {
  b_lm <- coef(t3$models[[nm]])["log_lagged_gdppc"]
  b_fe <- coef(t3_feols[[nm]])["log_lagged_gdppc"]
  if (abs(b_lm - b_fe) > 1e-8) {
    stop(sprintf("feols/lm disagree on %s INCOME: %.6f vs %.6f",
                 nm, b_fe, b_lm), call. = FALSE)
  }
}

build_dual_t3 <- function(fits) {
  doubled <- list(); vcovs <- list()
  for (nm in names(fits)) {
    doubled[[paste0(nm, " (IID)")]]           <- fits[[nm]]
    doubled[[paste0(nm, " (cluster:iso3)")]]  <- fits[[nm]]
    vcovs <- c(vcovs, list("iid"), list(~iso3))
  }
  list(models = doubled, vcov = vcovs)
}

coef_rename <- c(
  "log_lagged_gdppc" = "INCOME (log)",
  "currency_regime" = "FX regime",
  "polity"          = "Polity",
  "candidate"       = "Candidate",
  "debt"            = "DEBT/GDP",
  "exp"             = "EXP/GDP",
  "year"            = "Year",
  "(Intercept)"     = "(Intercept)"
)

t3_args <- build_dual_t3(t3_feols)

ms_t3 <- modelsummary::modelsummary(
  t3_args$models,
  vcov = t3_args$vcov,
  stars = TRUE,
  gof_omit = "AIC|BIC|Log.Lik|F|RMSE|Std.Errors",
  coef_rename = coef_rename,
  output = "html"
)
writeLines(as.character(ms_t3),
           file.path(out_dir, "01_table3_six_models.html"))

# ---- Table 2: precompute clustered vcov matrices for each logit
# sandwich::vcovCL aligns the cluster vector with the fit's na.action.
t2_clustered_vcov <- lapply(t2$models, function(fit) {
  sandwich::vcovCL(fit, cluster = df_t2$iso3, type = "HC0")
})

t2_doubled <- list()
t2_vcovs   <- list()
for (nm in names(t2$models)) {
  t2_doubled[[paste0(nm, " (IID)")]]           <- t2$models[[nm]]
  t2_doubled[[paste0(nm, " (cluster:iso3)")]]  <- t2$models[[nm]]
  t2_vcovs <- c(t2_vcovs, list("iid"), list(t2_clustered_vcov[[nm]]))
}

ms_t2 <- modelsummary::modelsummary(
  t2_doubled,
  vcov = t2_vcovs,
  stars = TRUE,
  gof_omit = "AIC|BIC|Log.Lik|F|RMSE|Deviance|Std.Errors",
  coef_rename = coef_rename,
  output = "html"
)
writeLines(as.character(ms_t2),
           file.path(out_dir, "01_table2_seven_models.html"))

# ---------------------------------------------------------------------------
# Summary line
# ---------------------------------------------------------------------------

all_repl <- c(cmp_t3$replicated, cmp_t2$replicated)
n_repl   <- sum(all_repl == TRUE,  na.rm = TRUE)
n_drift  <- sum(all_repl == FALSE, na.rm = TRUE)
n_review <- sum(is.na(all_repl))

cat(sprintf(
  "\nReplication summary: replicated=%d  drift=%d  needs-review=%d  total=%d\n",
  n_repl, n_drift, n_review, length(c(cmp_t3$replicated, cmp_t2$replicated))
))

cat("\nWritten:\n")
for (f in c(
  "01_table3_vs_paper.md", "01_table3_vs_paper.html",
  "01_table2_vs_paper.md", "01_table2_vs_paper.html",
  "01_table3_six_models.html", "01_table2_seven_models.html"
)) {
  cat("  ", file.path(out_dir, f), "\n")
}
