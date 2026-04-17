# scripts/01c_compare_to_paper.R
# Build the side-by-side comparison of our replicated coefficients against
# the published Metrick & Schmelzing (2024) Tables 2 and 3, render the
# comparison tables, and emit the full six- and seven-column regression
# tables with iid + cluster-iso3 SEs side-by-side via modelsummary.
#
# Replication threshold (per project memo):
#   replicated == TRUE  iff  |delta_coef| <= 0.02 AND |delta_n| <= 5
#
# Status taxonomy (cleaner than a single TRUE/FALSE):
#   replicated   - both coef and N within tolerance       (green)
#   partial      - one of coef/N within, the other not    (amber)
#   drift        - both outside tolerance                  (red)
#   needs review - paper value missing for either coef/N   (amber)
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
# All values keyed by hand from the print copy of Metrick & Schmelzing (2024)
# NBER w29281 (revised Feb 2024). Values were verified internally by the
# partition-sum check (canonical_N + candidate_N == combined_N) for both
# Table 3 specifications. INCOME = log lagged real GDP per capita.

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
    coef = -0.2306, se = 0.0416, n = 203L,
    src  = "M&S (2024) Table 3 col 3, paper p. 32"
  ),
  "Canonical, + fiscal" = list(
    # Table 3, col 4: canonical, + fiscal. Paper p. 32.
    coef = -0.2815, se = 0.0585, n = 146L,
    src  = "M&S (2024) Table 3 col 4, paper p. 32"
  ),
  "Candidate, no fiscal" = list(
    # Table 3, col 5: candidate (Candidate=1), no fiscal. Paper p. 32.
    coef = -0.3621, se = 0.0791, n = 131L,
    src  = "M&S (2024) Table 3 col 5, paper p. 32"
  ),
  "Candidate, + fiscal" = list(
    # Table 3, col 6: candidate, + fiscal. Paper p. 32.
    coef = -0.2783, se = 0.0680, n = 115L,
    src  = "M&S (2024) Table 3 col 6, paper p. 32"
  )
)

# Table 2 (paper p. 29): seven logits. N = 273 across all columns
# (complete-case sample on Table 3 regressors). The `sign` field is now
# redundant with sign(coef) but is retained as a belt-and-braces audit
# trail of which signs the paper text discusses explicitly (cols 1-2,
# guarantees and lending, p. 30 narrative).
paper_table2 <- list(
  "guarantees_d"        = list(
    # Table 2, col 1: guarantees indicator. Paper p. 29.
    coef = 0.609,  se = 0.226, n = 273L, sign = +1L,
    src  = "M&S (2024) Table 2 col 1, paper p. 29"
  ),
  "lending_d"           = list(
    # Table 2, col 2: lending indicator. Paper p. 29.
    coef = 0.642,  se = 0.201, n = 273L, sign = +1L,
    src  = "M&S (2024) Table 2 col 2, paper p. 29"
  ),
  "capital_injections_d"= list(
    # Table 2, col 3: capital injections indicator. Paper p. 29.
    coef = 0.255,  se = 0.196, n = 273L, sign = +1L,
    src  = "M&S (2024) Table 2 col 3, paper p. 29"
  ),
  "restructuring_d"     = list(
    # Table 2, col 4: restructuring indicator. Paper p. 29.
    coef = -0.454, se = 0.200, n = 273L, sign = -1L,
    src  = "M&S (2024) Table 2 col 4, paper p. 29"
  ),
  "asset_management_d"  = list(
    # Table 2, col 5: asset management indicator. Paper p. 29.
    coef = 0.408,  se = 0.247, n = 273L, sign = +1L,
    src  = "M&S (2024) Table 2 col 5, paper p. 29"
  ),
  "rules_d"             = list(
    # Table 2, col 6: rules indicator. Paper p. 29.
    coef = -0.014, se = 0.261, n = 273L, sign = -1L,
    src  = "M&S (2024) Table 2 col 6, paper p. 29"
  ),
  "other_d"             = list(
    # Table 2, col 7: other indicator. Paper p. 29.
    coef = -0.148, se = 0.382, n = 273L, sign = -1L,
    src  = "M&S (2024) Table 2 col 7, paper p. 29"
  )
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

classify_status <- function(delta_coef, delta_n,
                            tol_coef = 0.02, tol_n = 5L) {
  if (is.na(delta_coef) || is.na(delta_n)) {
    return(list(replicated = NA, status = "needs review"))
  }
  coef_ok <- abs(delta_coef) <= tol_coef
  n_ok    <- abs(delta_n)    <= tol_n
  if (coef_ok && n_ok)        list(replicated = TRUE,  status = "replicated")
  else if (coef_ok || n_ok)   list(replicated = FALSE, status = "partial")
  else                        list(replicated = FALSE, status = "drift")
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

  cls <- classify_status(delta_coef, delta_n)

  notes_v <- character()
  if (length(warns) > 0L)        notes_v <- c(notes_v, "glm-warning")
  if (!isTRUE(converged))        notes_v <- c(notes_v, "non-converged")
  if (is.na(pcoef))              notes_v <- c(notes_v, "paper-coef-NA")
  if (is.na(pn))                 notes_v <- c(notes_v, "paper-n-NA")
  if (identical(cls$status, "partial")) {
    coef_ok <- abs(delta_coef) <= 0.02
    notes_v <- c(notes_v, if (coef_ok) "N-only-drift" else "coef-only-drift")
  }
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
    replicated  = cls$replicated,
    status      = cls$status,
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
# Now that paper coefficients are populated for all 13 specs, we sign-check
# every spec by comparing sign(ours_coef) to sign(paper_coef). The previous
# `paper_table2[[*]]$sign` field is retained as a redundant audit trail
# (paper-text-documented signs vs hand-keyed table values) but is no longer
# the only source for the check.

cmp_all <- dplyr::bind_rows(cmp_t3, cmp_t2)

flips <- character()
for (i in seq_len(nrow(cmp_all))) {
  pcoef <- cmp_all$paper_coef[i]
  if (!is.na(pcoef) && sign(cmp_all$ours_coef[i]) != sign(pcoef)) {
    flips <- c(flips, sprintf("[%s / %s] ours=%+.4f vs paper=%+.4f",
                              cmp_all$table[i], cmp_all$spec[i],
                              cmp_all$ours_coef[i], pcoef))
  }
}

# Cross-check: documented paper-text signs must agree with hand-keyed
# Table 2 coefficient signs. A mismatch here means either the keyed
# coefficient is wrong or the sign annotation is wrong.
text_sign_mismatches <- character()
for (spec in names(paper_table2)) {
  entry <- paper_table2[[spec]]
  if (!is.na(entry$sign) && !is.na(entry$coef) &&
      sign(entry$coef) != sign(entry$sign)) {
    text_sign_mismatches <- c(text_sign_mismatches,
      sprintf("paper_table2[['%s']]: coef=%+.4f vs declared sign=%+d",
              spec, entry$coef, entry$sign))
  }
}
if (length(text_sign_mismatches) > 0L) {
  stop(paste(c("paper_table2 internal sign inconsistency:",
               paste0("  ", text_sign_mismatches)), collapse = "\n"),
       call. = FALSE)
}

n_sign_checked <- sum(!is.na(cmp_all$paper_coef))
if (length(flips) > 0L) {
  msg <- paste(c("INCOME sign flipped on:", paste0("  ", flips)),
               collapse = "\n")
  stop(msg, call. = FALSE)
} else {
  cat(sprintf("INCOME sign-flip guard: PASSED (%d specs checked).\n",
              n_sign_checked))
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
      delta_n    = format_int(delta_n)
    ) |>
    dplyr::select(spec, target_coef, ours_coef, ours_se, ours_n,
                  paper_coef, paper_se, paper_n,
                  delta_coef, delta_se, delta_n, status, notes)
}

# Color scheme for the rendered HTML tables. The status taxonomy is:
#   replicated   - both |Δβ| ≤ 0.02 AND |ΔN| ≤ 5  (green)
#   partial      - one criterion within tolerance, the other not  (amber)
#   drift        - both criteria outside tolerance  (red)
#   needs review - paper coef or N missing  (amber)
status_color <- function(status) {
  switch(status,
         "replicated"   = "#d4edda",  # soft green
         "partial"      = "#fff3cd",  # amber
         "needs review" = "#fff3cd",  # amber
         "drift"        = "#fde2e2",  # soft red
         NA)
}

caption_with_legend <- function(title) {
  paste0(
    title,
    " — green = replicated (|Δβ| ≤ 0.02 AND |ΔN| ≤ 5); ",
    "amber = partial (one of coef/N within tolerance, the other not, ",
    "or paper value missing); red = drift (both outside tolerance)."
  )
}

write_md_table <- function(df, path, title) {
  rendered <- format_for_render(df)
  md <- knitr::kable(rendered, format = "pipe",
                     caption = caption_with_legend(title), align = "l")
  writeLines(as.character(md), path)
}

write_html_table <- function(df, path, title) {
  rendered <- format_for_render(df)
  kt <- kableExtra::kbl(rendered, format = "html",
                        caption = caption_with_legend(title)) |>
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))
  for (i in seq_len(nrow(df))) {
    bg <- status_color(df$status[i])
    if (!is.na(bg)) kt <- kt |> kableExtra::row_spec(i, background = bg)
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

summarize_status <- function(cmp) {
  c(
    replicated = sum(cmp$status == "replicated"),
    partial    = sum(cmp$status == "partial"),
    drift      = sum(cmp$status == "drift"),
    review     = sum(cmp$status == "needs review"),
    total      = nrow(cmp)
  )
}
sum_t3  <- summarize_status(cmp_t3)
sum_t2  <- summarize_status(cmp_t2)
sum_all <- summarize_status(cmp_all)

cat("\n=== Replication summary ===\n")
cat(sprintf("  Table 3: replicated=%d  partial=%d  drift=%d  review=%d  (total %d)\n",
            sum_t3["replicated"], sum_t3["partial"], sum_t3["drift"],
            sum_t3["review"],  sum_t3["total"]))
cat(sprintf("  Table 2: replicated=%d  partial=%d  drift=%d  review=%d  (total %d)\n",
            sum_t2["replicated"], sum_t2["partial"], sum_t2["drift"],
            sum_t2["review"],  sum_t2["total"]))
cat(sprintf("  Overall: replicated=%d  partial=%d  drift=%d  review=%d  (total %d)\n",
            sum_all["replicated"], sum_all["partial"], sum_all["drift"],
            sum_all["review"], sum_all["total"]))

if (sum_all["review"] > 0L) {
  stop(sprintf("INVARIANT VIOLATED: %d row(s) still NA on `replicated`. ",
               sum_all["review"]),
       "Every spec must have a non-NA replicated flag after Prompt 3b.",
       call. = FALSE)
}

cat("\nWritten:\n")
for (f in c(
  "01_table3_vs_paper.md", "01_table3_vs_paper.html",
  "01_table2_vs_paper.md", "01_table2_vs_paper.html",
  "01_table3_six_models.html", "01_table2_seven_models.html"
)) {
  cat("  ", file.path(out_dir, f), "\n")
}
