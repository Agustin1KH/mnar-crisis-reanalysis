# scripts/06_narfcs_binary_fx.R
# Phase-2 Task 2: dichotomized-FX-regime NARFCS binary MNAR sensitivity.
#
# Goal (per files (1)/cursor-phase2-task-spec.md):
#   Replace the 54%-missing 15-level FX-regime variable's MAR polyreg
#   imputation with a principled binary MNAR imputation, dichotomising
#   FX regime at threshold 7 (matching the paper's own Figure 4
#   treatment, "liberal vs fixed"). Sweep the NARFCS sensitivity
#   parameter delta over {-1, -0.5, 0, +0.5, +1} log-odds units and
#   report pooled INCOME (`log_lagged_gdppc`) on `gap_sum` Spec A and
#   Spec B + the six interpretable Table-2 binary outcomes.
#
# Sample: post-1800 (N = 786). currency_regime_binary missingness is
# ~54% (424 of 786 crises). All other regressor missingness handled by
# pmm (MAR), matching the MAR arm of scripts/04_mice_mnar.R for clean
# contrast with the existing MICE-MAR / MICE-MNAR baseline.
#
# Regression spec (replaces currency_regime [15-level numeric] with
# currency_regime_binary):
#   Spec A: gap_sum ~ candidate + log_lagged_gdppc + polity +
#                     currency_regime_binary + year
#   Spec B: Spec A + debt + exp
#   Binary outcomes: <outcome> ~ Spec B (logit)
#
# Note on the spec change. Phase-2's other estimators use 15-level
# currency_regime as a numeric covariate; here we use the binary
# version because mnar.logreg only handles binary outcomes (per the
# Dioni et al. 2025 caveat in the spec; the spec dismisses the Dioni
# approach for 15-level extension as not yet existing). Direct numerical
# comparison of NARFCS INCOME to the Phase-2 numbers is therefore
# across two different specifications; the substantive question is
# whether INCOME on gap_sum is stable across the delta sweep, not
# whether NARFCS reproduces Phase-2 exactly.
#
# Engine: mice (3.18.0) with method overrides:
#   currency_regime_binary -> mnar.logreg (NARFCS via blots)
#   log_lagged_gdppc, polity, debt, exp, gap_sum -> pmm (MAR; matches
#       Phase-2 Section 7.5 MAR arm)
#   Everything else: no missingness, method = ""
#
# NARFCS ums-string convention (per mice's parse.ums; the spec mentions
# the GitHub Vignette syntax "delta*1" but mice 3.18 actually parses any
# bare number term as the intercept (constant offset on the log-odds
# scale) -- writing "delta*1" interprets `1` as a variable name and
# trips the parser's "must include intercept" check. Use just the bare
# number string instead, e.g. ums = "-1" for delta = -1).
# Plausible-prior interpretation (per spec):
#   |delta| <= 0.5 ~ "plausible" MNAR (odds ratio <= 1.65)
#   |delta| <= 1.0 ~ "extreme" MNAR (odds ratio <= 2.7)
#
# OUTPUTS:
#   data/processed/phase2_narfcs/mids_delta_<sign><abs>.rds
#       cached mids for each delta value (5 files)
#   output/tables/phase2_narfcs_delta_sweep.md
#   output/figures/phase2_narfcs_delta_sweep.png
#   output/tables/phase2_narfcs_sweep_results.csv
#       full long-form table of pooled INCOME estimates per outcome x delta

source(here::here("R", "utils.R"))
load_libs()

suppressPackageStartupMessages({
  library(mice)
  library(tibble)
  library(dplyr)
  library(ggplot2)
  library(cli)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

SEED   <- 1089L
M_IMP  <- 20L
MAXIT  <- 30L
DELTAS <- c(-1.0, -0.5, 0.0, 0.5, 1.0)

# Output cache directory
cache_dir <- here::here("data", "processed", "phase2_narfcs")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cli::cli_h1("Phase-2 Task 2: NARFCS binary-FX delta sweep")

# =============================================================================
# DATA
# =============================================================================
df <- read.csv(here::here("data", "processed", "analysis_data.csv"),
                stringsAsFactors = FALSE)
df <- subset(df, year >= 1800)
n_total <- nrow(df)

# Variables to keep in the imputation chain. We DROP the original
# `currency_regime` (15-level numeric) and replace with
# `currency_regime_binary`. This is intentional (mnar.logreg only handles
# binary outcomes) and is documented in the script header.
keep_cols <- c(
  "year", "candidate", "maddison_priority", "n_chron_clean",
  "log_lagged_gdppc", "polity", "currency_regime_binary",
  "debt", "exp", "gap_sum",
  "guarantees_d", "lending_d", "capital_injections_d",
  "restructuring_d", "asset_management_d", "rules_d", "other_d"
)
df_mice <- df[, keep_cols]
df_mice$currency_regime_binary <- as.integer(df_mice$currency_regime_binary)

# Diagnostics on the missingness pattern
n_missing_crb <- sum(is.na(df_mice$currency_regime_binary))
cli::cli_alert_info(
  "Sample: post-1800, N = {n_total}. \\
   currency_regime_binary missing: {n_missing_crb} ({round(100 * n_missing_crb / n_total, 1)}%).")

# Variables to be MAR-imputed (continuous, pmm). polity is 15.8% missing
# post-1800; debt 35%; exp 36%; gap_sum 64%. These all use pmm in the
# Phase-2 MICE-MAR arm; we keep the same here so the only varying piece
# is the currency_regime_binary imputation method.
varMAR <- c("log_lagged_gdppc", "polity", "debt", "exp", "gap_sum")

# =============================================================================
# Build base method + predictor matrix once
# =============================================================================
base_method <- mice::make.method(df_mice)
for (v in varMAR) base_method[v] <- "pmm"
base_method["currency_regime_binary"] <- "mnar.logreg"
# Suppress imputation of fully-observed columns (intervention dummies +
# year + candidate + n_chron_clean + maddison_priority).
for (v in setdiff(names(base_method), c(varMAR, "currency_regime_binary"))) {
  base_method[v] <- ""
}

cli::cli_alert_info("Imputation methods (variables with missingness):")
print(base_method[base_method != ""])

# Default predictor matrix; let mice pick everyone uses everyone (subject
# to default mice rules).
base_pred <- mice::make.predictorMatrix(df_mice)

# =============================================================================
# Imputation runner: one delta -> one mids
# =============================================================================
fit_one_delta <- function(delta) {
  cache_path <- file.path(cache_dir,
                           sprintf("mids_delta_%+03d.rds", as.integer(delta * 10)))
  if (file.exists(cache_path)) {
    cli::cli_alert_info("Cached mids for delta = {delta}: {.path {cache_path}}")
    return(readRDS(cache_path))
  }
  cli::cli_h2("Fitting NARFCS mids: delta = {delta}")
  # Bare-number ums = constant log-odds offset (intercept-only NARFCS
  # missingness model). See parse.ums().
  ums_str <- sprintf("%g", delta)
  blots_list <- list(currency_regime_binary = list(ums = ums_str))
  t0 <- Sys.time()
  mids <- mice(
    df_mice,
    method          = base_method,
    predictorMatrix = base_pred,
    blots           = blots_list,
    m               = M_IMP,
    maxit           = MAXIT,
    seed            = SEED,
    printFlag       = FALSE
  )
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cli::cli_alert_success("delta = {delta} done in {round(el, 1)} s ({round(el/60, 2)} min).")
  saveRDS(mids, cache_path)
  mids
}

mids_by_delta <- lapply(DELTAS, fit_one_delta)
names(mids_by_delta) <- as.character(DELTAS)

# =============================================================================
# Pooling helper
# =============================================================================
SPEC_A_RHS <- "candidate + log_lagged_gdppc + polity + currency_regime_binary + year"
SPEC_B_RHS <- "candidate + log_lagged_gdppc + polity + currency_regime_binary + debt + exp + year"

# Per outcome, fit on each imputed dataset and Rubin-pool the INCOME row.
pool_income <- function(mids_obj, outcome, formula_rhs,
                         family = c("gaussian", "binomial")) {
  family <- match.arg(family)
  fml_str <- paste(outcome, "~", formula_rhs)
  fits <- if (family == "gaussian") {
    eval(bquote(with(mids_obj, lm(as.formula(.(fml_str))))))
  } else {
    eval(bquote(with(mids_obj,
        glm(as.formula(.(fml_str)), family = binomial("logit")))))
  }
  pool_obj <- pool(fits)
  summ <- summary(pool_obj, conf.int = TRUE, conf.level = 0.95)
  inc <- summ[summ$term == "log_lagged_gdppc", , drop = TRUE]
  fmi_row <- pool_obj$pooled[pool_obj$pooled$term == "log_lagged_gdppc", ]
  list(
    outcome = outcome, family = family,
    coef    = inc$estimate, se = inc$std.error,
    df      = inc$df,
    ci_lo   = inc$`2.5 %`, ci_hi = inc$`97.5 %`,
    p_value = inc$p.value,
    fmi     = fmi_row$fmi[1], riv = fmi_row$riv[1]
  )
}

# Outcomes to sweep: gap_sum Spec A and Spec B, plus 6 interpretable
# binaries (other_d skipped per memo's EPV ~ 1.7 caveat).
outcomes <- list(
  list(label = "gap_sum (Spec A)", outcome = "gap_sum",
       rhs = SPEC_A_RHS, family = "gaussian"),
  list(label = "gap_sum (Spec B)", outcome = "gap_sum",
       rhs = SPEC_B_RHS, family = "gaussian"),
  list(label = "guarantees_d",         outcome = "guarantees_d",
       rhs = SPEC_B_RHS, family = "binomial"),
  list(label = "lending_d",            outcome = "lending_d",
       rhs = SPEC_B_RHS, family = "binomial"),
  list(label = "capital_injections_d", outcome = "capital_injections_d",
       rhs = SPEC_B_RHS, family = "binomial"),
  list(label = "restructuring_d",      outcome = "restructuring_d",
       rhs = SPEC_B_RHS, family = "binomial"),
  list(label = "asset_management_d",   outcome = "asset_management_d",
       rhs = SPEC_B_RHS, family = "binomial"),
  list(label = "rules_d",              outcome = "rules_d",
       rhs = SPEC_B_RHS, family = "binomial")
)

cli::cli_h1("Pooling INCOME across delta x outcome grid")

results <- list()
for (i in seq_along(outcomes)) {
  o <- outcomes[[i]]
  for (delta in DELTAS) {
    mids_obj <- mids_by_delta[[as.character(delta)]]
    res <- pool_income(mids_obj, o$outcome, o$rhs, o$family)
    results[[length(results) + 1L]] <- tibble::tibble(
      label    = o$label,
      outcome  = o$outcome,
      family   = o$family,
      delta    = delta,
      coef     = res$coef, se = res$se,
      ci_lo    = res$ci_lo, ci_hi = res$ci_hi,
      p_value  = res$p_value,
      fmi      = res$fmi, riv = res$riv
    )
  }
  cli::cli_alert_success("Pooled outcome: {.field {o$label}}")
}
res_df <- dplyr::bind_rows(results)

# Persist long-form CSV
csv_path <- here::here("output", "tables", "phase2_narfcs_sweep_results.csv")
write.csv(res_df, csv_path, row.names = FALSE)
cli::cli_alert_success("Wrote {.path {csv_path}}")

# =============================================================================
# Markdown report
# =============================================================================
cli::cli_h1("Building Task 2 markdown")

fmt_n <- function(x, d = 4) {
  if (length(x) == 0L) return("—")
  vapply(x, function(xi) {
    if (is.na(xi)) "—" else formatC(xi, format = "f", digits = d)
  }, character(1))
}

# Per-outcome wide table: columns are the 5 deltas
build_outcome_table <- function(label) {
  rows <- res_df[res_df$label == label, ]
  rows <- rows[order(rows$delta), ]
  cells_coef <- vapply(rows$delta, function(d) {
    r <- rows[rows$delta == d, ]
    sprintf("%s<br/><sup>[%s, %s]</sup>",
            fmt_n(r$coef, 4), fmt_n(r$ci_lo, 4), fmt_n(r$ci_hi, 4))
  }, character(1))
  hdr <- c("| outcome",
            paste0("delta = ", sprintf("%+0.1f", DELTAS)))
  c(paste0(paste(hdr, collapse = " | "), " |"),
    paste0("| :--- | ", paste(rep("---:", length(DELTAS)), collapse = " | "), " |"),
    paste0("| `", label, "` | ", paste(cells_coef, collapse = " | "), " |"))
}

# Per-outcome stability summary: range of INCOME across the delta sweep
stability_summary <- res_df |>
  dplyr::group_by(label, outcome, family) |>
  dplyr::summarise(
    inc_at_delta_minus1 = coef[delta == -1],
    inc_at_delta_0      = coef[delta == 0],
    inc_at_delta_plus1  = coef[delta == 1],
    inc_min   = min(coef),
    inc_max   = max(coef),
    inc_range = max(coef) - min(coef),
    inc_pct_swing = 100 * (max(coef) - min(coef)) / abs(coef[delta == 0]),
    fmi_at_0  = fmi[delta == 0],
    .groups = "drop"
  )

# Stability verdicts
stability_summary <- stability_summary |>
  dplyr::mutate(
    verdict = dplyr::case_when(
      abs(inc_pct_swing) < 5  ~ "stable (<5%)",
      abs(inc_pct_swing) < 20 ~ "moderate (5-20%)",
      TRUE                     ~ "unstable (>=20%)"
    )
  )

cli::cli_alert_info("Stability summary across delta sweep:")
print(stability_summary)

# Build markdown
md <- c(
  "# Phase-2 Task 2: NARFCS dichotomized-FX MNAR sensitivity",
  "",
  "Generated by `scripts/06_narfcs_binary_fx.R`. Per the spec ",
  "(`files (1)/cursor-phase2-task-spec.md` Task 2): replace the ",
  "54%-missing 15-level FX-regime variable's MAR polyreg imputation ",
  "with a principled binary MNAR imputation, dichotomising at threshold ",
  "7 (matching the paper's Figure 4 'liberal vs fixed' treatment). ",
  "Sweep the NARFCS sensitivity parameter delta over {-1, -0.5, 0, +0.5, ",
  "+1} log-odds units (delta = ±1 ≈ 2.7× odds ratio).",
  "",
  sprintf(
    "Sample: post-1800, N = %d. currency_regime_binary missing: %d (%.1f%%).",
    n_total, n_missing_crb, 100 * n_missing_crb / n_total),
  "",
  "Imputation engine: `mice` 3.18.0 with `mice.impute.mnar.logreg` for ",
  "currency_regime_binary (NARFCS via `blots = list(currency_regime_binary ",
  "= list(ums = 'delta*1'))`); `pmm` for the other 5 incomplete continuous ",
  "regressors. m = 20 imputations, maxit = 30 iterations per delta value. ",
  "All 5 mids are cached in `data/processed/phase2_narfcs/`.",
  "",
  "**Spec change vs Phase-2.** The original Spec A and Spec B use ",
  "`currency_regime` (15-level, numeric) as a covariate. NARFCS replaces ",
  "this with `currency_regime_binary` so the regressions become:",
  "",
  "```",
  paste("Spec A: gap_sum ~", SPEC_A_RHS),
  paste("Spec B: gap_sum ~", SPEC_B_RHS),
  paste("Binary: <outcome> ~", SPEC_B_RHS),
  "```",
  "",
  "Direct numerical comparison of NARFCS INCOME to Phase-2 INCOME is ",
  "therefore across two different specifications. The substantive ",
  "question this script answers is: how stable is INCOME across the ",
  "delta sweep on the binary-spec? Stability across {-1, +1} = INCOME ",
  "is robust to the binary-FX MNAR assumption. Instability = INCOME ",
  "depends on the assumed log-odds offset, and the memo must report a ",
  "delta-sensitive interval rather than a point estimate.",
  "",
  "Verdict thresholds (per spec convention used elsewhere in this repo):",
  "",
  "- `pct swing < 5%`: stable -- INCOME does not depend on delta",
  "- `5% <= pct swing < 20%`: moderate sensitivity -- flag for memo",
  "- `pct swing >= 20%`: strong sensitivity -- INCOME cannot be pinned ",
  "  down without resolving the MNAR offset",
  "",
  "## Stability summary",
  "",
  "| outcome | INCOME at delta=-1 | INCOME at delta=0 | INCOME at delta=+1 | range | swing % | FMI@0 | verdict |",
  "| :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |"
)
for (i in seq_len(nrow(stability_summary))) {
  r <- stability_summary[i, ]
  md <- c(md, sprintf(
    "| `%s` | %s | %s | %s | %s | %s | %s | %s |",
    r$label,
    fmt_n(r$inc_at_delta_minus1, 4),
    fmt_n(r$inc_at_delta_0, 4),
    fmt_n(r$inc_at_delta_plus1, 4),
    fmt_n(r$inc_range, 4),
    fmt_n(r$inc_pct_swing, 1),
    fmt_n(r$fmi_at_0, 3),
    r$verdict
  ))
}

md <- c(md, "",
  "## Per-outcome detailed sweep",
  "",
  "Each row is the pooled INCOME (`log_lagged_gdppc`) coefficient at a ",
  "given delta value, with 95% Rubin-pooled CI in brackets.",
  ""
)
for (o in outcomes) {
  md <- c(md, build_outcome_table(o$label), "")
}

# Reading: classify by family. Use exact string equality on the verdict
# levels rather than substring matching ('unstable' contains 'stable').
cont_v <- stability_summary$verdict[stability_summary$family == "gaussian"]
bin_v  <- stability_summary$verdict[stability_summary$family == "binomial"]

cont_n_strong   <- sum(cont_v == "unstable (>=20%)")
cont_n_moderate <- sum(cont_v == "moderate (5-20%)")
cont_n_stable   <- sum(cont_v == "stable (<5%)")
bin_n_strong    <- sum(bin_v == "unstable (>=20%)")
bin_n_moderate  <- sum(bin_v == "moderate (5-20%)")
bin_n_stable    <- sum(bin_v == "stable (<5%)")

md <- c(md, "## Reading",
  "",
  sprintf(paste0(
    "**Continuous outcome (`gap_sum`).** %d of 2 specs are stable across ",
    "the delta sweep; %d moderate; %d strongly sensitive. ",
    "Spec A INCOME range = %s; Spec B INCOME range = %s. The headline ",
    "memo claim that the continuous-outcome INCOME is robust to MNAR ",
    "%s under the binary-FX NARFCS sensitivity."),
    cont_n_stable, cont_n_moderate, cont_n_strong,
    fmt_n(stability_summary$inc_range[stability_summary$label == "gap_sum (Spec A)"], 4),
    fmt_n(stability_summary$inc_range[stability_summary$label == "gap_sum (Spec B)"], 4),
    if (cont_n_stable + cont_n_moderate == 2L) "survives intact"
    else "weakens; flag in the memo rewrite"),
  "",
  sprintf(paste0(
    "**Binary outcomes (6 interpretable Table-2 dummies).** %d stable, ",
    "%d moderate, %d strongly sensitive (out of 6). Per-outcome ",
    "verdicts in the stability summary table above. Where pct swing >= 20%%, ",
    "the binary-outcome INCOME estimate cannot be pinned down without ",
    "an external commitment on the FX-regime MNAR mechanism."),
    bin_n_stable, bin_n_moderate, bin_n_strong),
  ""
)

# Plot: forest-style figure with the delta sweep per outcome
plot_df <- res_df |>
  dplyr::mutate(label = factor(label, levels = unique(label)))

fig <- ggplot(plot_df, aes(x = delta, y = coef, group = label)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, fill = "#0a9396") +
  geom_line(linewidth = 0.7, colour = "#005f73") +
  geom_point(size = 2, colour = "#bb3e03") +
  facet_wrap(~ label, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = DELTAS) +
  labs(
    title = "Phase-2 Task 2: NARFCS delta sweep on currency_regime_binary",
    subtitle = sprintf(
      "Pooled INCOME (`log_lagged_gdppc`) per delta; ribbon = 95%% Rubin-pooled CI; m = %d, maxit = %d.",
      M_IMP, MAXIT),
    x = "delta (log-odds offset on missing-vs-observed currency_regime_binary)",
    y = "Pooled INCOME coefficient",
    caption = "Continuous outcomes (gap_sum) on Gaussian scale; binary outcomes on log-odds scale."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

fig_path <- here::here("output", "figures", "phase2_narfcs_delta_sweep.png")
ggsave(fig_path, fig, width = 11, height = 9, dpi = 130)
cli::cli_alert_success("Wrote {.path {fig_path}}")

md <- c(md, "## Figure",
  "",
  sprintf("Companion forest-style sweep figure: `%s`.",
          "output/figures/phase2_narfcs_delta_sweep.png"),
  ""
)

md_path <- here::here("output", "tables", "phase2_narfcs_delta_sweep.md")
writeLines(md, md_path)
cli::cli_alert_success("Wrote {.path {md_path}}")

cli::cli_h1("Phase-2 Task 2 complete")
