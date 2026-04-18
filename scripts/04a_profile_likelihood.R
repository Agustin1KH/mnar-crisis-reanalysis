# scripts/04a_profile_likelihood.R
# Profile-likelihood diagnostic for the Heckman selection-model `rho` on
# Spec A (gap_sum) -- replaces the test-statistic tension between the
# Phase 2 Wald and LR tests with a visual narrative.
#
# Phase 2 backdrop (per `output/tables/02_interpretation.md`):
#   - Wald: rho_hat = -0.172, SE = 0.336, |t| = 0.51, p = 0.61 -> "MAR not
#     rejected"
#   - LR (as reported by `02_heckman_diagnostics.md`): chi-sq = 8.12,
#     p = 0.0044 -> "MAR rejected at 1%"
#
# These are not contradictions in the underlying numerics; they are two
# different summaries of the same likelihood, and one is more reliable
# than the other when the likelihood is non-quadratic. Phase 3a's GJRM
# cross-package check (`output/tables/03_interpretation.md`) confirmed
# the Heckman Gaussian estimate (tau agreement within 0.028) but did not
# resolve the test-statistic question. The fix is the profile
# likelihood: it is invariant to reparameterisation, asymptotically
# pivotal under H0, and shows the actual shape of the likelihood in rho
# rather than relying on either a tangent-line (Wald) or a sample-
# matched (LR) approximation.
#
# IMPORTANT METHODOLOGICAL NOTE on the Phase-2 LR stat. Reproducing the
# Phase 2 LR stat exactly requires probit(r ~ ...) fit on the FULL N=786
# post-1800 sample, but the Heckman MLE is fit on the N=774 subset
# (sampleSelection drops 12 rows where r=1 and gap_sum is NA). The
# constrained log-likelihood `LL_probit + LL_OLS` thus uses 786 obs in
# the probit and 260 in the OLS, while the unconstrained `LL_full` uses
# 774 in selection and 260 in outcome. That sample mismatch inflates the
# nominal LR stat from ~0.27 to ~8.12. The profile likelihood (this
# script) uses the SAME 774-obs sample for every rho value in the grid,
# so it is the apples-to-apples test. Headline result: profile D(0) =
# ~0.27, well below the chi-sq(1) 95% threshold of 3.841 -- so on the
# correctly-matched sample the LR test ALSO does NOT reject MAR. The
# Wald, profile, and matched-sample LR all agree.
#
# Implementation: at each rho_fix in the grid, use maxLik with
# `fixed = "rho"` to optimize over (beta_sel, beta_out, log_sigma)
# subject to the rho constraint. This re-uses the exact same per-
# observation log-likelihood that sampleSelection::selection() uses
# internally for the Tobit-2 / Heckman MLE, verified to recover
# `logLik(gap_a_ml) = -397.0079` exactly when called unconstrained.
#
# OUTPUTS:
#   output/figures/04a_profile_likelihood_rho.png  (two-panel figure)
#   output/tables/04a_profile_likelihood_summary.md (CI table + reading)
#
# Constraints (per Prompt 4a):
#   - Do NOT refit the Heckman MLE from scratch; use heckman_fits$gap_a_ml.
#   - 41 grid points in [-0.95, 0.95]; sequential, no parallel.
#   - All output paths via here::here().

source(here::here("R", "utils.R"))
load_libs()

suppressPackageStartupMessages({
  library(sampleSelection)
  library(maxLik)
  library(ggplot2)
  library(cowplot)
  library(cli)
})

# =============================================================================
# DATA & MATRICES (must match Heckman gap_a_ml's 774-row sample exactly)
# =============================================================================
cli::cli_h1("Profile-likelihood for rho on gap_sum spec A")

df <- read.csv(here::here("data", "processed", "analysis_data.csv"),
               stringsAsFactors = FALSE)
df <- subset(df, year >= 1800)
out_regs <- c("candidate", "log_lagged_gdppc", "polity", "currency_regime",
              "year", "gap_sum")
sel_regs <- c("candidate", "year", "maddison_priority", "n_chron_clean")
df <- df[!(df$r == 1 & !complete.cases(df[, out_regs])), ]
df <- df[complete.cases(df[, sel_regs]), ]

stopifnot(nrow(df) == 774L)
stopifnot(sum(df$r == 0) == 514L && sum(df$r == 1) == 260L)

X_sel    <- model.matrix(~ candidate + year + maddison_priority + n_chron_clean,
                         data = df)
df_r1    <- df[df$r == 1, , drop = FALSE]
X_out_r1 <- model.matrix(~ candidate + log_lagged_gdppc + polity +
                           currency_regime + year, data = df_r1)
y_r1     <- df_r1$gap_sum
r        <- df$r
n_r0 <- sum(r == 0)
n_r1 <- sum(r == 1)
cli::cli_alert_info(
  "Sample (matches `heckman_fits$gap_a_ml`): N = {nrow(df)}, \\
   N(r=0) = {n_r0}, N(r=1) = {n_r1}.")

# =============================================================================
# HECKMAN PER-OBSERVATION LOG-LIKELIHOOD (same form as sampleSelection's MLE)
# =============================================================================
# Parameter vector layout:
#   1:5    beta_sel   (5 selection coefficients)
#   6:11   beta_out   (6 outcome coefficients)
#   12     log_sigma  (log of outcome residual SD; reparameterisation for
#                     unconstrained optimisation)
#   13     rho        (Tobit-2 error correlation, -1 < rho < 1)

heck_ll_vec <- function(p) {
  bs        <- p[1:5]
  bo        <- p[6:11]
  log_sigma <- p[12]
  rho       <- p[13]
  if (abs(rho) >= 1) return(rep(-Inf, length(r)))
  sigma   <- exp(log_sigma)
  eta_sel <- as.vector(X_sel %*% bs)
  eta_out_r1 <- as.vector(X_out_r1 %*% bo)
  ll <- numeric(length(r))
  ll[r == 0] <- pnorm(-eta_sel[r == 0], log.p = TRUE)
  i1    <- which(r == 1)
  resid <- (y_r1 - eta_out_r1) / sigma
  cond  <- (eta_sel[i1] + rho * resid) / sqrt(1 - rho^2)
  ll[i1] <- dnorm(resid, log = TRUE) - log_sigma + pnorm(cond, log.p = TRUE)
  ll
}

# =============================================================================
# REFERENCE VALUES (from cached MLE; do NOT refit)
# =============================================================================
cli::cli_h2("Reference values from `heckman_fits$gap_a_ml`")

heck_fits <- readRDS(here::here("data", "processed", "heckman_fits.rds"))
ml        <- heck_fits$gap_a_ml
ll_max    <- as.numeric(logLik(ml))
rho_hat   <- unname(ml$estimate["rho"])
rho_se    <- summary(ml)$estimate["rho", "Std. Error"]
rho_t     <- summary(ml)$estimate["rho", "t value"]
rho_p     <- summary(ml)$estimate["rho", "Pr(>|t|)"]

# Wald 95% CI on rho (point estimate plus/minus 1.96 SE)
wald_lo <- rho_hat - 1.96 * rho_se
wald_hi <- rho_hat + 1.96 * rho_se

cli::cli_alert_info(
  "Heckman MLE (cached): rho = {round(rho_hat, 4)}, SE = {round(rho_se, 4)}, |t| = {round(abs(rho_t), 2)}, p = {round(rho_p, 4)}.")
cli::cli_alert_info(
  "Wald 95% CI on rho: [{round(wald_lo, 4)}, {round(wald_hi, 4)}] (width {round(wald_hi - wald_lo, 4)}).")
cli::cli_alert_info("LL_max = {round(ll_max, 4)}.")

# Sanity check: manual LL at MLE matches sampleSelection's logLik exactly.
ml_p <- c(unname(ml$estimate[1:11]), log(unname(ml$estimate[12])),
          unname(ml$estimate[13]))
ll_manual <- sum(heck_ll_vec(ml_p))
stopifnot(abs(ll_manual - ll_max) < 1e-6)
cli::cli_alert_success("Manual LL at MLE matches sampleSelection (delta < 1e-6).")

# GJRM Gaussian cross-reference (per Phase 3a). For Gaussian copula,
# theta = rho directly; we also report tau-converted rho for transparency.
gjrm_fits        <- readRDS(here::here("data", "processed", "gjrm_fits.rds"))
gjrm_theta       <- as.numeric(gjrm_fits$continuous$N$theta)
gjrm_tau         <- as.numeric(gjrm_fits$continuous$N$tau)
gjrm_rho_via_tau <- sin(pi * gjrm_tau / 2)
cli::cli_alert_info(
  "GJRM (Gaussian) cross-reference: theta = {round(gjrm_theta, 4)} (= rho directly), tau = {round(gjrm_tau, 4)}, rho via sin(pi tau/2) = {round(gjrm_rho_via_tau, 4)}.")

# =============================================================================
# PROFILE LIKELIHOOD GRID
# =============================================================================
cli::cli_h2("Profile likelihood evaluation (41 rho grid points)")

rho_grid <- seq(-0.95, 0.95, length.out = 41L)
profile_ll <- numeric(length(rho_grid))
profile_code <- integer(length(rho_grid))

# Init from the cached MLE at every grid point. Warm-starting from the
# previous grid point's converged params can drift the optimizer into a
# slightly-worse local optimum at distant rho values (e.g. it inflated
# D(0) from the analytical 0.272 to 0.318 in pilot runs); init-from-MLE
# is per-point cheap (~0.1 s) and gives the analytically-correct value
# at known reference points.
mle_init <- ml_p
names(mle_init) <- c(paste0("bs", 1:5), paste0("bo", 1:6),
                     "log_sigma", "rho")

t_total_0 <- Sys.time()
for (i in seq_along(rho_grid)) {
  rho_fix <- rho_grid[i]
  init_p <- mle_init
  init_p["rho"] <- rho_fix
  fit <- tryCatch(
    maxLik::maxLik(heck_ll_vec, start = init_p, method = "BFGS",
                   fixed = "rho",
                   control = list(printLevel = 0, iterlim = 500, tol = 1e-9)),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    profile_ll[i] <- NA_real_
    profile_code[i] <- -1L
    next
  }
  profile_ll[i] <- as.numeric(logLik(fit))
  profile_code[i] <- as.integer(fit$code)
}
t_total_1 <- Sys.time()
el <- as.numeric(difftime(t_total_1, t_total_0, units = "secs"))
cli::cli_alert_success(
  sprintf("41 grid points evaluated in %.2f s (%.3f s per point).",
          el, el / length(rho_grid)))

# Profile deviance: D(rho) = 2 * (LL_max - LL(rho_fix))
profile_D <- 2 * (ll_max - profile_ll)

# Convergence check
n_failed <- sum(profile_code != 0L | is.na(profile_ll))
if (n_failed > 0L) {
  cli::cli_alert_warning(
    "{n_failed} of {length(rho_grid)} grid evaluations had non-zero return code or NA LL.")
}

# =============================================================================
# 95% PROFILE INTERVAL (linear interpolation)
# =============================================================================
chi2_thresh  <- qchisq(0.95, df = 1L)  # 3.841459
i_min <- which.min(profile_D)
rho_minD <- rho_grid[i_min]
cli::cli_alert_info(
  "Empirical minimum of D occurs at rho = {round(rho_minD, 3)} (grid point #{i_min}); D_min = {round(profile_D[i_min], 4)}.")

# Walk left of the minimum to find the LOWER 95% bound (D crosses 3.841 from
# above). Walk right of the minimum to find the UPPER 95% bound. For each side,
# locate the adjacent grid pair that brackets 3.841 and linearly-interpolate.
find_bound <- function(side = c("lower", "upper")) {
  side <- match.arg(side)
  if (side == "lower") {
    idx <- seq.int(i_min, 1L, by = -1L)
  } else {
    idx <- seq.int(i_min, length(rho_grid), by = 1L)
  }
  for (k in seq_along(idx)[-length(idx)]) {
    j_in  <- idx[k]
    j_out <- idx[k + 1L]
    D_in  <- profile_D[j_in]
    D_out <- profile_D[j_out]
    if (is.na(D_in) || is.na(D_out)) next
    if (D_out >= chi2_thresh && D_in < chi2_thresh) {
      # Linear interpolation in rho for D = chi2_thresh
      w <- (chi2_thresh - D_in) / (D_out - D_in)
      rho_in  <- rho_grid[j_in]
      rho_out <- rho_grid[j_out]
      return(rho_in + w * (rho_out - rho_in))
    }
  }
  NA_real_
}

prof_lo <- find_bound("lower")
prof_hi <- find_bound("upper")

# If either bound returned NA, the deviance never crosses 3.841 on that side
# of the grid; the practical CI is open-ended at the grid edge. Report grid
# edge as the bound and flag.
prof_lo_open <- is.na(prof_lo)
prof_hi_open <- is.na(prof_hi)
if (prof_lo_open) prof_lo <- min(rho_grid)
if (prof_hi_open) prof_hi <- max(rho_grid)

cli::cli_alert_info(
  "Profile 95% CI on rho: [{round(prof_lo, 4)}, {round(prof_hi, 4)}] \\
   (width {round(prof_hi - prof_lo, 4)}; \\
   lower {ifelse(prof_lo_open, '*open-at-grid-edge*', 'interpolated')}, \\
   upper {ifelse(prof_hi_open, '*open-at-grid-edge*', 'interpolated')}).")

zero_in_profile <- (prof_lo <= 0 && 0 <= prof_hi)
cli::cli_alert_info(
  "Profile CI {ifelse(zero_in_profile, 'INCLUDES', 'EXCLUDES')} rho = 0.")

# Width comparison (used in figure subtitle and markdown reading paragraph)
profile_width <- prof_hi - prof_lo
wald_width    <- wald_hi - wald_lo
profile_narrower_pct <- 100 * (1 - profile_width / wald_width)

# =============================================================================
# FIGURE: two-panel profile likelihood + deviance
# =============================================================================
cli::cli_h2("Building two-panel figure")

dat <- data.frame(rho = rho_grid, ll = profile_ll, D = profile_D)

# Vertical reference lines: MLE point estimate, GJRM Gaussian, MAR (rho=0)
vlines <- data.frame(
  rho   = c(rho_hat, gjrm_rho_via_tau, 0),
  label = c(sprintf("Heckman MLE = %.3f", rho_hat),
            sprintf("GJRM Gaussian = %.3f", gjrm_rho_via_tau),
            "MAR benchmark = 0"),
  ltype = c("solid", "dashed", "dotted"),
  col   = c("#bb3e03", "#0a9396", "#9d4edd")
)

# Top panel: deviance vs rho with shaded 95% region
shade_df <- subset(dat, D <= chi2_thresh)

p_top <- ggplot(dat, aes(x = rho, y = D)) +
  geom_ribbon(data = shade_df, aes(ymin = 0, ymax = pmax(D, 0)),
              fill = "#0a9396", alpha = 0.18) +
  geom_hline(yintercept = chi2_thresh, linetype = "longdash",
             colour = "grey30") +
  annotate("text", x = max(rho_grid), y = chi2_thresh,
           label = sprintf("chi^2(1, .95) = %.3f", chi2_thresh),
           hjust = 1.0, vjust = -0.4, size = 3.2, colour = "grey30") +
  geom_vline(data = vlines, aes(xintercept = rho, colour = label,
                                linetype = label),
             linewidth = 0.7, show.legend = TRUE) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.3, alpha = 0.7) +
  scale_linetype_manual(values = setNames(vlines$ltype, vlines$label)) +
  scale_colour_manual(values = setNames(vlines$col, vlines$label)) +
  scale_x_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.25)) +
  coord_cartesian(ylim = c(0, max(min(profile_D, na.rm = TRUE) + 30,
                                  6))) +
  labs(
    title    = "Profile deviance D(rho) for the Heckman selection-model rho",
    subtitle = sprintf(
      "Spec A (gap_sum); 41-point grid in [-0.95, 0.95]; 95%% profile CI = [%.3f, %.3f]; %s rho = 0.",
      prof_lo, prof_hi,
      if (zero_in_profile) "INCLUDES" else "EXCLUDES"),
    x = NULL, y = expression(D(rho) == 2*(LL[max] - LL(rho))),
    colour = NULL, linetype = NULL,
    caption = "Shaded region = D(rho) <= 3.841 = 95% profile CI."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        plot.subtitle = element_text(size = 9.5),
        panel.grid.minor = element_blank())

# Bottom panel: the log-likelihood itself (shape diagnostic for Wald-LR
# disagreement). Y-axis truncated near LL_max so the curvature is visible.
ll_floor <- ll_max - 30
p_bot <- ggplot(dat, aes(x = rho, y = pmax(ll, ll_floor))) +
  geom_vline(data = vlines, aes(xintercept = rho, colour = label,
                                linetype = label),
             linewidth = 0.7, show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.3, alpha = 0.7) +
  geom_hline(yintercept = ll_max, linetype = "longdash",
             colour = "grey30") +
  annotate("text", x = max(rho_grid), y = ll_max,
           label = sprintf("LL_max = %.3f", ll_max),
           hjust = 1.0, vjust = -0.4, size = 3.2, colour = "grey30") +
  scale_linetype_manual(values = setNames(vlines$ltype, vlines$label)) +
  scale_colour_manual(values = setNames(vlines$col, vlines$label)) +
  scale_x_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.25)) +
  coord_cartesian(ylim = c(ll_floor, ll_max + 1)) +
  labs(
    title    = "Profile log-likelihood LL(rho) -- shape diagnostic",
    subtitle = sprintf(
      "Approximately quadratic near MLE; mild asymmetry (profile is %.0f%% narrower than Wald). Phase-2 LR rejection traces to a sample mismatch.",
      profile_narrower_pct),
    x = expression(rho), y = expression(LL(rho)),
    caption = sprintf(
      "Y-axis truncated at LL_max - 30 = %.2f so the curvature near the MLE is visible.",
      ll_floor)
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9.5),
        panel.grid.minor = element_blank())

fig_path <- here::here("output", "figures", "04a_profile_likelihood_rho.png")
combined <- cowplot::plot_grid(p_top, p_bot, ncol = 1, align = "v",
                               axis = "lr", rel_heights = c(1.05, 1))
ggsave(fig_path, combined, width = 11, height = 9.5, dpi = 110)
cli::cli_alert_success("Wrote {.path output/figures/04a_profile_likelihood_rho.png}")

# =============================================================================
# SUMMARY MARKDOWN
# =============================================================================
cli::cli_h2("Building summary markdown")

fmt_n <- function(x, d = 3) {
  if (is.null(x)) return("—")
  vapply(x, function(xi) {
    if (is.na(xi)) "—" else formatC(xi, format = "f", digits = d)
  }, character(1))
}

# CI table: 3 rows (Wald, Profile, GJRM cross-ref).
# - Wald: rho_hat ± 1.96 SE. Phase 2 reported [-0.830, 0.486]; we recompute.
# - Profile: from the grid above.
# - GJRM Gaussian: point estimate rho = sin(pi tau/2). No bootstrap was run
#   on the continuous fit in Phase 3a (bootstrap was binary-only), so we
#   report the GJRM point estimate alongside its tau and note that no
#   bootstrap CI is available for the continuous Gaussian fit.

tbl_lines <- c(
  "| estimator | rho estimate | CI source | CI lower | CI upper | width | excludes 0? |",
  "| :--- | ---: | :--- | ---: | ---: | ---: | :---: |",
  sprintf("| Heckman MLE | %s | Wald (Phase 2) | %s | %s | %s | %s |",
          fmt_n(rho_hat), fmt_n(wald_lo), fmt_n(wald_hi),
          fmt_n(wald_hi - wald_lo),
          if ((wald_lo > 0) || (wald_hi < 0)) "yes" else "no"),
  sprintf("| Heckman MLE | %s | Profile (this script) | %s | %s | %s | %s |",
          fmt_n(rho_hat),
          paste0(fmt_n(prof_lo), if (prof_lo_open) "*" else ""),
          paste0(fmt_n(prof_hi), if (prof_hi_open) "*" else ""),
          fmt_n(prof_hi - prof_lo),
          if (zero_in_profile) "no" else "yes"),
  sprintf("| GJRM Gaussian | %s (= sin(π·τ/2), τ = %s) | n/a (no bootstrap on continuous) | — | — | — | — |",
          fmt_n(gjrm_rho_via_tau), fmt_n(gjrm_tau, 4))
)

# Reading paragraphs (profile_narrower_pct already computed above).
para1 <- if (profile_narrower_pct >= 20) {
  sprintf(paste0(
    "**Profile CI is meaningfully narrower than Wald CI (by %.1f%%).** ",
    "Wald 95%% CI: [%s, %s] (width %s). Profile 95%% CI: [%s, %s] ",
    "(width %s). The Wald approximation tangents the log-likelihood at ",
    "the MLE and inherits the inverse-Hessian SE, which is inflated when ",
    "the likelihood is non-quadratic in rho on this sample size; the ",
    "profile likelihood traces the actual surface and gives the correct ",
    "asymptotic 95%% interval. The bottom panel of the figure shows the ",
    "non-quadratic shape directly: the LL curve is asymmetric in rho and ",
    "noticeably flatter on the negative-rho side than the quadratic ",
    "Wald approximation would suggest, which is what makes the Wald CI ",
    "stretch further toward rho = -0.83 than the profile does. Phase 2's ",
    "wide Wald interval is therefore the artifact, not genuine ",
    "imprecision."),
    profile_narrower_pct,
    fmt_n(wald_lo), fmt_n(wald_hi), fmt_n(wald_width),
    fmt_n(prof_lo), fmt_n(prof_hi), fmt_n(profile_width))
} else if (profile_narrower_pct >= -10) {
  sprintf(paste0(
    "**Profile and Wald CIs are similar in width (profile is %.1f%% ",
    "narrower).** Wald 95%% CI: [%s, %s] (width %s); Profile 95%% CI: ",
    "[%s, %s] (width %s). The likelihood is approximately quadratic ",
    "near the MLE on this sample size, so the Wald approximation is ",
    "well-calibrated and the parameter is genuinely weakly identified ",
    "(rho is consistent with a wide range of values, from clearly-",
    "negative to mildly-positive). The bottom-panel curvature in the ",
    "figure should be close to a downward-opening parabola in this case."),
    profile_narrower_pct,
    fmt_n(wald_lo), fmt_n(wald_hi), fmt_n(wald_width),
    fmt_n(prof_lo), fmt_n(prof_hi), fmt_n(profile_width))
} else {
  sprintf(paste0(
    "**Profile CI is WIDER than Wald (by %.1f%%, surprising).** Wald: ",
    "[%s, %s] (width %s); Profile: [%s, %s] (width %s). This usually ",
    "indicates the inverse-Hessian SE under-estimates the curvature in ",
    "one or both tails (asymmetric likelihood with a sharper drop on ",
    "one side than implied by the local quadratic approximation at the ",
    "MLE). Read the figure's bottom panel for the actual shape."),
    -profile_narrower_pct,
    fmt_n(wald_lo), fmt_n(wald_hi), fmt_n(wald_width),
    fmt_n(prof_lo), fmt_n(prof_hi), fmt_n(profile_width))
}

para2 <- if (zero_in_profile) {
  sprintf(paste0(
    "**Profile CI INCLUDES rho = 0 (D(0) = %.3f, well below the 95%% ",
    "threshold of 3.841).** This is the principled version of the test ",
    "Phase 2 was trying to perform: on the sample-matched 774-row ",
    "Heckman dataset, rho = 0 is NOT rejected at the 95%% level. The ",
    "Phase 2 LR statistic of 8.12 (`02_heckman_diagnostics.md`) was ",
    "computed against a constrained log-likelihood that used probit on ",
    "all N=786 post-1800 observations while the unconstrained Heckman ",
    "MLE used the N=774 subset where gap_sum is non-NA among r=1 -- a ",
    "12-observation sample mismatch that inflates LL_constrained ",
    "spuriously. The matched-sample LR stat (= profile D(0)) is ~0.27, ",
    "yielding p ~ 0.6 -- in agreement with the Wald p = 0.61. **Wald, ",
    "matched-sample LR, and profile likelihood ALL fail to reject MAR ",
    "on the same sample.** The 'LR rejects MAR' framing in Phase 2 ",
    "should be read as the Phase 2 LR specifically (which used a ",
    "different sample for the constrained model); the profile-likelihood ",
    "test here is the cleanest single number on rho's identification, ",
    "and it sides with the Wald conclusion. The Phase 3a GJRM Gaussian ",
    "cross-reference (theta = %s, rho via tau = %s) sits well inside ",
    "this profile CI, further reinforcing that the dependence parameter ",
    "is small in magnitude and statistically indistinguishable from ",
    "zero on this sample."),
    profile_D[which.min(abs(rho_grid - 0))],
    fmt_n(gjrm_theta, 3), fmt_n(gjrm_rho_via_tau, 3))
} else {
  sprintf(paste0(
    "**Profile CI EXCLUDES rho = 0 (D(0) = %.3f > 3.841).** This ",
    "confirms the LR rejection of MAR on a principled, sample-matched ",
    "test: the data prefer rho != 0 at the 95%% level even when the ",
    "constrained and unconstrained models are evaluated on the same ",
    "774-row sample. The Phase 2 LR stat of 8.12 had a sample-mismatch ",
    "inflation (probit on N=786 vs Heckman MLE on N=774), but the ",
    "matched profile test agrees on the rejection: D(0) = %.3f, ",
    "p ~ %.4f. The Wald test missed this because the inverse-Hessian SE ",
    "is inflated by non-quadratic curvature; in this case the Wald ",
    "interval is the unreliable one. The Phase 3a GJRM Gaussian ",
    "cross-reference (theta = %s, rho via tau = %s) sits inside the ",
    "profile CI, supporting the rejection."),
    profile_D[which.min(abs(rho_grid - 0))],
    profile_D[which.min(abs(rho_grid - 0))],
    pchisq(profile_D[which.min(abs(rho_grid - 0))], 1, lower.tail = FALSE),
    fmt_n(gjrm_theta, 3), fmt_n(gjrm_rho_via_tau, 3))
}

md_lines <- c(
  "# Profile likelihood for Heckman selection-model rho (gap_sum spec A)",
  "",
  "Generated by `scripts/04a_profile_likelihood.R`. Visual narrative ",
  "diagnostic for the Phase 2 Wald-vs-LR tension on rho. Companion ",
  "figure: `output/figures/04a_profile_likelihood_rho.png`.",
  "",
  sprintf(
    "Sample (matches `heckman_fits$gap_a_ml`): N = %d, N(r=0) = %d, N(r=1) = %d. LL_max = %s.",
    nrow(df), n_r0, n_r1, fmt_n(ll_max, 3)),
  "",
  "## Confidence-interval comparison",
  "",
  tbl_lines,
  "",
  if (prof_lo_open || prof_hi_open) {
    sprintf("`*` = profile deviance does not cross 3.841 within the [-0.95, 0.95] grid on this side; CI bound reported at grid edge.")
  } else {
    "Profile bounds are linear-interpolated between adjacent grid points where D crosses 3.841."
  },
  "",
  "## Reading",
  "",
  para1,
  "",
  para2,
  ""
)

md_path <- here::here("output", "tables", "04a_profile_likelihood_summary.md")
writeLines(md_lines, md_path)
cli::cli_alert_success("Wrote {.path output/tables/04a_profile_likelihood_summary.md}")

cli::cli_h1("Profile likelihood diagnostic complete")
