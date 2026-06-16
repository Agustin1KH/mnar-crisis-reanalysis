# Kim & Lee (2025): identifying-assumption note for the Phase-2 application

Standalone methods note. Written **before** the Kim–Lee implementation
in `scripts/07_kimlee_noexcl.R` and the resulting INCOME estimates in
`output/tables/phase2_kimlee_estimates.md`, so the reader has the
identifying-assumption story in hand before reading the numbers.

**Reference.** Kim, D. & Lee, Y. J. (2025). "Point-Identifying
Semiparametric Sample Selection Models with No Excluded Variable."
arXiv:2502.05353v1 (submitted 7 Feb 2025; downloaded for this note from
https://arxiv.org/abs/2502.05353 on 2026-04-18). Section 2
("Identification") and Section 3 ("Estimators and asymptotic
properties") are the relevant sections; specifically, the linear-index
failure case is documented in Section 2.2 (in v1, around eq. (4) /
discussion at the end of the section), and the empirical nonlinearity
test is proposed in Section 3 (Remark in the practical-implementation
subsection, after the two-step procedure is introduced).

## What Kim & Lee (2025) requires

Their identification result for the latent-outcome slope `β₀` requires
**at least one continuously distributed covariate** in the regressor
set X **and certain nonlinearity in the conditional selection
probability** `p₀(X) := E[D | X]`. The nonlinearity is what does the
identifying work in place of an exclusion restriction: when `p₀(·)` is
nonlinear in X, the partial derivatives `∂p₀/∂Xₖ` and `∂p₀/∂Xⱼ` are
not proportional in general, which lets the sieve-based partial linear
regression in Step 2 separate the identifying variation in X from the
selection-correction term `λ(p̂)`. Their Proposition 1 (single-X case)
and the Multivariate-X corollary (Proposition 2 in v1) make this
formal.

**The case where Kim–Lee identification fails.** The paper explicitly
states that when `p₀(X) = F_ε(X·γ)` (a single-index probit with a
linear index), the partial-derivative ratio collapses to
`(∂p₀/∂xₖ) / (∂p₀/∂xⱼ) = γₖ/γⱼ`, a constant in X. In that case the
Kim–Lee partial-derivative argument provides no separating variation
beyond what a standard linear-index probit already gives. The
linear-index probit is **the canonical Heckman selection-equation
form**, so the paper is explicitly noting that Kim–Lee identification
does not improve on Heckman in the strict-linear-probit case.

## What our selection equation looks like (Kim–Lee version)

Per the task spec, the **Kim–Lee selection equation drops the chronology
count entirely** (that's the whole point of the no-exclusion estimator):

```
r ~ candidate + year + maddison_priority      [Kim & Lee, this script]
```

vs. the Phase-2 Heckman selection equation, which keeps the shadow:

```
r ~ candidate + year + maddison_priority + SHADOW_VAR    [Heckman et al.]
```

Either way, the relevant regressor set for Kim & Lee's identification
condition is the X-set in the OUTCOME equation; the selection equation
just needs to deliver a `p₀(X)` that is nonlinear in X. Three
regressors here:

| regressor | type | values |
| :--- | :--- | :--- |
| `candidate` | binary | {0, 1} |
| `maddison_priority` | binary | {0, 1} |
| `year` | numeric, integer-typed in source | range 1800-2019 (220 unique values on the post-1800 sample) |

Two are discrete and binary. The only "continuous" regressor is
`year`, which is integer-typed in the source data and has 220 unique
values across the 786-row sample. By the convention in Kim & Lee's
setup (`X` includes "at least one continuously distributed covariate"),
`year` *is* a continuously distributed covariate in the practical
sense — it has many distinct values relative to the sample size, no
point mass at any single value above 5%, and its conditional density
on the selection event is well-approximated by a smooth function. So
Kim & Lee's *necessary* identifying condition (point-identification
requires ≥ 1 continuous X-element) is technically satisfied, but only
barely (one continuous element out of three regressors).

The bigger problem is **nonlinearity**. Our selection equation (in
either the Heckman or the Kim–Lee variant) enters all regressors
**linearly through a single index** — exactly the `p₀(X) = F_ε(X·γ)`
failure case that Kim & Lee explicitly call out. Without explicitly
nonlinearizing the selection equation (e.g., by adding a polynomial in
`year` or a smooth via `mgcv::s(year)`), the conditional selection
probability has no nonlinearity beyond what the probit link function
provides on a single linear index. The probit link is monotone and
injective in its argument, so it does not by itself break the
partial-derivative proportionality.

We have decided **not to artificially add a smooth on year to the
selection equation for the Kim–Lee application** (per the user's
review-followup decision recorded in the Phase-2 run history).
Rationale: artificially nonlinearizing the selection equation to force
Kim–Lee identification would be tuning a model to produce the result
Kim–Lee can deliver, and a reviewer would rightly object. We use the
selection equation as-is and report the Kim–Lee numbers honestly with
this caveat front-and-centre.

## Empirical nonlinearity test (per Section 3 of Kim & Lee)

Kim & Lee propose a direct empirical test for the requisite
nonlinearity (Section 3 Remark): fit the Step-1 selection model as a
sieve probit with polynomial-basis expansions of X, then test the
high-order polynomial coefficients for joint significance. A rejection
of the joint zero hypothesis confirms nonlinearity.

In `scripts/07_kimlee_noexcl.R` we run a small version of this test
before reporting the Kim–Lee point estimates: fit the Step-1 selection
probit with a 5-df B-spline expansion of `year` (the only candidate
continuous variable), and report the joint Wald-style chi-square test
for the four high-order spline coefficients (everything beyond the
linear term). A failure to reject is direct evidence that our selection
equation does not satisfy the Kim–Lee identifying nonlinearity, and
the Kim–Lee point estimates should be interpreted accordingly. A
rejection reverses that reading: there *is* enough year-curvature in
selection to give Kim–Lee separating identification, and the resulting
`β̂_INCOME` is more interpretable as a no-exclusion-restriction
robustness check.

## Conclusion: what Kim & Lee provides for our application

**Kim & Lee (2025) functions here as a cross-check with partially-
overlapping identifying assumptions, not as a clean no-exclusion
replacement for Heckman.** Specifically:

- The paper's headline contribution is point identification *without
  an exclusion restriction*. We do not need that contribution for the
  continuous outcome (`gap_sum`), where the exclusion restriction is
  in-sample-defensible across all three Phase-2 shadows and the 03a
  falsification gives `|Δ INCOME|` between 0.0029 and 0.0096 (well
  under the 0.030 stability threshold — see
  `output/tables/phase2_falsification_trajectory.md`).
- We *do* need it for the binary outcomes, where 03b reports
  in-sample violations on 4 of 7 outcomes under every shadow (see
  `output/tables/phase2_exclusion_across_shadows.md`).
- But Kim & Lee's identification mechanism — partial-derivative
  separation between continuous X-elements in a nonlinear selection
  process — does not directly apply to our setup, because our
  selection equation is essentially linear-index probit on a discrete +
  one roughly-continuous-but-noisy regressor set. The empirical
  nonlinearity test in `scripts/07_kimlee_noexcl.R` will tell us
  whether the year regressor carries enough curvature to provide the
  separating variation; in the failure case, Kim & Lee gives an
  estimator whose identification rests on the same linear-index
  bivariate-normal-style assumptions as Heckman, just with a different
  parametrization in Step 2.

**Reading rules for the resulting Kim & Lee numbers in
`output/tables/phase2_kimlee_estimates.md`:**

1. If the nonlinearity test rejects (year-spline high-order coefs
   jointly significant) AND the Kim–Lee `β̂_INCOME` agrees with
   Heckman ML to within ±0.02 on `gap_sum`: this is the
   robustness-strengthens reading. The continuous-outcome INCOME is
   identified the same way under two different identifying-assumption
   sets, which is a cleaner robustness story than Heckman alone.

2. If the nonlinearity test rejects AND Kim–Lee disagrees with Heckman
   by more than ±0.05 on a binary outcome: this is the
   reportable-finding case the spec flagged. The exclusion-restriction-
   driven identification in Heckman was doing load-bearing work that
   the no-exclusion estimator rejects, and the binary-outcome
   conclusions in the paper need a delta-aware interpretation.

3. If the nonlinearity test fails to reject (no curvature in year-on-
   selection): Kim–Lee identification is not formally satisfied, the
   resulting point estimates are best read as "what a partial-linear-
   regression with a sieve-corrected selection-propensity term would
   say", and they cannot be cited as a no-exclusion robustness check
   in the strict sense. They remain a useful directional diagnostic
   (does the Heckman INCOME survive a different parametric correction
   of the selection bias?), but the no-exclusion claim must be
   qualified.

In any of the three readings, the Kim–Lee number is one more
parametric estimator on top of Heckman / GJRM / MICE-MNAR, not a
no-assumption ground truth. The substantive memo conclusions remain
anchored to the cross-method agreement we already have on the
continuous outcome and the cross-method disagreement we already have
on the binary outcomes.

## What this means for the memo rewrite

The Kim & Lee (2025) row in the Table 2 / Table 3 robustness chain
should be added with an asterisk pointing to this note. The asterisk
text should be roughly:

> Kim & Lee (2025) is implemented per the identifying-assumption note
> at `output/tables/07_kimlee_identifying_assumption_note.md`. The
> selection equation is essentially linear-index probit, which is the
> failure case for the paper's no-exclusion identification result. We
> report Kim–Lee numbers as a directional cross-check whose
> identifying assumptions partially overlap with Heckman's; agreement
> strengthens the robustness story, disagreement is a finding to
> investigate but not a clean refutation of either approach.

The memo's "load-bearing assumptions, named explicitly" paragraph
(Section 5 in the current draft) needs the Kim–Lee identifying
assumption added to the named list, with the linear-index-probit
caveat explicit. The honest framing is: Heckman uses bivariate
normality plus an exclusion restriction; Kim & Lee uses
nonlinearity-of-selection plus at least one continuous X-regressor;
under our specification both assumption sets are partially defensible
and partially questionable, so cross-method agreement is what the
robustness story rests on, not the validity of any single method's
identifying assumption.
