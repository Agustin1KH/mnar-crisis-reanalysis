# MNAR Re-analysis of Metrick–Schmelzing (2024)

Bayesian / selection-model re-estimation of the crisis-response regressions in
Metrick & Schmelzing, "Banking-Crisis Interventions Across Time and Space"
(NBER w29281, rev. Feb 2024).

The paper's Table 2 and Table 3 regressions use only 260–273 of 910 crises —
the complete-case subsample after dropping any crisis with a missing regressor.
That subsample has **median year 1996** vs **1904 for the dropped crises**, so
the paper's quantitative claims effectively speak to post-Bretton-Woods data
despite a narrative spanning 33 AD to 2019. The missingness is not random:
pre-modern and smaller-country crises are systematically under-observed.
Complete-case estimation under MNAR is biased.

This repo:

1. Replicates the paper's Tables 2, 3, and Figure 5 cleanly.
2. Re-estimates the headline regressions under a selection model
   (Heckman / Type II Tobit) using established R libraries.
3. Stress-tests the result with copula alternatives and MNAR multiple imputation.
4. Reports sensitivity across assumed error correlation ρ and exclusion restrictions.

## Structure

```
.
├── R/                   helper functions
├── scripts/             sequential analysis scripts (run in order)
│   ├── 00_etl.R         build processed/analysis_data.csv from raw kit
│   ├── 01_replicate_table3.R    match paper's OLS / logit numbers
│   ├── 02_heckman_baseline.R    sampleSelection Heckman (frequentist MLE)
│   ├── 03_gjrm_copula.R         GJRM robustness over copulas
│   ├── 04_mice_mnar.R           Heckman-based MI with miceMNAR
│   └── 05_sensitivity.R         ρ-grid and exclusion-restriction swaps
├── data/
│   ├── raw/             drop original replication kit files here (gitignored)
│   └── processed/       clean CSVs produced by 00_etl.R
├── output/
│   ├── tables/          rendered tables (modelsummary)
│   └── figures/         ggplot figures
└── memo/                LaTeX memo + bibliography
```

## Setup

```r
install.packages(c("renv", "here"))
renv::init()   # first time only
renv::restore()  # subsequent pulls
```

Core packages:

- **readxl, dplyr, tidyr, stringr** — ETL
- **fixest, sandwich, lmtest** — OLS / logit with clustered and robust SEs
- **sampleSelection** — Heckman Type II Tobit (Toomet & Henningsen 2008)
- **GJRM** — generalized joint regression models with copula flexibility
- **miceMNAR** — Heckman-based multiple imputation under MNAR
- **modelsummary, kableExtra** — publication-quality tables

## Data

Place the original replication kit files in `data/raw/`:

```
data/raw/
├── Interventions master_excel__Sep2023_JWAU__MS.xlsx
├── Unique Crisis & Unique Interventions.xlsx
├── Conversions_country_code.xlsx
└── mpd2020.xlsx
```

Run `scripts/00_etl.R` once. It produces `data/processed/analysis_data.csv`:
910 crisis-level rows with outcomes, covariates, selection indicator `R`,
and exclusion-restriction variables (`n_chron_clean`, `maddison_priority`).

## References

- Metrick, A. & Schmelzing, P. (2024), NBER w29281.
- Heckman, J. (1979), *Econometrica* 47(1).
- Toomet, O. & Henningsen, A. (2008), *J. Stat. Software* 27(7).
- Marra, G. & Radice, R. (2017), *GJRM* package.
- Galimard, J-E. et al. (2016, 2018), *Statist. Med.* & *BMC Med. Res. Meth.*
- Little, R.J.A. & Rubin, D.B. (2019), *Statistical Analysis with Missing Data*, 3e.
- van Hasselt, M. (2011), *J. Econometrics* 165(2).
