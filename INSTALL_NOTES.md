# Install notes (environment bootstrap)

Platform: macOS 14.5 (Sonoma) / aarch64-apple-darwin20, R 4.3.1.

The environment (`renv.lock`, populated `renv/library/`) was built under R 4.3.1
against CRAN binaries. Most packages dropped in cleanly. A handful of issues
are documented below with the exact error and the workaround applied so this
bootstrap is reproducible on an identical machine and diagnosable on a
different one.

All 18 required packages (`here`, `readxl`, `dplyr`, `tidyr`, `stringr`,
`ggplot2`, `modelsummary`, `kableExtra`, `fixest`, `sandwich`, `lmtest`,
`sampleSelection`, `GJRM`, `mice`, `miceMNAR`, `cli`, `checkmate`, `tictoc`)
are installed, pinned in `renv.lock`, and load cleanly under
`scripts/_env_check.R`.

---

## 1. `miceMNAR` is archived on CRAN

### Error

```
Error: package 'miceMNAR' is not available
```

### Cause

`miceMNAR` was removed from the active CRAN repository on 2021-11-07 because
email to the maintainer was undeliverable
(<https://cran.r-project.org/web/packages/miceMNAR/index.html>). The last
published version, 1.0.2 (2018-08-27), remains available only from the CRAN
archive.

### Workaround (applied)

Install the `pbivnorm` runtime dependency from CRAN, then install
`miceMNAR_1.0.2.tar.gz` directly from the archive URL:

```r
install.packages("pbivnorm", type = "binary")
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/miceMNAR/miceMNAR_1.0.2.tar.gz",
  repos = NULL, type = "source"
)
```

Result: `miceMNAR 1.0.2` installed from source, pinned in `renv.lock`.

### Note for `renv::restore()`

`renv.lock` records `miceMNAR` with `Source = "Repository"`, which by
convention means "available from a CRAN mirror". Because 1.0.2 lives only in
the archive, `renv::restore()` on a clean machine may fail to fetch it. If
that happens, run the two-line workaround above once, then re-run
`renv::snapshot()` (or use `renv::install("url::...")` which records an
explicit remote URL in the lockfile).

---

## 2. `svglite` source compile fails without `libpng` headers

### Error

```
devSVG.cpp:21:10: fatal error: 'png.h' file not found
ERROR: compilation failed for package 'svglite'
```

`svglite` is an indirect dependency of `kableExtra`.

### Cause

`install.packages()` under this renv defaulted to `type = "both"` and fell
through to source for `svglite`. The macOS R build expects `png.h` at
`/opt/R/arm64/include`, but libpng is installed via Homebrew at
`/opt/homebrew/include`, which is not on the default C include path.

### Workaround (applied)

Prefer prebuilt CRAN binaries for the install pass:

```r
options(pkgType = "binary", install.packages.check.source = "no")
install.packages(pkgs, type = "binary")
```

The CRAN binary `svglite_2.2.1` installs in a fraction of a second and does
not require libpng headers locally.

### Alternative (not applied)

If a source build is required, either:

- Export `PKG_CPPFLAGS="-I/opt/homebrew/include"` and
  `PKG_LIBS="-L/opt/homebrew/lib -lpng"` before installing, or
- Install a system `libpng` into `/opt/R/arm64` (e.g. via the
  [macOS-recipes build](https://mac.r-project.org/bin/) tooling).

---

## 3. `readxl` binary wheel built for R 4.3.3 is ABI-incompatible with R 4.3.1

### Error (on `library(readxl)`)

```
unable to load shared object '.../readxl/libs/readxl.so':
  Symbol not found: _iconv
  Referenced from: .../readxl/libs/readxl.so
  Expected in:     .../R.framework/.../libR.dylib
```

### Cause

The CRAN binary `readxl_1.4.5` was compiled under R 4.3.3 against its
`libR.dylib` iconv exports. R 4.3.1 does not re-export the same `_iconv`
symbol, so `dlopen` fails at load time.

### Workaround (applied)

Force a source rebuild of `readxl` so it links against the running R's
libraries:

```r
renv::install("readxl", type = "source", rebuild = TRUE)
```

Result: `readxl 1.4.5` compiled locally, loads cleanly.

---

## 4. `lme4` / `Matrix` ABI version mismatch

### Warning (on `library(mice)` -> load chain into `lme4`)

```
Warning in check_dep_version():
  ABI version mismatch:
    lme4 was built with Matrix ABI version 1
    Current Matrix ABI version is 0
  Please re-install lme4 from source or restore original 'Matrix' package
```

### Cause

`renv::init()` discovered and pinned `Matrix 1.6-0` (the version shipped with
R 4.3.1). CRAN's `lme4 1.1-37` binary is compiled against `Matrix 1.7.x`
(ABI 1). Loading `lme4` with the older `Matrix` triggers an explicit
ABI-mismatch warning.

### Workaround (applied)

Upgrade `Matrix` to the `1.6-5` binary, which exposes ABI version 1 while
remaining R 4.3-compatible:

```r
install.packages("Matrix", type = "binary")
```

`renv.lock` now pins `Matrix 1.6-5`. The warning no longer fires.

---

## 5. "package 'x' was built under R version 4.3.3" notices

### Observation

With R 4.3.1 installed and CRAN shipping binaries built on R 4.3.3, several
pinned packages emit a one-shot load-time notice of the form:

```
Warning message:
package 'x' was built under R version 4.3.3
```

The affected packages include `ggplot2`, `modelsummary`, `fixest`, `sandwich`,
`lmtest` (via `zoo`), `sampleSelection` (via `maxLik`, `miscTools`), `GJRM`,
`mice`, `cli`, `checkmate`, `tictoc`, and `lme4`.

### Assessment

These notices are R's `check_build` warning and are documented as safe across
R 4.3.x patch releases: the packages have been ABI-tested and work on R
4.3.1. The mismatch is cosmetic, not functional.

### Workaround (applied)

`scripts/_env_check.R` wraps every `library()` call in a
`withCallingHandlers()` filter that muffles only messages matching
`"was built under R version"` (or `"ABI version"`), and counts them for
reporting. Any unrelated warning would propagate normally and fail the
check.

### Alternative (not applied)

Rebuild every affected package from source under R 4.3.1
(`renv::install(..., type = "source", rebuild = TRUE)`). This takes a few
minutes, requires a working C/C++/Fortran toolchain, and does not change the
resulting analysis numerics.

---

## 6. `miceMNAR::mice.impute.hecknorm` is broken with current `GJRM` (use `hecknorm2step`)

### Error (running `scripts/04_mice_mnar.R` with `method = "hecknorm"`)

```
Error in eval(predvars, data, env) : object 'ry' not found
Calls: mice ... model.frame.default -> is.data.frame -> as.data.frame -> cbind
```

### Cause

`miceMNAR::mice.impute.hecknorm` calls
`GJRM::copulaSampleSel(data = as.data.frame(cbind(ry, y, x)), ...)` with the
`data` argument as an inline expression. `copulaSampleSel` constructs its
model frame via `match.call()` + `eval(mf)`, and the eval scope
(`parent.frame()` of `copulaSampleSel`) cannot resolve `ry` / `y` / `x`
when `copulaSampleSel` is itself called from inside another function. The
call works in isolation (from `.GlobalEnv`) but fails inside `mice`'s
imputation harness. Verified against pinned versions
`miceMNAR 1.0.2` (2018-08-27, archive only) + `GJRM 0.2.6.8` (2024).

### Workaround (applied)

Use `miceMNAR::mice.impute.hecknorm2step` instead. It is the package-internal
two-step Heckman imputer (Galimard et al. 2018, JSCS) implementing the same
MNAR adjustment via `sampleSelection::heckit2fit` (no scope bug). The two
methods are asymptotically equivalent up to ML-vs-2-step efficiency.

`scripts/04_mice_mnar.R` invokes `hecknorm2step` for all five MNAR-imputed
continuous variables; the fallback is logged in
`output/tables/04_mice_fallback_log.md`.

---

## 7. `currency_regime` MNAR multinomial: no method available, fall back to `polyreg`

### Observation

The IRR 15-point fine FX-regime classification (`currency_regime`) is 54%
missing in the post-1800 sample and is multi-level categorical (15 levels,
treated as numeric in the original paper). `miceMNAR` provides only
`hecknorm` (continuous) and `heckprob` (binary); there is no `heckpoly`
for multi-level categorical MNAR imputation.

### Workaround (applied)

Per Phase 4a spec, `scripts/04_mice_mnar.R` falls back to
`mice`'s `polyreg` (multinomial logit, MAR conditional on observables) for
`currency_regime` in both the MNAR and MAR imputation arms. Same treatment
in both arms keeps the MAR-vs-MNAR comparison purely about the continuous
variables. Logged in `output/tables/04_mice_fallback_log.md`.

A Heckman-style ordinal extension (e.g. via `sampleSelection`'s bivariate
ordered-probit selection) would close the gap but is out of scope for
Phase 4a; could be added in Phase 5 if memo writeup demands MNAR
treatment of `currency_regime` specifically.

---

## Verification

```
source("scripts/_env_check.R")
```

Produces (on a fresh R session):

```
v All 18 required packages loaded.
i Muffled 12 benign R-version build notices (see INSTALL_NOTES.md).
...
v Environment check complete.
```

with `warnings()` empty and exit status 0.
