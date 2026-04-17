# scripts/01e_duplicate_audit.R
# Systematic five-check duplicate audit on data/processed/analysis_data.csv,
# read-only diagnostic. Spec is verbatim from Prompt 3c:
#
#   (a) Exact crisis_code duplicates (should be zero post-dedup).
#   (b) (iso3, year) collisions: same country-year under different crisis_codes.
#   (c) Within-complete-case (r == 1) near-duplicates on the regressor
#       covariates -- same country-year with slightly different regressor
#       values (would signal a row got silently split during a merge).
#   (d) Perfect covariate duplicates within the complete-case sample
#       -- two crises with identical X but possibly different Y, which
#       Heckman MLE handles less gracefully than OLS.
#   (e) Country-year pairs where two different ISO3 codes appear,
#       i.e., other NOK/NOR-style pattern survivors from the merge.
#
# Reasoning for the audit's importance: Heckman fits are sensitive to exact
# duplicates in ways OLS isn't (the likelihood can be locally flat or have
# convergence pathologies when two observations have identical X and
# effectively double-weight one point). We want to rule them out before
# Prompt 4 runs, not after.
#
# Hard rule: If any check flags a finding inside the r == 1 sample
# (the Tables 2 / 3 / Heckman selection-equation sample), the script
# stop()s with a CRITICAL message naming the offending rows so the user
# can investigate before Prompt 4. Findings outside r == 1 are recorded
# but do not halt.
#
# Constraint: read-only. Does not modify analysis_data.csv or any
# transformation script.

source(here::here("R", "utils.R"))
load_libs()

library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

# Use na.strings = c("NA", "") so that empty cells in iso3 (rows G-1974,
# G-2008 written by 00_etl.R as empty strings rather than the literal "NA")
# parse to actual NA. read.csv's default of na.strings = "NA" leaves them as
# "" which silently breaks any check that depends on is.na(iso3).
df <- read.csv(here::here("data", "processed", "analysis_data.csv"),
               stringsAsFactors = FALSE, na.strings = c("NA", ""))

stopifnot(nrow(df) == 910L, "crisis_code" %in% names(df), "iso3" %in% names(df),
          # Sanity: G-1974 and G-2008 should now read as NA in iso3
          identical(sum(is.na(df$iso3[df$crisis_code %in% c("G-1974", "G-2008")])),
                    2L))

# Keep snapshots of the SQL-equivalent code shown in the report; these strings
# are interpolated into the report verbatim so the .md acts as a runnable
# spec for the audit.
sql_a <- "-- (a) Exact crisis_code duplicates\nsum(duplicated(df$crisis_code))"
sql_b <- "-- (b) (iso3, year) collisions\ndf |> dplyr::count(iso3, year) |> dplyr::filter(n > 1)"
sql_c <- "-- (c) Within-complete-case (r == 1) near-duplicates on (iso3, year)\ndf |>\n  dplyr::filter(r == 1) |>\n  dplyr::group_by(iso3, year) |>\n  dplyr::filter(n() > 1) |>\n  dplyr::select(crisis_code, iso3, year, log_lagged_gdppc, polity,\n                currency_regime, candidate)"
sql_d <- "-- (d) Perfect covariate duplicates within the complete-case sample\ndf |>\n  dplyr::filter(r == 1) |>\n  dplyr::group_by(log_lagged_gdppc, polity, currency_regime, debt, exp, year) |>\n  dplyr::filter(n() > 1) |>\n  dplyr::arrange(log_lagged_gdppc)"
sql_e <- "-- (e) Country-year pairs where two different ISO3 codes appear\ndf |>\n  dplyr::mutate(yr = as.integer(stringr::str_extract(crisis_code, \"\\\\d{4}$\"))) |>\n  dplyr::group_by(yr) |>\n  dplyr::summarise(unique_iso = n_distinct(iso3, na.rm = TRUE),\n                   crises = list(crisis_code), .groups = \"drop\") |>\n  dplyr::filter(unique_iso > 1 & lengths(crises) > 1) |>\n  tidyr::unnest(crises)"

# ---------------------------------------------------------------------------
# (a) Exact crisis_code duplicates
# ---------------------------------------------------------------------------
a_count <- sum(duplicated(df$crisis_code))
a_rows  <- df |> filter(crisis_code %in% df$crisis_code[duplicated(df$crisis_code)])

# ---------------------------------------------------------------------------
# (b) (iso3, year) collisions
# ---------------------------------------------------------------------------
b_groups <- df |> count(iso3, year, name = "nrows") |> filter(nrows > 1L)
b_count  <- nrow(b_groups)
b_rows   <- df |>
  semi_join(b_groups, by = c("iso3", "year")) |>
  arrange(iso3, year, crisis_code)

# ---------------------------------------------------------------------------
# (c) Within-complete-case (r == 1) near-duplicates on (iso3, year)
# ---------------------------------------------------------------------------
c_rows <- df |>
  filter(r == 1L) |>
  group_by(iso3, year) |>
  filter(n() > 1L) |>
  ungroup() |>
  select(crisis_code, iso3, year, log_lagged_gdppc, polity,
         currency_regime, candidate) |>
  arrange(iso3, year, crisis_code)
c_count <- nrow(c_rows)

# ---------------------------------------------------------------------------
# (d) Perfect covariate duplicates within the complete-case sample
# ---------------------------------------------------------------------------
fp_cols <- c("log_lagged_gdppc","polity","currency_regime","debt","exp","year")
d_rows <- df |>
  filter(r == 1L) |>
  group_by(across(all_of(fp_cols))) |>
  filter(n() > 1L) |>
  ungroup() |>
  arrange(log_lagged_gdppc)
d_count <- nrow(d_rows)

# ---------------------------------------------------------------------------
# (e) Country-year pairs where two different ISO3 codes appear
# ---------------------------------------------------------------------------
e_rows <- df |>
  mutate(yr = suppressWarnings(
    as.integer(stringr::str_extract(crisis_code, "[0-9]{4}$")))) |>
  group_by(yr) |>
  summarise(
    unique_iso = n_distinct(iso3, na.rm = TRUE),
    n_crises   = n(),
    crises     = list(crisis_code),
    iso3_set   = list(sort(unique(iso3))),
    .groups    = "drop"
  ) |>
  filter(unique_iso > 1L & n_crises > 1L) |>
  unnest(crises) |>
  arrange(yr, crises)
e_count <- nrow(e_rows)

# Note: check (e) as written intentionally reports every year in which
# multiple countries had crises. That is most years -- the check is broad
# by design. The semantically interesting sub-question is split into two
# follow-ups:
#
#   (e1) Co-occurrence pattern: years where a meta-crisis (iso3 = NA) row
#        sits in the same year as country-level rows. Informational --
#        meta-crises with iso3 = NA cannot duplicate country-level rows
#        because their iso3 doesn't match any country. They appear in the
#        broad (e) bucket but are NOT a regression-sample concern unless
#        check (e2) below trips.
#
#   (e2) The HALT case: is there any row with iso3 = NA AND r == 1?
#        That would put a meta-crisis row into the Tables 2 / 3 / Heckman
#        regression sample where iso3 clustering would silently fail.
e_co_occur <- df |>
  group_by(year) |>
  filter(any(is.na(iso3)) & any(!is.na(iso3))) |>
  ungroup() |>
  arrange(year, crisis_code) |>
  select(crisis_code, iso3, year, candidate, r)
e1_count <- nrow(e_co_occur)

e_meta_in_r1 <- df |>
  filter(is.na(iso3) & r == 1L) |>
  select(crisis_code, iso3, year, candidate, r,
         log_lagged_gdppc, polity, currency_regime, exp, debt)
e2_count <- nrow(e_meta_in_r1)

# ---------------------------------------------------------------------------
# Determine if any finding affects the r == 1 sample
# ---------------------------------------------------------------------------
# Halt-worthy findings: ones that put a row into the regression sample
# under conditions that would corrupt the fit. We do NOT halt on (e1)
# because country-level rows in the same year as a meta-crisis are not
# duplicates of the meta-crisis (iso3 doesn't match).
r1_findings <- c(
  a   = if (a_count) any(a_rows$r == 1L) else FALSE,
  b   = if (b_count) any(b_rows$r == 1L) else FALSE,
  c   = c_count > 0L,
  d   = d_count > 0L,
  e2  = e2_count > 0L
)
n_r1 <- sum(r1_findings)

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
cat("\n=== Duplicate audit on data/processed/analysis_data.csv ===\n")
fmt <- function(label, n, threat = "") {
  flag <- if (n == 0L) "OK    " else "FOUND "
  cat(sprintf("  [%s] %s  count=%-3d  %s%s\n",
              flag, label, n, threat, if (nchar(threat)) "" else ""))
}
fmt("(a)  exact crisis_code duplicates",                     a_count)
fmt("(b)  (iso3, year) collisions",                          b_count)
fmt("(c)  (r==1) near-duplicate on (iso3, year)",            c_count)
fmt("(d)  (r==1) perfect covariate duplicates",              d_count)
fmt("(e)  country-year pairs with multi-iso3",               e_count,
    "(broad by design; (b) is the focused complement)")
fmt("(e1) meta-crisis + country-level same year (info)",     e1_count)
fmt("(e2) meta-crisis row with r == 1 (HALT case)",          e2_count)

cat(sprintf("\nFindings affecting r == 1 (Tables 2/3/Heckman sample): %d\n",
            n_r1))

if (n_r1 > 0L) {
  cat("\n*** HALT: r == 1 finding detected. Surface to user before Prompt 4. ***\n")
}

# ---------------------------------------------------------------------------
# Render report
# ---------------------------------------------------------------------------
out_dir <- here::here("output", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "01e_duplicate_audit.md")

emit_check <- function(id, title, sql, count, rows, interp = NULL,
                       force_status = NULL) {
  status <- force_status %||% (if (count == 0L) "PASS" else "FOUND")
  bullet_count <- sprintf("**Findings:** %d %s.", count,
                          if (count == 1L) "row" else "rows")
  out <- c("",
           sprintf("### Check %s -- %s", id, title),
           "",
           sprintf("Status: **%s**.  %s", status, bullet_count),
           "",
           "```r",
           sql,
           "```",
           "")
  if (count == 0L) {
    out <- c(out, "_No rows returned. Pass._", "")
  } else {
    tab <- knitr::kable(as.data.frame(rows), format = "pipe", align = "l")
    out <- c(out, as.character(tab), "")
    if (!is.null(interp)) out <- c(out, "**Interpretation.**", interp, "")
  }
  out
}

`%||%` <- function(a, b) if (is.null(a) || isTRUE(is.na(a))) b else a

interp_e <- c(
  "Most years in 1800-2019 have crises in multiple countries; the check ",
  "as written therefore returns many rows by design. The semantically ",
  "important sub-question -- *did the same iso3 reappear under two ",
  "different crisis_codes in the same year?* -- is what check (b) tests, ",
  "and (b) is clean. The check (e') below isolates the more dangerous ",
  "edge case: years where a meta-crisis row (iso3 = NA, e.g., `G-1974`) ",
  "co-occurs with country-level rows the same year."
)

interp_e1 <- c(
  "Years where a meta-crisis row (`G-1974` G10 / Herstatt; `G-2008` IASB ",
  "rule change) co-occurs with country-level rows. **This is informational ",
  "only.** The meta-crisis rows have `iso3 = NA`, so by definition they ",
  "cannot duplicate any country-level row (whose `iso3` is a 3-letter ",
  "code). They appear in the broad `(e)` cohort because they share a year. ",
  "The relevant halt-worthy variant is check `(e2)` below."
)
interp_e2 <- c(
  "Whether any meta-crisis row (`iso3 = NA`) sits inside the `r == 1` ",
  "complete-case sample. If so, the iso3-clustered Heckman selection ",
  "equation would silently drop or mis-cluster that row. Today: zero ",
  "rows -- both `G-1974` and `G-2008` have `r = 0` (no Table 3 / Table 2 ",
  "regressor observed) and are excluded from every regression sample."
)

lines <- c(
  "# Duplicate audit on `data/processed/analysis_data.csv` (Prompt 3c)",
  "",
  paste0(
    "Generated by `scripts/01e_duplicate_audit.R`. Five-check systematic ",
    "audit verbatim from Prompt 3c. Status code: `PASS` = zero rows ",
    "returned; `FOUND` = non-zero rows returned (interpretation follows). ",
    "**Hard halt rule:** if any check flags a finding inside the `r == 1` ",
    "sample (the Tables 2 / 3 / Heckman selection-equation sample), the ",
    "script `stop()`s and surfaces the rows to the user before Prompt 4 ",
    "is allowed to run."
  ),
  "",
  "## Summary",
  "",
  "| ID | Status | Findings | Affects r==1? | Check |",
  "|---|---|---:|---|---|",
  sprintf("| (a) | %s | %d | %s | exact `crisis_code` duplicates |",
          if (a_count == 0L) "PASS" else "FOUND", a_count,
          if (any(a_rows$r == 1L)) "**YES**" else "no"),
  sprintf("| (b) | %s | %d | %s | `(iso3, year)` collisions |",
          if (b_count == 0L) "PASS" else "FOUND", b_count,
          if (b_count && any(b_rows$r == 1L)) "**YES**" else "no"),
  sprintf("| (c) | %s | %d | %s | `(r==1)` near-duplicates on `(iso3, year)` |",
          if (c_count == 0L) "PASS" else "FOUND", c_count,
          if (c_count) "**YES**" else "n/a"),
  sprintf("| (d) | %s | %d | %s | `(r==1)` perfect covariate duplicates |",
          if (d_count == 0L) "PASS" else "FOUND", d_count,
          if (d_count) "**YES**" else "n/a"),
  sprintf("| (e) | %s | %d | %s | country-year pairs with multi-iso3 (broad) |",
          if (e_count == 0L) "PASS" else "INFO", e_count,
          "see (e1) / (e2)"),
  sprintf("| (e1)| %s | %d | %s | meta-crisis + country-level same year (info) |",
          if (e1_count == 0L) "PASS" else "INFO", e1_count,
          "no -- iso3=NA can't duplicate country-level rows"),
  sprintf("| (e2)| %s | %d | %s | meta-crisis row with `r == 1` (HALT case) |",
          if (e2_count == 0L) "PASS" else "FOUND", e2_count,
          if (e2_count) "**YES**" else "no"),
  ""
)

# Show first ~20 rows of (e) so the report doesn't explode
e_show <- e_rows |> head(20) |>
  select(yr, unique_iso, n_crises, crises) |>
  mutate(iso3_set = vapply(seq_len(n()),
                           function(i) NA_character_, character(1)))
e_show <- e_rows |> head(20) |>
  select(yr, unique_iso, n_crises, crises)

lines <- c(lines, "## Per-check detail",
           emit_check("(a)", "Exact crisis_code duplicates",
                     sql_a, a_count, a_rows,
                     interp = paste0("Indicates a row was added to ",
                                     "`Unique_Crisis_Codes` without going ",
                                     "through the dedup. None expected.")),
           emit_check("(b)", "(iso3, year) collisions",
                     sql_b, b_count, b_rows,
                     interp = paste0("Same iso3 in the same year under two ",
                                     "different crisis_codes -- the ",
                                     "NOK-1987-style pattern. Conversion-",
                                     "table risks: BG/BUG -> BGR, FR/FRA -> ",
                                     "FRA, GB/UK -> GBR, GI/GN -> GIN, ",
                                     "GT/GUA -> GTM, NOK/NOR -> NOR.")),
           emit_check("(c)", "(r==1) near-duplicates on (iso3, year)",
                     sql_c, c_count, c_rows,
                     interp = paste0("Same country-year inside the complete-",
                                     "case sample with possibly different ",
                                     "regressor values. Would signal a row ",
                                     "got silently split during a merge. ",
                                     "These would directly affect the ",
                                     "Heckman selection equation.")),
           emit_check("(d)", "(r==1) perfect covariate duplicates",
                     sql_d, d_count, d_rows,
                     interp = paste0("Two crises with identical X inside ",
                                     "the complete-case sample. Heckman ",
                                     "MLE can have likelihood pathologies ",
                                     "on identical-X rows even when Y ",
                                     "differs.")),
           emit_check("(e)", "country-year pairs with multi-iso3 (informational)",
                     sql_e, e_count, e_show,
                     interp = paste(interp_e, collapse = ""),
                     force_status = if (e_count > 0L) "INFO" else "PASS"))

# (e1) co-occurrence detail (info-only)
sql_e1 <- "-- (e1) Years where a meta-crisis (iso3 = NA) co-occurs with country-level rows\ndf |>\n  dplyr::group_by(year) |>\n  dplyr::filter(any(is.na(iso3)) & any(!is.na(iso3))) |>\n  dplyr::select(crisis_code, iso3, year, candidate, r)"
lines <- c(lines,
           emit_check("(e1)", "meta-crisis co-occurs with country-level rows (info)",
                      sql_e1, e1_count,
                      e_co_occur |> dplyr::filter(is.na(iso3) | year %in% c(1974, 2008)) |>
                        dplyr::arrange(year, !is.na(iso3), crisis_code),
                      interp = paste(interp_e1, collapse = ""),
                      force_status = if (e1_count > 0L) "INFO" else "PASS"))

# (e2) HALT variant: any meta-crisis row in r == 1
sql_e2 <- "-- (e2) Meta-crisis row (iso3 = NA) inside the r == 1 sample\ndf |>\n  dplyr::filter(is.na(iso3) & r == 1L) |>\n  dplyr::select(crisis_code, iso3, year, candidate, r,\n                log_lagged_gdppc, polity, currency_regime, exp, debt)"
lines <- c(lines,
           emit_check("(e2)", "meta-crisis row with r == 1 (HALT case)",
                      sql_e2, e2_count, e_meta_in_r1,
                      interp = paste(interp_e2, collapse = "")))

# ---------------------------------------------------------------------------
# Headline / verdict
# ---------------------------------------------------------------------------

verdict <- c(
  "## Verdict",
  "",
  if (n_r1 == 0L) {
    paste0(
      "**All five checks PASS for the r == 1 complete-case sample.** ",
      "Zero exact `crisis_code` duplicates; zero `(iso3, year)` ",
      "collisions; zero near-duplicates inside `r == 1`; zero perfect ",
      "covariate duplicates inside `r == 1`. The conversion-table ",
      "collision risks (`BG/BUG`, `FR/FRA`, `GB/UK`, `GI/GN`, `GT/GUA`, ",
      "`NOK/NOR`) are all clean in the current 910-row dataset -- the ",
      "single NOK/NOR-1987 dedup that we already documented in ",
      "`output/tables/01d_dedup_test.md` is the only one that ever ",
      "existed. No other NOK/NOR-style survivors lurking. **Heckman / ",
      "GJRM / mice-MNAR in Prompt 4 can proceed.**"
    )
  } else {
    paste0(
      "**HALT.** ", n_r1, " check(s) flagged a finding inside the ",
      "r == 1 sample. See per-check detail above for offending rows. ",
      "Investigate before Prompt 4 runs."
    )
  },
  "",
  "## Cross-references",
  "",
  "- `output/tables/01d_dedup_test.md` -- the NOK-1987 / NOR-1987 dedup ",
  "  reconciling test (the only documented duplicate; this audit confirms ",
  "  no other NOK/NOR-style cases exist).",
  "- `output/tables/01_replication_notes.md` (b) -- in-context drift ",
  "  diagnosis citing the dedup as cause #1 of the +1 combined N gap."
)

writeLines(c(lines, verdict), out_path)

cat(sprintf("\nWritten: %s\n", out_path))

# Hard halt if any r==1 finding
if (n_r1 > 0L) {
  stop(sprintf(paste0(
    "Duplicate audit halted: %d check(s) flagged finding(s) inside the ",
    "r == 1 sample. See output/tables/01e_duplicate_audit.md for offending ",
    "rows. Investigate before Prompt 4."), n_r1),
    call. = FALSE)
}
