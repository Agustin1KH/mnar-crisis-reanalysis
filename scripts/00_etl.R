# scripts/00_etl.R
# Build the clean crisis-level analysis dataset from the original
# Metrick-Schmelzing replication kit files.
#
# Inputs  (place in data/raw/):
#   - Interventions master_excel__Sep2023_JWAU__MS.xlsx
#   - Unique Crisis & Unique Interventions.xlsx
#   - Conversions_country_code.xlsx
#
# Output:
#   - data/processed/analysis_data.csv (910 rows, ~28 cols)
#
# The output has:
#   * Paper's regressors: log_lagged_gdppc, polity, currency_regime, exp, debt
#   * Outcomes: gap_sum and the seven intervention dummies (guarantees_d etc.)
#   * Selection indicator r (= 1 iff all Table 3 regressors observed)
#   * Exclusion-restriction candidates: n_chron_clean, maddison_priority
#   * Era dummy and iso3

source(here::here("R", "utils.R"))
load_libs()

# ---- Paths ----
raw_dir <- here::here("data", "raw")
out_path <- here::here("data", "processed", "analysis_data.csv")

master_path  <- file.path(raw_dir, "Interventions master_excel__Sep2023_JWAU__MS.xlsx")
crisis_path  <- file.path(raw_dir, "Unique Crisis & Unique Interventions.xlsx")
conv_path    <- file.path(raw_dir, "Conversions_country_code.xlsx")

stopifnot(file.exists(master_path), file.exists(crisis_path), file.exists(conv_path))

# ---- Load ----
# Master list has a 4-row banner; actual header is row 5.
master <- readxl::read_excel(master_path, sheet = "Master list", skip = 4)

crisis <- readxl::read_excel(crisis_path, sheet = "Unique_Crisis_Codes")

conv <- readxl::read_excel(conv_path, sheet = "conversions_country_code_1") |>
  dplyr::rename(master_prefix = before, iso3 = after)

# ---- Build chronology coverage from master list ----
# The master file's "Crisis Code" uses 2-3 letter country prefixes (IT, UK, ESP...);
# the unique-crisis file uses ISO3 (ITA, GBR, ESP...). Bridge via `conv`.

chron_col <- "B/V/X, L/V, R/R or S/T?"

master_keyed <- master |>
  dplyr::mutate(
    master_prefix = stringr::str_extract(`Crisis Code`, "^[A-Z]{2,3}"),
    master_year   = as.numeric(stringr::str_extract(`Crisis Code`, "\\d{3,4}"))
  ) |>
  dplyr::left_join(conv, by = "master_prefix") |>
  dplyr::mutate(new_key = paste0(iso3, "-", master_year)) |>
  dplyr::filter(!is.na(iso3), !is.na(master_year))

# Parse chronology flags per row, then union within each crisis key
master_keyed$chrons <- lapply(master_keyed[[chron_col]], parse_chronologies)

chron_by_key <- master_keyed |>
  dplyr::filter(lengths(chrons) > 0) |>
  dplyr::group_by(new_key) |>
  dplyr::summarise(chrons = list(unique(unlist(chrons))), .groups = "drop")

chron_lookup <- setNames(chron_by_key$chrons, chron_by_key$new_key)

has_chron <- function(crisis_code, tag) {
  flags <- chron_lookup[[crisis_code]]
  if (is.null(flags)) 0L else as.integer(tag %in% flags)
}

crisis <- crisis |>
  dplyr::mutate(
    chron_BVX = vapply(`Crisis Code`, has_chron, integer(1), tag = "B/V/X"),
    chron_LV  = vapply(`Crisis Code`, has_chron, integer(1), tag = "L/V"),
    chron_RR  = vapply(`Crisis Code`, has_chron, integer(1), tag = "R/R"),
    chron_ST  = vapply(`Crisis Code`, has_chron, integer(1), tag = "S/T"),
    n_chronologies = chron_BVX + chron_LV + chron_RR + chron_ST,
    # Fallback: any crisis with Candidate = 0 is canonical by construction;
    # if the chronology tag wasn't populated in the master, impute 1.
    n_chron_clean = dplyr::if_else(
      n_chronologies == 0L & Candidate == 0L, 1L, n_chronologies
    ),
    # Phase-2 alternative shadow variables. Both strip Laeven-Valencia (LV)
    # from the count, so each ranges over {0, 1, 2, 3} rather than {0, 1, 2,
    # 3, 4}. The "raw" version makes no fallback imputation; the "cleaned"
    # version applies the same canonical-crisis fallback rule used by
    # n_chron_clean (i.e. if zero-count and Candidate == 0, impute 1).
    # Used by the SHADOW_VAR pipeline switch in scripts/{02,03,03a,03b,04,
    # 05} to test whether stripping LV resolves the in-sample exclusion-
    # restriction violation flagged in output/tables/03b_insample_exclusion_test.md.
    # See files (1)/cursor-phase2-task-spec.md.
    n_chron_3chron_raw = chron_BVX + chron_RR + chron_ST,
    n_chron_3chron = dplyr::if_else(
      n_chron_3chron_raw == 0L & Candidate == 0L, 1L, n_chron_3chron_raw
    )
  )

# ---- Maddison-priority dummy ----
crisis <- crisis |>
  dplyr::mutate(
    iso3 = stringr::str_extract(`Crisis Code`, "^[A-Z]{3}"),
    maddison_priority = as.integer(iso3 %in% maddison_priority_iso())
  )

# ---- Phase-2 Task-2 dichotomized FX regime ---------------------------------
# `Currency Regime` is the 15-level Reinhart-Rogoff IRR fine classification
# (1 = no separate legal tender, ..., 15 = freely floating). It is ~54%
# missing on the post-1800 sample; the original ETL keeps it as numeric and
# the paper treats it as such. For Phase-2 Task 2 (NARFCS binary MNAR
# sensitivity, see files (1)/cursor-phase2-task-spec.md), we add a binary
# dichotomization at threshold 7, matching the paper's own Figure 4
# treatment ("liberal vs fixed"; >=7 corresponds to managed floating and
# above). Missingness pattern is identical to `Currency Regime` itself
# (carries through one-to-one). Used by scripts/06_narfcs_binary_fx.R via
# mice::mnar.logreg.
crisis <- crisis |>
  dplyr::mutate(
    currency_regime_binary = dplyr::case_when(
      is.na(`Currency Regime`)   ~ NA_integer_,
      `Currency Regime` >= 7L    ~ 1L,
      TRUE                        ~ 0L
    )
  )

# ---- Era bucketing (matches Table 1 in the paper) ----
era_breaks <- c(-Inf, 1799, 1869, 1913, 1945, 1971, 1999, Inf)
era_labels <- c("pre1800", "1800_1869", "1870_1913", "1914_1945",
                "1946_1971", "1972_1999", "2000_2019")
crisis <- crisis |>
  dplyr::mutate(era = cut(year, breaks = era_breaks, labels = era_labels))

# ---- Selection indicator: all Table 3 regressors observed ----
table3_regs <- c("Log Lagged gdppc", "Polity", "Currency Regime", "Exp", "Debt")
crisis <- crisis |>
  dplyr::mutate(
    r = as.integer(rowSums(is.na(dplyr::across(dplyr::all_of(table3_regs)))) == 0)
  )

# ---- Diagnostic: confirm the MNAR pattern ----
cat("\n=== Selection diagnostic ===\n")
complete <- crisis |> dplyr::filter(r == 1)
dropped  <- crisis |> dplyr::filter(r == 0)
cat(sprintf("Complete-case sample: %d / %d (%.1f%%)\n",
            nrow(complete), nrow(crisis), 100 * nrow(complete) / nrow(crisis)))
cat(sprintf("  median year, kept:    %.0f\n", median(complete$year, na.rm = TRUE)))
cat(sprintf("  median year, dropped: %.0f\n", median(dropped$year, na.rm = TRUE)))

# ---- Select and rename, save ----
out <- crisis |>
  dplyr::select(
    crisis_code = `Crisis Code`,
    where = Where,
    iso3,
    year,
    era,
    candidate = Candidate,
    r,
    log_lagged_gdppc = `Log Lagged gdppc`,
    polity = Polity,
    currency_regime = `Currency Regime`,
    exp = Exp,
    debt = Debt,
    gap_sum = `Gap Sum`,
    guarantees_d = GUARANTEES_d,
    lending_d = LENDING_d,
    capital_injections_d = CAPITAL_INJECTIONS_d,
    restructuring_d = RESTRUCTURING_d,
    asset_management_d = ASSET_MANAGEMENT_d,
    rules_d = RULES_d,
    noi_d = NOI_d,
    other_d = OTHER_d,
    # Chron columns are lowercase in the canonical analysis_data.csv
    # (matches the convention used by scripts_py/* and the historical
    # downstream R output). Internally they're stored uppercase per the
    # paper's notation; we lowercase them on write.
    chron_bvx = chron_BVX, chron_lv = chron_LV,
    chron_rr  = chron_RR,  chron_st = chron_ST,
    n_chronologies, n_chron_clean,
    n_chron_3chron, n_chron_3chron_raw,
    maddison_priority,
    currency_regime_binary
  )

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

# ---- Phase-2 safety rail: backup + diff existing columns --------------------
# Before overwriting analysis_data.csv, snapshot the current version (if any)
# so we can verify the additive Phase-2 changes (n_chron_3chron and
# n_chron_3chron_raw) don't perturb any pre-existing column. The backup is
# written once and then left alone; if it already exists, we don't overwrite
# it (the snapshot is the pre-Phase-2 ground truth).
backup_path <- file.path(dirname(out_path), "analysis_data_pre_phase2.csv")
prev_existed <- file.exists(out_path)
if (prev_existed && !file.exists(backup_path)) {
  file.copy(out_path, backup_path, overwrite = FALSE)
  cat(sprintf("Backed up pre-Phase-2 CSV: %s\n", backup_path))
}

write.csv(out, out_path, row.names = FALSE)
cat(sprintf("\nSaved: %s   (%d rows, %d cols)\n", out_path, nrow(out), ncol(out)))

# ---- Diff check on existing columns ----------------------------------------
# If a pre-Phase-2 backup is present, re-read the freshly written CSV and
# the backup, then compare every column the backup has. Any difference is a
# bug (the only intended Phase-2 changes are the two additive shadow
# columns + their invariants). Halt loudly if so, except for known pre-
# existing drifts in the ETL code (documented in EXPECTED_DIFFS_KNOWN).
#
# EXPECTED_DIFFS_KNOWN: pre-existing drifts between the on-disk
# pre-Phase-2 CSV and the current ETL code that pre-date the SHADOW_VAR
# work and do not affect any analysis number. Each entry must be
# accompanied by a one-line reason. If you add an entry, also note it in
# memo/memo.tex if the drift is reader-visible.
EXPECTED_DIFFS_KNOWN <- c(
  where = "UTF-8 encoding diff for 'Cote d'Ivoire' (row 636); cosmetic, no analysis impact.",
  iso3  = "stringr::str_extract NA encoded as empty string in pre-Phase-2 CSV vs 'NA' literal now (rows 547, 818); maddison_priority membership test handles both identically.",
  era   = "Era label format changed from '1800-69' to '1800_1869' style (cosmetic; cut() bins are unchanged so all groupings are identical)."
)

if (file.exists(backup_path)) {
  prev <- read.csv(backup_path, stringsAsFactors = FALSE)
  curr <- read.csv(out_path,    stringsAsFactors = FALSE)
  shared_cols <- intersect(names(prev), names(curr))
  diffs <- character()
  for (col in shared_cols) {
    p <- prev[[col]]; c <- curr[[col]]
    # Compare with NA-aware equality.
    same <- length(p) == length(c) &&
            all((is.na(p) & is.na(c)) | (!is.na(p) & !is.na(c) & p == c))
    if (!isTRUE(same)) diffs <- c(diffs, col)
  }
  new_cols <- setdiff(names(curr), names(prev))
  dropped  <- setdiff(names(prev), names(curr))
  cat(sprintf("Diff vs pre-Phase-2 backup: %d shared cols checked, %d differ; new cols added: %s; dropped: %s.\n",
              length(shared_cols), length(diffs),
              if (length(new_cols)) paste(new_cols, collapse = ", ") else "(none)",
              if (length(dropped))  paste(dropped,  collapse = ", ") else "(none)"))

  expected_known <- intersect(diffs, names(EXPECTED_DIFFS_KNOWN))
  unexpected     <- setdiff(diffs, names(EXPECTED_DIFFS_KNOWN))

  if (length(expected_known) > 0L) {
    cat(sprintf("\nKnown pre-existing drifts (allow-listed; not Phase-2 caused):\n"))
    for (col in expected_known) {
      cat(sprintf("  - %-6s : %s\n", col, EXPECTED_DIFFS_KNOWN[[col]]))
    }
  }

  if (length(unexpected) > 0L) {
    stop(sprintf(
      "Phase-2 ETL safety rail tripped: %d existing column(s) changed unexpectedly: %s. ",
      length(unexpected), paste(unexpected, collapse = ", ")
    ), "These are not in EXPECTED_DIFFS_KNOWN and may be a Phase-2 bug. ",
       "Investigate before proceeding. Restore from ",
       "data/processed/analysis_data_pre_phase2.csv if you need the original.",
       call. = FALSE)
  }
  if (length(dropped) > 0L) {
    stop(sprintf(
      "Phase-2 ETL safety rail tripped: %d existing column(s) dropped: %s. ",
      length(dropped), paste(dropped, collapse = ", ")
    ), "The Phase-2 ETL change is meant to be additive only.",
    call. = FALSE)
  }

  required_new <- c("n_chron_3chron_raw", "n_chron_3chron",
                    "currency_regime_binary")
  missing_new <- setdiff(required_new, new_cols)
  if (length(missing_new) > 0L) {
    stop(sprintf(
      "Phase-2 ETL safety rail tripped: required new column(s) missing: %s.",
      paste(missing_new, collapse = ", ")
    ), call. = FALSE)
  }
}

# ---- Contractual invariants (additive: assertions only) ---------------------
# These guard the dataset shape that every downstream script relies on.
# They run on the freshly written `out` table, so any future change to the
# transformation logic above must keep these invariants intact or update them
# deliberately (and update scripts/01c_compare_to_paper.R accordingly).
checkmate::assert_data_frame(out, nrows = 910L)
checkmate::assert_true(sum(out$r) >= 260L && sum(out$r) <= 280L)
checkmate::assert_true(!anyDuplicated(out$crisis_code))
checkmate::assert_integerish(out$year, lower = 33L, upper = 2019L, any.missing = FALSE)
checkmate::assert_integerish(out$candidate, lower = 0L, upper = 1L, any.missing = FALSE)
checkmate::assert_integerish(out$n_chron_clean, lower = 0L, upper = 4L, any.missing = FALSE)
# Phase-2 additive invariants. Both shadows strip LV (max count = 3) and
# either apply or skip the canonical-crisis fallback rule used by
# n_chron_clean. Range {0, 1, 2, 3} for both; no missingness.
checkmate::assert_integerish(out$n_chron_3chron_raw, lower = 0L, upper = 3L,
                             any.missing = FALSE)
checkmate::assert_integerish(out$n_chron_3chron, lower = 0L, upper = 3L,
                             any.missing = FALSE)
# currency_regime_binary inherits missingness from `Currency Regime` itself
# (~54% missing post-1800); allow NA. Range {0, 1, NA}.
checkmate::assert_integerish(out$currency_regime_binary, lower = 0L,
                             upper = 1L, any.missing = TRUE)
cat(sprintf("Invariants OK: nrow=%d, sum(r)=%d, year=[%d,%d]\n",
            nrow(out), sum(out$r), min(out$year), max(out$year)))
