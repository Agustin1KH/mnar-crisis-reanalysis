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
    )
  )

# ---- Maddison-priority dummy ----
crisis <- crisis |>
  dplyr::mutate(
    iso3 = stringr::str_extract(`Crisis Code`, "^[A-Z]{3}"),
    maddison_priority = as.integer(iso3 %in% maddison_priority_iso())
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
    chron_BVX, chron_LV, chron_RR, chron_ST,
    n_chronologies, n_chron_clean,
    maddison_priority
  )

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_path, row.names = FALSE)
cat(sprintf("\nSaved: %s   (%d rows, %d cols)\n", out_path, nrow(out), ncol(out)))

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
cat(sprintf("Invariants OK: nrow=%d, sum(r)=%d, year=[%d,%d]\n",
            nrow(out), sum(out$r), min(out$year), max(out$year)))
