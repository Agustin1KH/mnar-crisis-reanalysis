# R/utils.R
# Shared helper functions for the MNAR re-analysis.

#' Load common libraries with a consistent conflict policy.
#' Call once at the top of any analysis script.
load_libs <- function() {
  suppressPackageStartupMessages({
    library(here)
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(modelsummary)
  })
  options(
    dplyr.summarise.inform = FALSE,
    readr.show_col_types = FALSE
  )
}

#' Parse the "B/V/X, L/V, R/R or S/T?" column of the master spreadsheet
#' into a character vector of chronology flags. Handles the "B/V//X" typo.
parse_chronologies <- function(cell) {
  if (is.na(cell)) return(character(0))
  cell <- stringr::str_replace_all(cell, "B/V//X", "B/V/X")
  parts <- stringr::str_trim(stringr::str_split(cell, ";")[[1]])
  parts[parts %in% c("B/V/X", "L/V", "R/R", "S/T")]
}

#' The 8 "leading DM" countries per Schmelzing (2020) BoE SWP 845.
#' Used to flag whether a crisis occurred in a Maddison-priority economy.
maddison_priority_iso <- function() {
  c("ITA", "GBR", "NLD", "FRA", "DEU", "ESP", "USA", "JPN")
}

#' Short diagnostic: given a data frame and a set of regressor columns,
#' print the complete-case sample size and the median year of dropped vs kept.
report_selection <- function(df, regs, year_col = "year") {
  complete <- tidyr::drop_na(df, dplyr::all_of(regs))
  dropped <- dplyr::anti_join(df, complete, by = "crisis_code")
  cat(sprintf(
    "Complete-case: %d / %d (%.1f%%)\n  kept median year: %.0f\n  dropped median year: %.0f\n",
    nrow(complete), nrow(df), 100 * nrow(complete) / nrow(df),
    median(complete[[year_col]], na.rm = TRUE),
    median(dropped[[year_col]], na.rm = TRUE)
  ))
  invisible(list(complete = complete, dropped = dropped))
}

#' Phase-2 SHADOW_VAR resolver.
#'
#' Returns a list with `var` (the chosen shadow-variable column name),
#' `is_default` (whether the SHADOW_VAR env var was unset, in which case we
#' fall back to "n_chron_clean"), and `subdir` (the relative path under
#' `output/tables/` to namespace Phase-2 outputs into; empty string when
#' is_default = TRUE so legacy outputs don't move). Also validates the
#' chosen column is one of the three expected shadows.
#'
#' Usage in each consumer script:
#'   sv  <- resolve_shadow_var()
#'   subdir <- sv$subdir   # "" or e.g. "phase2_shadow_n_chron_3chron_raw"
#'   sel_eq <- as.formula(paste0("r ~ candidate + year + maddison_priority + ", sv$var))
resolve_shadow_var <- function(allowed = c("n_chron_clean",
                                            "n_chron_3chron",
                                            "n_chron_3chron_raw")) {
  raw <- Sys.getenv("SHADOW_VAR", unset = "")
  is_default <- !nzchar(raw)
  var <- if (is_default) "n_chron_clean" else raw
  if (!var %in% allowed) {
    stop(sprintf("SHADOW_VAR='%s' is not in {%s}; refusing to proceed.",
                 var, paste(allowed, collapse = ", ")), call. = FALSE)
  }
  subdir <- if (is_default) "" else paste0("phase2_shadow_", var)
  list(var = var, is_default = is_default, subdir = subdir)
}

#' Phase-2 output-path helper. Concatenates `output/tables/<subdir>/<file>`
#' (or `output/tables/<file>` when subdir is ""), creating the directory
#' if needed. Pass either a single filename or a relative path.
phase2_output_path <- function(file, subdir = "", kind = c("tables", "figures"),
                                ensure_dir = TRUE) {
  kind <- match.arg(kind)
  parts <- c("output", kind)
  if (nzchar(subdir)) parts <- c(parts, subdir)
  dir <- do.call(here::here, as.list(parts))
  if (ensure_dir) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, file)
}

#' Cache-path counterpart for phase2_output_path: Phase-2 RDS / cached
#' fits go under `data/processed/<subdir>/`. Same defaulting behaviour
#' (subdir "" -> the existing flat layout).
phase2_cache_path <- function(file, subdir = "", ensure_dir = TRUE) {
  parts <- c("data", "processed")
  if (nzchar(subdir)) parts <- c(parts, subdir)
  dir <- do.call(here::here, as.list(parts))
  if (ensure_dir) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, file)
}
