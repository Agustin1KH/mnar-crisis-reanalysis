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
