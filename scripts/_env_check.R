# scripts/_env_check.R
# Environment diagnostic for the MNAR re-analysis.
#
# Sources R/utils.R, attaches every required package, and prints sessionInfo()
# plus a tidy tibble of pinned versions. Exits cleanly (zero errors, zero
# warnings) when the renv library is healthy.
#
# Any "package 'x' was built under R version y" messages come from using binary
# CRAN builds under a slightly older R release; they are benign and are filtered
# here. See INSTALL_NOTES.md for the platform context.

required_pkgs <- c(
  "here", "readxl", "dplyr", "tidyr", "stringr", "ggplot2",
  "modelsummary", "kableExtra", "fixest", "sandwich", "lmtest",
  "sampleSelection", "GJRM", "mice", "miceMNAR",
  "cli", "checkmate", "tictoc"
)

is_benign_build_warning <- function(msg) {
  grepl("was built under R version", msg, fixed = TRUE) ||
    grepl("ABI version", msg, fixed = TRUE)
}

quiet_library <- function(pkg) {
  muffled <- character()
  res <- tryCatch(
    {
      withCallingHandlers(
        suppressPackageStartupMessages(
          library(pkg, character.only = TRUE, quietly = TRUE)
        ),
        warning = function(w) {
          msg <- conditionMessage(w)
          if (is_benign_build_warning(msg)) {
            muffled[length(muffled) + 1L] <<- msg
            invokeRestart("muffleWarning")
          }
        }
      )
      list(ok = TRUE, error = NULL)
    },
    error = function(e) list(ok = FALSE, error = conditionMessage(e))
  )
  res$muffled <- muffled
  res
}

source(file.path("R", "utils.R"))

load_status <- lapply(required_pkgs, quiet_library)
names(load_status) <- required_pkgs

failed <- names(load_status)[!vapply(load_status, `[[`, logical(1), "ok")]
if (length(failed) > 0L) {
  for (pkg in failed) {
    cli::cli_alert_danger("{.pkg {pkg}} failed to load: {load_status[[pkg]]$error}")
  }
  stop(sprintf("Package load failures: %s", paste(failed, collapse = ", ")),
       call. = FALSE)
}

cli::cli_alert_success("All {length(required_pkgs)} required packages loaded.")

muffled_counts <- vapply(load_status, function(s) length(s$muffled), integer(1))
if (any(muffled_counts > 0L)) {
  cli::cli_alert_info(
    "Muffled {sum(muffled_counts)} benign R-version build notice{?s} \\
    (see INSTALL_NOTES.md)."
  )
}

pkg_versions <- tibble::tibble(
  package = required_pkgs,
  version = vapply(required_pkgs,
                   function(p) as.character(utils::packageVersion(p)),
                   character(1)),
  built   = vapply(required_pkgs,
                   function(p) {
                     b <- utils::packageDescription(p)$Built
                     if (is.null(b)) NA_character_ else b
                   },
                   character(1))
)

cat("\n")
cli::cli_h2("Pinned package versions")
print(pkg_versions, n = nrow(pkg_versions))

cat("\n")
cli::cli_h2("sessionInfo()")
print(sessionInfo())

cli::cli_alert_success("Environment check complete.")
invisible(NULL)
