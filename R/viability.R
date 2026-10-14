#' Check whether a dataset can support a dose-response meta-analysis
#'
#' The commonest reason a dose-response meta-analysis fails is not the model;
#' it is that the included studies never could have supported one. This runs
#' the checks that matter before any fitting happens.
#'
#' Two-stage methods need at least three exposure levels per study, because a
#' study with two points carries no information about curvature and the
#' within-study covariance cannot be reconstructed from it (Greenland and
#' Longnecker, 1992). One-stage methods keep such studies (Crippa et al.,
#' 2019), which is precisely their selling point, so a dataset that fails the
#' two-stage check is not useless: it is a dataset that argues for one-stage.
#'
#' @param data A data frame, one row per exposure level per study.
#' @param study Column identifying the study.
#' @param dose Column holding the exposure level, already converted to a common
#'   unit across studies.
#' @param cases,n Columns holding the event count and the group size. Optional;
#'   when supplied they are checked for completeness, because the covariance
#'   reconstruction needs them.
#' @param min_levels Minimum exposure levels per study for a two-stage fit.
#'   Three is the conventional floor and the default.
#'
#' @return An object of class `dr_viability`: a list with `ok`, a per-study
#'   `summary` data frame, and a character vector of `problems`.
#'
#' @examples
#' d <- data.frame(
#'   study = c("A", "A", "A", "B", "B"),
#'   dose  = c(0, 10, 20, 0, 15),
#'   cases = c(10, 20, 31, 5, 9),
#'   n     = c(100, 100, 100, 50, 50)
#' )
#' dr_viability(d)
#'
#' @export
dr_viability <- function(data, study = "study", dose = "dose",
                         cases = "cases", n = "n", min_levels = 3L) {
  stopifnot(is.data.frame(data))
  for (col in c(study, dose)) {
    if (!col %in% names(data)) {
      stop("column '", col, "' not found in data", call. = FALSE)
    }
  }

  studies <- unique(data[[study]])
  levels_per_study <- vapply(
    studies,
    function(s) length(unique(data[[dose]][data[[study]] == s])),
    integer(1)
  )
  has_reference <- vapply(
    studies,
    function(s) any(data[[dose]][data[[study]] == s] == 0, na.rm = TRUE),
    logical(1)
  )

  summary <- data.frame(
    study = studies,
    levels = levels_per_study,
    has_reference = has_reference,
    two_stage_ok = levels_per_study >= min_levels,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  problems <- character(0)
  thin <- summary$study[!summary$two_stage_ok]
  if (length(thin)) {
    problems <- c(problems, paste0(
      length(thin), " study(ies) with fewer than ", min_levels,
      " exposure levels (", paste(thin, collapse = ", "),
      "): a two-stage fit drops these. Consider one-stage."
    ))
  }
  if (any(!summary$has_reference)) {
    problems <- c(problems, paste0(
      "study(ies) without a zero-dose reference level (",
      paste(summary$study[!summary$has_reference], collapse = ", "),
      "): check the reference category before pooling."
    ))
  }
  for (col in c(cases, n)) {
    if (col %in% names(data) && anyNA(data[[col]])) {
      problems <- c(problems, paste0(
        "missing values in '", col,
        "': the within-study covariance cannot be reconstructed without it."
      ))
    }
  }

  structure(
    list(ok = length(problems) == 0L, summary = summary, problems = problems),
    class = "dr_viability"
  )
}

#' @export
print.dr_viability <- function(x, ...) {
  cat("Dose-response viability\n\n")
  print(x$summary)
  if (x$ok) {
    cat("\nNo blocking problems found.\n")
  } else {
    cat("\nProblems:\n")
    for (p in x$problems) cat(" - ", p, "\n", sep = "")
  }
  invisible(x)
}
