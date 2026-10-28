#' Check whether a set of studies can support a dose-response meta-analysis
#'
#' The commonest reason a dose-response meta-analysis fails is not the model;
#' it is that the included studies never could have supported one. [easydr()]
#' runs these checks on the prepared analysis frame and carries the answer on
#' its result, so the reviewer sees it beside the fit rather than having to
#' remember to ask for it.
#'
#' Two-stage methods need at least three exposure levels per study, because a
#' study with two points carries no information about curvature and the
#' within-study covariance cannot be reconstructed from it (Greenland and
#' Longnecker, 1992). One-stage methods keep such studies (Crippa et al.,
#' 2019), which is precisely their selling point, so a dataset that fails the
#' two-stage check is not useless: it is a dataset that argues for
#' `proc = "1stage"`.
#'
#' @param data A data frame, one row per exposure level per study.
#' @param study Column identifying the study.
#' @param dose Column holding the exposure level, in a common unit across
#'   studies.
#' @param ref Optional logical column marking each study's reference category.
#'   When given, that is what `has_reference` reports, which is the reference
#'   the fit actually used. When `NULL`, a dose of zero is taken as the
#'   reference, which is only a guess.
#' @param cases,n Columns holding the event count and the denominator. Either
#'   may be `NULL`, which skips it: a table of reported effects carries no
#'   counts by design and has nothing to complain about.
#' @param min_levels Minimum exposure levels per study for a two-stage fit.
#'   Three is the conventional floor and the default.
#'
#' @return An object of class `dr_viability`: a list with `ok`, a per-study
#'   `summary` data frame, and a character vector of `problems`.
#'
#' @export
dr_viability <- function(data, study = "study", dose = "dose", ref = NULL,
                         cases = "cases", n = "n", min_levels = 3L) {
  stopifnot(is.data.frame(data))
  for (col in c(study, dose, ref)) {
    if (!col %in% names(data)) {
      stop("column '", col, "' not found in data", call. = FALSE)
    }
  }

  studies <- unique(data[[study]])
  rows_of <- function(s) which(data[[study]] == s)
  levels_per_study <- vapply(
    studies,
    function(s) length(unique(data[[dose]][rows_of(s)])),
    integer(1)
  )
  has_reference <- vapply(
    studies,
    function(s) {
      i <- rows_of(s)
      if (is.null(ref)) {
        any(data[[dose]][i] == 0, na.rm = TRUE)
      } else {
        any(as.logical(data[[ref]][i]), na.rm = TRUE)
      }
    },
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
      "study(ies) without a reference level (",
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
    list(ok = length(problems) == 0L, summary = summary, problems = problems,
         min_levels = min_levels),
    class = "dr_viability"
  )
}

#' Print a viability check
#'
#' @param x a `dr_viability` object, carried on the `viability` slot of an
#'   [easydr()] result.
#' @param ... ignored.
#'
#' @return `x`, invisibly. Called for the printed output.
#'
#' @export
print.dr_viability <- function(x, ...) {
  cat("Dose-response viability\n\n")
  print(x$summary, row.names = FALSE, right = FALSE)
  if (x$ok) {
    cat("\nNo blocking problems found.\n")
  } else {
    cat("\nProblems:\n")
    for (p in x$problems) cat(" - ", p, "\n", sep = "")
  }
  invisible(x)
}

#' One line summarising the viability check, for the result print
#'
#' The per-study table is on the result for anyone who wants it. What the
#' printed fit needs is the one sentence that says whether it should have been
#' fitted at all.
#'
#' @keywords internal
#' @noRd
dr_viability_line <- function(v) {
  s  <- v$summary
  ok <- sum(s$two_stage_ok)
  nr <- sum(!s$has_reference)
  out <- paste0("Viability: ", ok, " of ", nrow(s), " studies carry ",
                v$min_levels, " or more dose levels")
  if (nr > 0L) {
    out <- paste0(out, "; ", nr, " without a reference category")
  }
  if (!v$ok) out <- paste0(out, ". See $viability for what to do about it")
  paste0(out, ".")
}
