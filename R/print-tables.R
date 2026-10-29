# What the plotting functions hand back, and how it prints.
#
# A reader who wants the numbers behind a picture should be able to print them.
# Before this, every one of those tables came back as a bare data.frame of full
# double precision, so the walkthrough wrapped each in round() and rebuilt the
# columns it wanted by hand. The values stay exact; only the printing rounds.

#' Per-study contributions behind a forest plot
#'
#' The table [easydr_forest()] returns, printed the way the figure reads it:
#' the study, what it brought, its own trend with a standard error, and its
#' weight. The stored values keep full precision, so `as.data.frame()` gives
#' the untouched numbers.
#'
#' @param x the `table` element of what [easydr_forest()] returned.
#' @param digits significant digits for the trend and its standard error.
#' @param ... ignored.
#'
#' @return `x`, invisibly. Called for the printed output.
#'
#' @examples
#' easydr.object <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
#'                         dose.unit = "g/day")
#' f <- easydr_forest(easydr.object)
#' f$table
#'
#' @exportS3Method print easydr_contributions
print.easydr_contributions <- function(x, digits = 3, ...) {
  tab <- data.frame(
    study  = x$study,
    levels = x$levels,
    `dose range` = paste0(signif(x$dose_min, 3), " to ", signif(x$dose_max, 3)),
    trend  = signif(x$est, digits + 2),
    SE     = signif(x$se, digits + 2),
    `weight %` = round(x$weight, 2),
    check.names = FALSE, stringsAsFactors = FALSE)
  print(tab, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' Egger's test for small-study effects
#'
#' The `egger` element of what [easydr_funnel()] returns. Printing it writes
#' the intercept with its standard error, the test and the warning about power,
#' so nothing has to be unlisted and rounded by hand.
#'
#' @param x the `egger` element of what [easydr_funnel()] returned.
#' @param digits digits for the intercept and its standard error.
#' @param ... ignored.
#'
#' @return `x`, invisibly. Called for the printed output.
#'
#' @examples
#' easydr.object <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
#'                         dose.unit = "g/day")
#' easydr_funnel(easydr.object)$egger
#'
#' @exportS3Method print easydr_egger
print.easydr_egger <- function(x, digits = 4, ...) {
  cat("Egger's test for small-study effects\n")
  if (!isTRUE(x$estimable)) {
    cat("  not estimable: ", x$reason, "\n", sep = "")
    return(invisible(x))
  }
  cat("  intercept ", format(round(x$intercept, digits), nsmall = digits),
      " (SE ", format(round(x$se, digits), nsmall = digits), ")\n", sep = "")
  cat("  t = ", format(round(x$t, 2), nsmall = 2), " on ", x$df, " df, ",
      dr_p_lab(x$p), "\n", sep = "")
  cat("  ", x$k, " studies", sep = "")
  if (isTRUE(x$underpowered)) {
    cat(": fewer than ", x$k_warn, ", so the test has little power and a null ",
        "result means very little", sep = "")
  }
  cat("\n")
  invisible(x)
}
