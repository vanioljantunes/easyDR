# Several shapes fitted to one table, held together.
#
# Comparing shapes is a normal part of a dose-response review: the protocol
# names a primary model and a secondary one, and the paper reports both. Doing
# that by hand means a list of fits, a loop and a layout, which is what the
# book used to show. This holds the fits instead, keeps them in the order they
# were named, and lets easydr_curve() draw them as panels.
#
# It is deliberately not an easydr_result. A forest, a funnel, a leave-one-out
# and a table each answer a question about one model, and quietly answering it
# for whichever shape happened to be named first would be worse than refusing.

#' Hold one fit per shape
#' @keywords internal
#' @noRd
dr_multi <- function(fits, shapes, call) {
  first <- fits[[1L]]
  structure(
    list(fits      = fits,
         shapes    = shapes,
         sm        = first$sm,
         sm_label  = first$sm_label,
         ratio     = first$ratio,
         dose_unit = first$dose_unit,
         reference = first$reference,
         counts    = first$counts,
         call      = call),
    class = "easydr_multi")
}

#' Refuse a multi fit, and say how to reach one of its fits
#'
#' The message names the shapes it holds and shows the expression that gets
#' one of them, because "takes a fitted result" on its own leaves the reader
#' with an object and no obvious way to open it.
#'
#' @keywords internal
#' @noRd
dr_require_single <- function(x, what) {
  if (inherits(x, "easydr_multi")) {
    stop("easydr: ", what, " reports one model, and this holds ",
         length(x$fits), " (", paste(x$shapes, collapse = ", "), "). ",
         "Pick the one you mean, such as x$fits$", x$shapes[1L],
         ". To see them side by side use ",
         "easydr_curve(x, mult.curves = TRUE).", call. = FALSE)
  }
  if (!inherits(x, "easydr_result")) {
    stop("easydr: ", what, " takes a fitted result from easydr().",
         call. = FALSE)
  }
  invisible(x)
}

#' A number, or NA when the engine will not give one
#' @keywords internal
#' @noRd
dr_try_num <- function(expr) {
  v <- try(suppressWarnings(as.numeric(expr)), silent = TRUE)
  if (inherits(v, "try-error") || length(v) != 1L || !is.finite(v)) {
    return(NA_real_)
  }
  v
}

#' One row per shape: what each model costs and what it buys
#'
#' The columns are the ones a protocol compares shapes on. Terms is how many
#' coefficients the shape spends. Departure from linearity is the Wald test on
#' the non-linear terms, which a one-term shape has nothing to report for.
#' Deviance is the fit against the saturated model. AIC and BIC are the fit
#' with the terms charged for, which is what lets shapes that are not nested
#' inside one another be compared at all.
#'
#' @keywords internal
#' @noRd
dr_multi_table <- function(x, digits = 3) {
  rows <- lapply(x$shapes, function(s) {
    f  <- x$fits[[s]]
    nl <- f$nonlinearity
    data.frame(
      Shape = s,
      Terms = f$shape$terms,
      `p nonlinearity` = if (isTRUE(nl$estimable)) dr_p(nl$p) else "-",
      Deviance = if (is.null(f$gof$deviance)) NA_real_ else
        round(f$gof$deviance, 2),
      `R-squared` = if (is.null(f$gof$r2)) NA_real_ else round(f$gof$r2, 3),
      AIC = round(dr_try_num(stats::AIC(f$fit)), 1),
      BIC = round(dr_try_num(stats::BIC(f$fit)), 1),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Print several shapes fitted to the same table
#'
#' @param x an object from [easydr()] called with more than one `shape`.
#' @param digits significant digits for reported effects.
#' @param ... ignored.
#'
#' @return `x`, invisibly. Called for the printed output.
#'
#' @examples
#' data(alcohol_cancer)
#' easydr.object <- easydr(alcohol_cancer, sm = "RR",
#'                         shape = c("linear", "quadratic", "rcs"),
#'                         dose.unit = "g/day")
#' easydr.object
#'
#' # One of them, to go on with.
#' easydr.object$fits$rcs
#'
#' @export
print.easydr_multi <- function(x, digits = 3, ...) {
  cat(length(x$fits), " shapes fitted to the same table (",
      x$sm_label, ")\n", sep = "")
  unit <- if (is.null(x$dose_unit)) "dose unit" else x$dose_unit
  cat(x$counts$studies, " studies, ", x$counts$categories,
      " dose categories, dose range ", x$counts$dose_min, " to ",
      x$counts$dose_max, " ", unit, "\n\n", sep = "")

  print(dr_multi_table(x, digits), row.names = FALSE, right = FALSE)

  cat("\nAIC and BIC charge for the extra terms, which is what lets shapes",
      "\nthat are not nested inside one another be compared. The lowest is",
      "\nthe shape the data supports, not the shape that fits closest.\n")
  cat("\nOne fit: x$fits$", x$shapes[1L], "\n", sep = "")
  cat("Side by side: easydr_curve(x, mult.curves = TRUE)\n")
  invisible(x)
}
