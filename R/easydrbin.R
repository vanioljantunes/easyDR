# The binary worker behind easydr(). It keeps its own signature because
# that signature is what dr_refit() calls back into, with the column and model
# arguments recorded on the fitted result. Not exported: easydr() is the one
# public entry point.

#' Fit a binary dose-response
#' @keywords internal
#' @noRd
easydrbin <- function(data,
                      dose = NULL,
                      studlab = NULL,
                      sm = "RR",
                      event.e = NULL,
                      n.e = NULL,
                      TE = NULL,
                      lower = NULL,
                      upper = NULL,
                      ref = NULL,
                      design = NULL,
                      method.tau = "REML",
                      proc = "2stage",
                      covariance = NULL,
                      shape = "linear",
                      knots = NULL,
                      powers = NULL,
                      grid = NULL,
                      at = NULL,
                      dose.unit = NULL,
                      level = 0.95,
                      random = TRUE,
                      common = FALSE) {
  call <- match.call()
  sm <- match.arg(sm, c("RR", "OR", "IRR"))

  dr_entry("bin", data, call,
           roles_args = list(dose = dose, studlab = studlab,
                             event.e = event.e, n.e = n.e, TE = TE,
                             lower = lower, upper = upper, ref = ref,
                             design = design),
           sm = sm, method.tau = method.tau, proc = proc,
           covariance = covariance, method.smd = NULL,
           shape = shape, knots = knots, powers = powers, grid = grid,
           at = at, dose.unit = dose.unit, level = level, common = common,
           random = random)
}

#' Translate method.tau into the engine's estimator name
#' @keywords internal
#' @noRd
dr_method_tau <- function(method.tau, common, random = TRUE) {
  for (nm in c("common", "random")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("easydr: `", nm, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  # easydr fits one model, so this pair chooses it rather than choosing what to
  # display. `meta` can show a common-effect and a random-effects row from one
  # object; a dose-response fit estimates the between-study covariance or it
  # does not, so asking for the common effect by either route gives the
  # common-effect fit, and the print says which model was fitted.
  if (isTRUE(common) || !isTRUE(random)) return("fixed")
  switch(toupper(method.tau),
         "REML"  = "reml",
         "ML"    = "ml",
         "MM"    = "mm",
         "FIXED" = "fixed",
         stop("easydr: method.tau must be one of ",
              paste(dQuote(c("REML", "ML", "mm", "fixed")), collapse = ", "),
              ".", call. = FALSE))
}

#' Check the confidence level
#' @keywords internal
#' @noRd
dr_check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("easydr: level must be a single number between 0 and 1, for example ",
         "0.95.", call. = FALSE)
  }
  invisible(TRUE)
}
