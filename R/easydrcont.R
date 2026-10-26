# The continuous worker behind easydr(). It keeps its own signature because
# that signature is what dr_refit() calls back into, with the column and model
# arguments recorded on the fitted result. Not exported: easydr() is the one
# public entry point.

#' Fit a continuous dose-response
#' @keywords internal
#' @noRd
easydrcont <- function(data,
                       dose = NULL,
                       studlab = NULL,
                       sm = "MD",
                       mean.e = NULL,
                       sd.e = NULL,
                       n.e = NULL,
                       ref = NULL,
                       method.tau = "REML",
                       proc = "2stage",
                       method.smd = "cohens",
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
  sm <- match.arg(sm, c("MD", "SMD"))
  method.smd <- match.arg(method.smd, c("cohens", "hedges", "glass"))

  dr_entry("cont", data, call,
           roles_args = list(dose = dose, studlab = studlab,
                             mean.e = mean.e, sd.e = sd.e, n.e = n.e,
                             ref = ref),
           sm = sm, method.tau = method.tau, proc = proc,
           covariance = NULL, method.smd = method.smd,
           shape = shape, knots = knots, powers = powers, grid = grid,
           at = at, dose.unit = dose.unit, level = level, common = common,
           random = random)
}
