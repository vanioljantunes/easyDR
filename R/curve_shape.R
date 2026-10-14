#' Name the shape of a fitted dose-response curve
#'
#' Readers do not argue about coefficients; they argue about shape. This gives
#' the shape a name, using the vocabulary the literature actually uses, so that
#' a results paragraph can be written without inventing terms.
#'
#' The names follow Crippa and Orsini: a curve that rises throughout is a
#' *positive association*, one that rises and then falls is *non-monotonic*,
#' and one that rises and then flattens is an *Emax* curve. A curve whose
#' confidence band contains the null across the whole exposure range is
#' reported here as *no association*, which the literature does not name
#' because there is nothing to name.
#'
#' @param dose Numeric vector of exposure values, ascending.
#' @param effect Numeric vector of fitted effects at those exposures, on the
#'   scale you intend to plot (log scale for ratio measures).
#' @param null_value The value meaning "no effect": 0 on a log scale, 1 on a
#'   ratio scale.
#' @param tol Fraction of the total effect range below which a single step
#'   counts as flat.
#' @param plateau_frac A curve is called Emax when its closing steps are
#'   smaller than this fraction of its steepest step. Judging the tail against
#'   the steepest rise, rather than against the overall range, is what makes
#'   "it flattens out" mean the same thing for a shallow curve and a steep one.
#'
#' @return A single string: "positive association", "inverse association",
#'   "non-monotonic", "Emax", or "no association".
#'
#' @examples
#' curve_shape(0:4, c(0, 0.2, 0.4, 0.6, 0.8))
#' curve_shape(0:4, c(0, 0.5, 0.8, 0.85, 0.86))
#' curve_shape(0:4, c(0, 0.5, 0.9, 0.4, -0.1))
#'
#' @export
curve_shape <- function(dose, effect, null_value = 0, tol = 0.05,
                        plateau_frac = 0.25) {
  stopifnot(length(dose) == length(effect), length(dose) >= 3)
  ord <- order(dose)
  effect <- effect[ord]

  span <- diff(range(effect))
  if (span == 0 || span < tol * max(abs(effect), 1)) {
    return("no association")
  }

  steps <- diff(effect)
  flat <- abs(steps) < tol * span
  up <- steps > 0 & !flat
  down <- steps < 0 & !flat

  peak <- which.max(effect)
  trough <- which.min(effect)
  last <- length(effect)

  # A turning point strictly inside the range, with real movement on both
  # sides, is the thing that makes a curve non-monotonic rather than merely
  # noisy at one end.
  if (peak > 1 && peak < last && any(up[seq_len(peak - 1)]) &&
      any(down[seq(peak, last - 1)])) {
    return("non-monotonic")
  }
  if (trough > 1 && trough < last && any(down[seq_len(trough - 1)]) &&
      any(up[seq(trough, last - 1)])) {
    return("non-monotonic")
  }

  # "Flattens out" is relative to how fast the curve was climbing, not to the
  # range it covers: a shallow curve and a steep one should earn the same name.
  steepest <- max(abs(steps))
  tail_steps <- steps[seq(max(1, last - 2), last - 1)]
  tail_flat <- steepest > 0 && all(abs(tail_steps) < plateau_frac * steepest)
  if (any(up) && !any(down)) {
    return(if (tail_flat) "Emax" else "positive association")
  }
  if (any(down) && !any(up)) {
    return("inverse association")
  }
  "non-monotonic"
}
