#' Name the shape of a fitted dose-response curve
#'
#' Readers do not argue about coefficients; they argue about shape. This gives
#' the shape a name, in the vocabulary dose-response reviews actually publish,
#' so that a results paragraph can be written without inventing terms.
#'
#' @section Where the names come from:
#' Published dose-response meta-analyses trichotomise the curve: it either
#' plateaus after an initial rise, is still rising at the maximum dose, or is
#' bell-shaped, with the effect falling away at the top (Leucht et al. 2020).
#' Those three names are used here as published, with direction made explicit
#' for the monotone case, the mirror of a bell added for exposures whose effect
#' falls and then returns, and a name for a curve that does nothing.
#'
#' `"plateau"` is deliberately not called Emax, even though easyDR can now fit
#' an Emax model with `shape = "emax"`. The distinction is the point: Emax
#' names a functional form that was fitted, and `"plateau"` names a shape that
#' was observed in a curve. A spline that happens to flatten has not fitted an
#' Emax model, and labelling it Emax would claim a model nobody chose.
#'
#' @section What the labels are judged on:
#' The labels describe how the effect **changes** from one dose to the next,
#' not where the curve sits, so a monotone curve earns the same name on the log
#' scale and on the ratio scale. A plateau is judged against the curve's own
#' steepest step rather than against its overall range, which is what separates
#' a genuine flattening from a curve that merely decelerates.
#'
#' [easydr()] passes the log scale and carries the answer on the `curve_label`
#' slot of its result.
#'
#' @param dose Numeric vector of exposure values, ascending.
#' @param effect Numeric vector of fitted effects at those exposures, on the
#'   scale you intend to plot (log scale for ratio measures).
#' @param tol Fraction of the total effect range below which a single step
#'   counts as flat.
#' @param plateau_frac A curve is called a plateau when its closing steps are
#'   smaller than this fraction of its steepest step.
#'
#' @return One of `"monotonic increasing"`, `"monotonic decreasing"`,
#'   `"plateau"`, `"bell-shaped"`, `"U-shaped"` or `"flat"`.
#'
#' @references
#' Leucht S, Crippa A, Siafis S, Patel MX, Orsini N, Davis JM (2020).
#' Dose-response meta-analysis of antipsychotic drugs for acute schizophrenia.
#' American Journal of Psychiatry 177, 342-353.
#'
#' Crippa A, Orsini N (2016). Dose-response meta-analysis of differences in
#' means. BMC Medical Research Methodology 16, 91.
#'
#' @keywords internal
#' @noRd
curve_shape <- function(dose, effect, tol = 0.05, plateau_frac = 0.25) {
  stopifnot(length(dose) == length(effect), length(dose) >= 3)
  ord <- order(dose)
  effect <- effect[ord]

  span <- diff(range(effect))
  if (span == 0 || span < tol * max(abs(effect), 1)) {
    return("flat")
  }

  steps <- diff(effect)
  flat <- abs(steps) < tol * span
  up <- steps > 0 & !flat
  down <- steps < 0 & !flat

  peak <- which.max(effect)
  trough <- which.min(effect)
  last <- length(effect)

  # A turning point strictly inside the range, with real movement on both
  # sides, is what makes a curve non-monotone rather than merely noisy at one
  # end. Which way it turns is the difference between a bell and a U, and the
  # two mean opposite things, so they are not collapsed into one name.
  if (peak > 1 && peak < last && any(up[seq_len(peak - 1)]) &&
      any(down[seq(peak, last - 1)])) {
    return("bell-shaped")
  }
  if (trough > 1 && trough < last && any(down[seq_len(trough - 1)]) &&
      any(up[seq(trough, last - 1)])) {
    return("U-shaped")
  }

  # "Flattens out" is relative to how fast the curve was climbing, not to the
  # range it covers: a shallow curve and a steep one should earn the same name.
  steepest <- max(abs(steps))
  tail_steps <- steps[seq(max(1, last - 2), last - 1)]
  tail_flat <- steepest > 0 && all(abs(tail_steps) < plateau_frac * steepest)
  if (any(up) && !any(down)) {
    return(if (tail_flat) "plateau" else "monotonic increasing")
  }
  if (any(down) && !any(up)) {
    return(if (tail_flat) "plateau" else "monotonic decreasing")
  }
  if (effect[last] >= effect[1L]) "bell-shaped" else "U-shaped"
}
