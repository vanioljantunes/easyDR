# The parametric bootstrap.
#
# Two things in this package cannot be given an interval in closed form.
#
# The Emax family profiles a parameter: ED50 is chosen by fitting over a grid
# and keeping the best. The coefficient's usual interval then conditions on
# that choice, as though ED50 had been known in advance, and it is too narrow
# by however much the data had to say about where the curve bends.
#
# A target dose is a non-linear function of the coefficients with no standard
# error of its own, and for a spline read off a grid there is no derivative
# for the delta method to use.
#
# Both are answered the same way. Every study's reported contrasts are drawn
# again from the covariance easyDR reconstructed for them, the whole fit is
# redone on the draw, including the profile, and the quantity is read off each
# time. The spread of those values is the interval.
#
# Drawing the contrasts rather than the pooled coefficients is what lets the
# profile run again: the basis of an Emax curve depends on ED50, so a draw of
# the coefficients would be a draw conditional on the very parameter whose
# uncertainty is the point. `Slist` holds the covariance of the reported
# contrasts, which does not depend on the basis, so a draw from it can be
# refitted under any shape.

#' Draw one study's contrasts from their reconstructed covariance
#' @keywords internal
#' @noRd
dr_mvn_draw <- function(mu, sigma) {
  p <- length(mu)
  # A reconstructed covariance can be singular, or very slightly asymmetric
  # from the arithmetic. The eigen decomposition handles both, and a negative
  # eigenvalue is clamped rather than failing the draw.
  s <- (sigma + t(sigma)) / 2
  e <- eigen(s, symmetric = TRUE)
  v <- pmax(e$values, 0)
  as.numeric(mu + e$vectors %*% (sqrt(v) * stats::rnorm(p)))
}

#' The design matrix a shape puts on a set of doses
#' @keywords internal
#' @noRd
dr_design <- function(shape, dose) {
  if (is.null(shape) || identical(shape$name, "linear")) {
    return(cbind(dose = dose))
  }
  as.matrix(shape$basis(dose))
}

#' Refit the whole model on a parametric draw of the reported contrasts
#'
#' The draw is taken from the fitted model, not from the observed contrasts.
#' Resampling around what was observed would add a second helping of the
#' within-study error the data already carries one realisation of, which makes
#' every interval too wide and, once the weights are re-estimated on the
#' noisier draw, shifts the centre as well. Both were visible against the
#' closed-form interval for a straight line, which is how this was caught.
#'
#' So each study's contrasts are drawn around its fitted values, with the
#' covariance the model says they have: the reconstructed within-study
#' covariance, plus the between-study variance carried through the design.
#'
#' Returns `NULL` for a draw the engine cannot fit, which a grid search over a
#' perturbed table is expected to produce occasionally. The caller counts how
#' many succeeded and refuses to report an interval built on too few.
#'
#' @keywords internal
#' @noRd
dr_boot_once <- function(x, reprofile) {
  frame <- x$frame
  Slist <- x$fit$Slist
  ids   <- unique(frame$id)
  b     <- as.numeric(stats::coef(x$fit))
  Psi   <- x$Psi
  if (is.null(Psi)) Psi <- matrix(0, length(b), length(b))
  Psi   <- as.matrix(Psi)

  for (k in seq_along(ids)) {
    rows <- which(frame$id == ids[k] & !frame$is_ref)
    S    <- if (k <= length(Slist)) Slist[[k]] else NULL
    if (is.null(S) || !length(rows)) next
    S <- as.matrix(S)
    if (length(rows) != nrow(S)) next
    # `y` holds what the study reported on its own scale, an arm mean or a
    # log effect, not a contrast. The engine derives the contrasts against
    # the study's own reference row. So the draw is built as that reference
    # value plus the contrast the model predicts: dropping the baseline
    # biases every draw towards zero, which the straight-line check caught.
    ref <- which(frame$id == ids[k] & frame$is_ref)
    if (!length(ref)) next
    base <- frame$y[ref[1L]]
    xr <- as.numeric(dr_design(x$shape, frame$dose[ref[1L]]))
    D  <- dr_design(x$shape, frame$dose[rows])
    D  <- sweep(D, 2L, xr, "-")
    mu <- base + as.numeric(D %*% b)
    frame$y[rows] <- dr_mvn_draw(mu, S + D %*% Psi %*% t(D))
  }

  shape <- if (isTRUE(reprofile) && !is.null(x$shape$profiled_params)) {
    # The unresolved shape, so the profile runs again on this draw. Holding
    # ED50 at the value the real data chose is exactly the narrowness this
    # bootstrap exists to correct.
    try(dr_shape(x$shape$name, frame$dose), silent = TRUE)
  } else {
    x$shape
  }
  if (inherits(shape, "try-error")) return(NULL)

  p <- x$provenance
  shape <- try(dr_profile_shape(shape, frame, sm = x$sm, intake = p$intake,
                                proc = p$proc, method = p$method,
                                covariance = p$covariance,
                                method.smd = p$method.smd), silent = TRUE)
  if (inherits(shape, "try-error")) return(NULL)

  fit <- try(suppressWarnings(
    dr_fit(frame, sm = x$sm, intake = p$intake, proc = p$proc,
           method = p$method, covariance = p$covariance,
           method.smd = p$method.smd, shape = shape)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  list(fit = fit, shape = shape)
}

#' The effect at a set of doses, against the reference, for one draw
#' @keywords internal
#' @noRd
dr_boot_predict <- function(d, doses, xref) {
  b <- as.numeric(stats::coef(d$fit))
  vapply(doses, function(dose) {
    as.numeric(dr_contrast(d$shape, dose, xref) %*% b)
  }, numeric(1))
}

#' @keywords internal
#' @noRd
dr_check_boot <- function(B) {
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 50) {
    stop("easydr: `B` is the number of bootstrap draws and must be at least ",
         "50. Several hundred is usual, because the interval is read off ",
         "their percentiles and too few makes it noise.", call. = FALSE)
  }
  invisible(NULL)
}

#' The fitted effect at a grid of doses, once per bootstrap draw
#'
#' @param x a fitted `easydr_result`.
#' @param doses the doses to read the curve at.
#' @param B how many draws.
#' @param reprofile re-run the profile on each draw, which is the whole point
#'   for a shape that has one.
#' @param seed a seed, so a reported interval can be reproduced.
#' @param quiet `FALSE` reports progress, since several hundred refits of a
#'   profiled shape is not quick.
#'
#' @return a list holding the matrix of draws, one row per usable draw and one
#'   column per dose, on the model scale.
#'
#' @keywords internal
#' @noRd
dr_boot_curve <- function(x, doses, B = 500L, reprofile = TRUE, seed = NULL,
                          quiet = TRUE) {
  dr_check_boot(B)
  if (!is.null(seed)) set.seed(seed)
  xref <- x$reference$curve
  out  <- matrix(NA_real_, nrow = B, ncol = length(doses))
  for (b in seq_len(B)) {
    d <- dr_boot_once(x, reprofile)
    if (is.null(d)) next
    pred <- try(dr_boot_predict(d, doses, xref), silent = TRUE)
    if (inherits(pred, "try-error") || !all(is.finite(pred))) next
    out[b, ] <- pred
    if (!isTRUE(quiet) && b %% 50L == 0L) {
      message("easydr: bootstrap draw ", b, " of ", B)
    }
  }
  out <- out[stats::complete.cases(out), , drop = FALSE]
  if (nrow(out) < 20L) {
    stop("easydr: only ", nrow(out), " of ", B, " bootstrap draws could be ",
         "fitted, which is too few to read an interval from. This usually ",
         "means the table is too small for the shape.", call. = FALSE)
  }
  list(draws = out, n = nrow(out), B = B, doses = doses)
}

#' Percentile limits of a bootstrap distribution
#' @keywords internal
#' @noRd
dr_boot_limits <- function(draws, level = 0.95) {
  a <- (1 - level) / 2
  lim <- apply(draws, 2L, stats::quantile, probs = c(a, 1 - a), names = FALSE)
  list(lb = lim[1L, ], ub = lim[2L, ])
}
