# The prediction interval. This is the one quantity easydr computes rather than
# reads from the engine, because predict.dosresmeta() returns confidence limits
# only (research decision D3, constitution Principle II exception).

#' Pointwise prediction interval for the pooled dose-response
#'
#' The confidence interval describes the pooled slope. The prediction interval
#' describes where a new study's dose-response would fall, so it adds the
#' between-study variance to the uncertainty of the pooled estimate:
#'
#' \deqn{x'b \pm t_{1-\alpha/2,\,k-p} \sqrt{x'(V_b + \Psi)x}}
#'
#' with \eqn{x} the dose contrast against the reference dose, \eqn{V_b} the
#' covariance of the fixed effects, \eqn{\Psi} the between-study covariance, and
#' \eqn{k-p} studies minus fitted coefficients. This is the standard
#' random-effects predictive distribution of Higgins, Thompson and Spiegelhalter
#' (2009), applied to the dose contrast.
#'
#' Returned on the linear predictor scale. The caller back-transforms.
#'
#' @references
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009). A re-evaluation of
#' random-effects meta-analysis. Journal of the Royal Statistical Society Series
#' A 172, 137-159.
#'
#' @keywords internal
#' @noRd
dr_prediction_interval <- function(fit, doses, xref, level = 0.95,
                                   single_study = FALSE, shape = NULL) {
  b <- as.numeric(coef(fit))
  k <- length(unique(fit$id))
  p <- length(b)

  # A prediction interval needs between-study variation to describe, and at
  # least one degree of freedom left over. Say so rather than invent one. This
  # check comes before the covariance is coerced, because a single-study fit
  # carries no between-study covariance at all.
  reason <- if (single_study) {
    "one study only, so there is no between-study variation to describe"
  } else if (is.null(fit$Psi)) {
    paste0("this is a common-effect fit, which assumes no between-study ",
           "variation, so there is no predictive distribution to report. ",
           "Refit with random = TRUE for a prediction interval")
  } else if (k - p < 1L) {
    paste0("a prediction interval needs more studies than fitted ",
           "coefficients; this fit has ", k, " stud",
           if (k == 1) "y" else "ies", " and ", p, " coefficient",
           if (p == 1) "" else "s")
  } else {
    NA_character_
  }
  if (!is.na(reason)) {
    return(list(estimable = FALSE, reason = reason, dose = doses,
                lb = rep(NA_real_, length(doses)),
                ub = rep(NA_real_, length(doses))))
  }
  # The contrast against the reference dose, in basis space: a scalar for a
  # linear fit, a row per dose for a curve.
  X <- if (is.null(shape) || identical(shape$name, "linear")) {
    matrix(doses - xref, ncol = 1L)
  } else {
    dr_contrast(shape, doses, xref)
  }

  Vb <- as.matrix(vcov(fit))
  Psi <- as.matrix(fit$Psi)
  V     <- Vb + Psi
  est   <- as.numeric(X %*% b)
  tcrit <- stats::qt(1 - (1 - level) / 2, df = k - p)
  sepi  <- sqrt(pmax(0, rowSums((X %*% V) * X)))

  list(estimable = TRUE,
       dose = doses,
       est  = est,
       lb   = est - tcrit * sepi,
       ub   = est + tcrit * sepi,
       df   = k - p,
       psi  = Psi)
}
