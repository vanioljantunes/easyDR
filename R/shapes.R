# Shape resolution and the design bases. The only file that decides what form the
# dose-response takes. Nothing here is exported.

#' Resolve a shape name into a shape object
#'
#' Returns the shape entity of the data model: its name, how many terms it
#' contributes, the dose levels a study needs to inform it, the parameters that
#' define it, and the basis function that builds its design columns. The basis
#' closes over the fitted parameters, so prediction later rebuilds the same
#' design rather than re-deriving knots from a prediction grid.
#'
#' @keywords internal
#' @noRd
dr_shape <- function(shape, dose, knots = NULL, powers = NULL, studies = NULL) {
  # An already-resolved shape passes straight through, keeping the parameters it
  # was fitted with. That is what lets a sensitivity refit carry the full fit's
  # knots instead of re-deriving them from a subset, which would silently change
  # the model between a fit and its own sensitivity analysis. The checks still
  # run against the subset, so a refit the shape no longer suits fails with a
  # reason rather than quietly fitting something else.
  if (dr_is_shape(shape)) {
    dr_check_shape(shape, dose, studies)
    return(shape)
  }
  shape <- match.arg(shape, c("linear", "rcs", "fp", "quadratic", "log"))
  out <- switch(shape,
                linear    = dr_shape_linear(),
                quadratic = dr_shape_quadratic(),
                log       = dr_shape_log(dose),
                rcs       = dr_shape_rcs(dose, knots),
                fp        = dr_shape_fp(dose, powers))
  dr_check_shape(out, dose, studies)
  out
}

#' Is this an already-resolved shape rather than a shape name?
#'
#' Recognised by structure rather than by class, so every constructor stays as
#' it is: a shape carries a name, a term count and the basis closure that holds
#' its fitted parameters.
#'
#' @keywords internal
#' @noRd
dr_is_shape <- function(x) {
  is.list(x) && !is.null(x$name) && !is.null(x$terms) && is.function(x$basis)
}

#' @keywords internal
#' @noRd
dr_shape_linear <- function() {
  list(name = "linear", terms = 1L, min_levels = 2L, knots = NULL,
       powers = NULL, scaling = NULL, label = "linear in dose",
       basis = function(d) matrix(d, ncol = 1L, dimnames = list(NULL, "dose")))
}

#' @keywords internal
#' @noRd
dr_shape_quadratic <- function() {
  list(name = "quadratic", terms = 2L, min_levels = 3L, knots = NULL,
       powers = NULL, scaling = NULL, label = "quadratic in dose",
       basis = function(d) {
         # The dose comes first so a Wald test on the remaining terms is a test
         # of departure from linearity (research decision D1).
         cbind(dose = d, dose_sq = d^2)
       })
}

#' @keywords internal
#' @noRd
dr_shape_log <- function(dose) {
  shift <- if (any(dose <= 0)) 1 else 0
  list(name = "log", terms = 1L, min_levels = 2L, knots = NULL, powers = NULL,
       scaling = if (shift > 0) list(shift = shift) else NULL,
       label = if (shift > 0) {
         paste0("logarithmic in dose, shifted by ", shift,
                " because a dose of zero is present")
       } else {
         "logarithmic in dose"
       },
       basis = function(d) matrix(log(d + shift), ncol = 1L,
                                  dimnames = list(NULL, "log_dose")))
}

#' Restricted cubic spline, Harrell's parameterization
#'
#' The first basis column is the dose itself and the rest are the nonlinear
#' components, which is what makes a Wald test on the non-first terms a test of
#' departure from linearity. A natural-spline basis spans the same space but
#' mixes the dose into every column, and testing one of those columns answers a
#' different question: on the bundled alcohol data the correct test gives
#' p = 0.052 while that substitution gives p = 3.8e-06 (research decision D1).
#'
#' Built with `Hmisc::rcspline.eval`, which is the exported routine behind
#' `rms::rcs` and returns the same basis. `rms::rcs` is formula-context
#' machinery and returns its input unchanged when called directly, so calling it
#' from a package would silently produce a one-column design.
#'
#' @keywords internal
#' @noRd
dr_shape_rcs <- function(dose, knots) {
  positive <- dose[dose > 0]
  if (length(positive) == 0L) positive <- dose

  k <- if (is.null(knots)) {
    as.numeric(stats::quantile(positive, c(0.10, 0.50, 0.90), names = FALSE))
  } else if (length(knots) == 1L) {
    n <- as.integer(knots)
    if (is.na(n) || n < 3L) {
      stop("easydr: a restricted cubic spline needs at least 3 knots. ",
           "`knots` was given as ", knots, ".", call. = FALSE)
    }
    as.numeric(stats::quantile(positive, seq(0.05, 0.95, length.out = n),
                               names = FALSE))
  } else {
    as.numeric(knots)
  }

  # V16: the knots must be distinct and inside the observed dose range.
  if (anyDuplicated(k) > 0L) {
    stop("easydr: the spline knots ", paste(signif(k, 4), collapse = ", "),
         " include repeated values. Knots must be distinct. The observed dose ",
         "range is ", signif(min(dose), 4), " to ", signif(max(dose), 4), ".",
         call. = FALSE)
  }
  outside <- k[k < min(dose) | k > max(dose)]
  if (length(outside) > 0L) {
    stop("easydr: the spline knots ", paste(signif(outside, 4), collapse = ", "),
         " fall outside the observed dose range, which is ",
         signif(min(dose), 4), " to ", signif(max(dose), 4),
         ". Place the knots inside the range, or give a knot count instead.",
         call. = FALSE)
  }
  if (length(k) < 3L) {
    stop("easydr: a restricted cubic spline needs at least 3 knots. ",
         length(k), " were given.", call. = FALSE)
  }

  list(name = "rcs", terms = length(k) - 1L, min_levels = length(k),
       knots = k, powers = NULL, scaling = NULL,
       label = paste0("restricted cubic spline, knots at ",
                      paste(signif(k, 4), collapse = ", ")),
       basis = function(d) {
         m <- as.matrix(Hmisc::rcspline.eval(d, knots = k, inclx = TRUE))
         colnames(m) <- c("dose", paste0("dose_nl", seq_len(ncol(m) - 1L)))
         m
       })
}

#' Fractional polynomial through the engine
#'
#' A negative power is undefined at a dose of zero, so the dose is shifted when
#' a zero is present and the shift is reported, because it changes how a power
#' reads.
#'
#' The engine's own scaling is switched off deliberately. Left on, it derives a
#' shift and a scale from whatever vector it is handed, so fitting on the study
#' doses and predicting on a grid would silently use two different transforms.
#' Fixing the shift here makes the basis a function of the fitted shape alone
#' (research decision D6, corrected during implementation).
#'
#' @keywords internal
#' @noRd
dr_shape_fp <- function(dose, powers) {
  p <- if (is.null(powers)) c(-1, 0) else as.numeric(powers)
  if (length(p) != 2L) {
    stop("easydr: `powers` must be a pair, for example c(-1, 0). ",
         length(p), " value", if (length(p) == 1L) " was" else "s were",
         " given.", call. = FALSE)
  }
  shift <- if (any(dose <= 0)) 1 - min(dose) else 0

  # A power of 3 on a dose reaching 60 puts design values near 2e5, which
  # ill-conditions the likelihood and makes one candidate look hundreds of AIC
  # better than the rest. Royston and Altman scale the dose for exactly this
  # reason. The scale is a fixed power of ten chosen from the fitting doses and
  # stored, so the basis stays a function of the fitted shape while the design
  # stays of order one.
  scale <- 10^round(log10(max(dose + shift, na.rm = TRUE)))
  if (!is.finite(scale) || scale <= 0) scale <- 1

  list(name = "fp", terms = 2L, min_levels = 3L, knots = NULL, powers = p,
       scaling = list(shift = shift, scale = scale),
       label = paste0("fractional polynomial, powers ",
                      paste(p, collapse = ", "),
                      if (shift > 0) {
                        paste0(", dose shifted by ", signif(shift, 4),
                               " because a dose of zero or below is present")
                      } else "",
                      if (scale != 1) {
                        paste0(", dose divided by ", format(scale, scientific = FALSE),
                               " to keep the design well conditioned")
                      } else ""),
       basis = function(d) {
         # fracpol needs at least two values to name its columns, so a single
         # dose is duplicated and the extra row dropped.
         one <- length(d) == 1L
         x   <- if (one) rep(d, 2L) else d
         m   <- as.matrix(dosresmeta::fracpol((x + shift) / scale, p = p,
                                              scaling = FALSE))
         if (one) m <- m[1L, , drop = FALSE]
         colnames(m) <- paste0("fp", seq_len(ncol(m)))
         m
       })
}

#' Rules V14, V15 and V17: can this shape be estimated from this data
#' @keywords internal
#' @noRd
dr_check_shape <- function(shape, dose, studies) {
  if (identical(shape$name, "linear")) return(invisible(TRUE))

  if (!is.null(studies)) {
    # V14: some study must carry enough distinct doses to inform the shape.
    levels_per_study <- tapply(dose, studies, function(x) length(unique(x)))
    best <- max(levels_per_study)
    if (best < shape$min_levels) {
      stop("easydr: ", shape$label, " needs a study with at least ",
           shape$min_levels, " distinct dose levels. The best any study offers ",
           "is ", best, ". Choose a simpler shape, or report the linear fit. ",
           "easydr does not substitute one shape for another.", call. = FALSE)
    }

    # V15: terms must leave degrees of freedom for heterogeneity and the
    # prediction interval.
    n_studies <- length(unique(studies))
    if (shape$terms >= n_studies) {
      stop("easydr: ", shape$label, " fits ", shape$terms, " terms, and there ",
           "are only ", n_studies, " stud", if (n_studies == 1L) "y" else "ies",
           ". The heterogeneity and the prediction interval would have no ",
           "degrees of freedom left.", call. = FALSE)
    }
  }

  # V17: the design must have full rank, which is the degenerate-basis case.
  B <- shape$basis(dose)
  if (any(!is.finite(B))) {
    stop("easydr: ", shape$label, " produces non-finite design values for the ",
         "observed doses. Choose different parameters.", call. = FALSE)
  }
  if (qr(B)$rank < ncol(B)) {
    stop("easydr: ", shape$label, " produces a rank-deficient design for the ",
         "observed doses, usually because they cluster too tightly for the ",
         "requested parameters. easydr stops rather than return an unstable ",
         "curve.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Contrast of each dose against the reference, in basis space
#'
#' The prediction interval and the curve both work on contrasts, so the
#' reference is subtracted in basis space rather than on the dose scale.
#'
#' @keywords internal
#' @noRd
dr_contrast <- function(shape, doses, xref) {
  B   <- shape$basis(c(xref, doses))
  ref <- B[1L, , drop = FALSE]
  out <- B[-1L, , drop = FALSE]
  sweep(out, 2L, as.numeric(ref), "-")
}
