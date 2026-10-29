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
dr_shape <- function(shape, dose, knots = NULL, powers = NULL,
                     ed50 = NULL, steepness = NULL, studies = NULL) {
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
  shape <- match.arg(shape, c("linear", "rcs", "fp", "quadratic", "cubic",
                              "log", "emax", "hyperbolic_emax"))
  out <- switch(shape,
                linear    = dr_shape_linear(),
                quadratic = dr_shape_quadratic(),
                cubic     = dr_shape_cubic(dose),
                log       = dr_shape_log(dose),
                rcs       = dr_shape_rcs(dose, knots),
                fp        = dr_shape_fp(dose, powers),
                emax             = dr_shape_emax(dose, ed50, steepness),
                hyperbolic_emax  = dr_shape_hyperbolic_emax(dose, ed50))
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
dr_shape_cubic <- function(dose) {
  # A dose reaching 60 gives a cubed design value near 2e5. That ill-conditions
  # the likelihood, which is the problem Royston and Altman scale a fractional
  # polynomial for, so the cubic is scaled the same way and by the same rule.
  # The scale is fixed from the fitting doses and stored, so the basis stays a
  # function of the fitted shape rather than of whatever grid it is handed.
  scale <- 10^round(log10(max(abs(dose), na.rm = TRUE)))
  if (!is.finite(scale) || scale <= 0) scale <- 1

  list(name = "cubic", terms = 3L, min_levels = 4L, knots = NULL,
       powers = NULL, scaling = list(shift = 0, scale = scale),
       label = paste0("cubic in dose",
                      if (scale != 1) {
                        paste0(", dose divided by ",
                               format(scale, scientific = FALSE),
                               " to keep the design well conditioned")
                      } else ""),
       basis = function(d) {
         # The dose comes first so a Wald test on the remaining terms is a test
         # of departure from linearity (research decision D1).
         z <- d / scale
         cbind(dose = z, dose_sq = z^2, dose_cu = z^3)
       })
}

#' @keywords internal
#' @noRd
dr_shape_log <- function(dose) {
  shift <- if (any(dose <= 0)) 1 else 0
  list(name = "log", terms = 1L, min_levels = 2L, knots = NULL, powers = NULL,
       # The coefficient of a one-term shape is per unit of whatever that term
       # is. Only the straight line has a term that is the dose itself, so only
       # it may be reported per unit of dose.
       slope_label = "Per 1 unit increase in log dose",
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
  # One power is a first-order fractional polynomial, which Royston and Altman
  # treat as the first candidate, and two is a second-order one. Three or more
  # is a model this package does not fit, and a silent truncation to the first
  # two would answer a different question.
  if (!length(p) %in% c(1L, 2L)) {
    stop("easydr: `powers` takes one value for a first-order fractional ",
         "polynomial, for example powers = -0.5, or a pair for a ",
         "second-order one, for example c(-1, 0). ", length(p),
         " values were given.", call. = FALSE)
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

  list(name = "fp", terms = length(p), min_levels = length(p) + 1L,
       knots = NULL, powers = p,
       slope_label = if (length(p) == 1L) {
         paste0("Per 1 unit increase in ",
                if (p[1L] == 0) "log(dose)" else paste0("dose^", p[1L]))
       } else NULL,
       scaling = list(shift = shift, scale = scale),
       label = paste0("fractional polynomial, ",
                      if (length(p) == 1L) "power " else "powers ",
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
         z <- (d + shift) / scale
         # dosresmeta::fracpol() computes second-order fractional polynomials
         # only, and says so, so a first-order one is built here from the
         # Royston and Altman definition: a single term z^p, with p = 0 read
         # as log(z).
         if (length(p) == 1L) {
           m <- matrix(if (p[1L] == 0) log(z) else z^p[1L], ncol = 1L)
           colnames(m) <- "fp1"
           return(m)
         }
         # fracpol needs at least two values to name its columns, so a single
         # dose is duplicated and the extra row dropped.
         one <- length(z) == 1L
         x   <- if (one) rep(z, 2L) else z
         m   <- as.matrix(dosresmeta::fracpol(x, p = p, scaling = FALSE))
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

# ---------------------------------------------------------------------------
# The Emax family
#
# Emax is not linear in its parameters, and the engine fits models that are.
# What rescues it is that for a fixed ED50 the curve is linear in the
# transformed dose d/(d + ED50): one coefficient, Emax, and an ordinary GLS.
# So ED50 is profiled, the same way Royston and Altman profile a fractional
# polynomial power, by fitting over a grid and keeping the best.
#
# Because effects in a dose-response meta-analysis are contrasts against a
# reference dose, E(reference) is zero by construction and the usual E0 term
# is absorbed. That leaves two parameters for the hyperbolic form, Emax and
# ED50, and three when the steepness is free.
#
# The profiled parameter is estimated from the same data as the coefficient,
# so a confidence interval that conditions on it is too narrow. That is why
# the Emax family carries a bootstrap, which re-profiles on every draw.
# ---------------------------------------------------------------------------

#' The hyperbolic Emax curve, steepness fixed at one
#'
#' `E(d) = Emax * d / (d + ED50)`.
#'
#' @keywords internal
#' @noRd
dr_shape_hyperbolic_emax <- function(dose, ed50 = NULL) {
  dr_emax_shape(dose, ed50 = ed50, steepness = 1, name = "hyperbolic_emax")
}

#' The Emax curve with a free steepness
#'
#' `E(d) = Emax * d^h / (d^h + ED50^h)`.
#'
#' @keywords internal
#' @noRd
dr_shape_emax <- function(dose, ed50 = NULL, steepness = NULL) {
  dr_emax_shape(dose, ed50 = ed50, steepness = steepness, name = "emax")
}

#' Both members of the family, since they differ only in what is profiled
#' @keywords internal
#' @noRd
dr_emax_shape <- function(dose, ed50, steepness, name) {
  obs <- dose[is.finite(dose) & dose > 0]
  if (length(obs) < 2L) {
    stop("easydr: the Emax family needs doses above zero, and this table ",
         "has ", length(obs), ".", call. = FALSE)
  }
  hyperbolic <- identical(name, "hyperbolic_emax")

  make <- function(theta, h) {
    theta <- as.numeric(theta)
    h     <- as.numeric(h)
    if (!is.finite(theta) || theta <= 0) {
      stop("easydr: `ed50` must be a dose above zero.", call. = FALSE)
    }
    if (!is.finite(h) || h <= 0) {
      stop("easydr: `steepness` must be above zero.", call. = FALSE)
    }
    lab <- if (hyperbolic) {
      paste0("hyperbolic Emax, ED50 = ", signif(theta, 4))
    } else {
      paste0("Emax, ED50 = ", signif(theta, 4), ", steepness = ",
             signif(h, 3))
    }
    list(name = name, terms = 1L, min_levels = 3L, knots = NULL,
         powers = NULL, scaling = NULL,
         ed50 = theta, steepness = h,
         # Emax itself, which is the coefficient, not a per-unit slope. The
         # curve has no constant slope, so none is claimed.
         slope_label = NULL,
         # The parameters the profile spent, which AIC and BIC have to charge
         # for: the engine only knows about the one coefficient it fitted.
         profiled_params = if (hyperbolic) 1L else 2L,
         label = lab,
         basis = function(d) {
           z <- ifelse(d > 0, d^h / (d^h + theta^h), 0)
           m <- cbind(z)
           colnames(m) <- "emax"
           m
         })
  }

  given_theta <- !is.null(ed50)
  given_h     <- hyperbolic || !is.null(steepness)
  if (given_theta && given_h) {
    return(make(ed50, if (hyperbolic) 1 else steepness))
  }

  # The grid. ED50 is a dose, so the candidates are doses: percentiles of the
  # observed range rather than a fixed sequence, which keeps the grid on the
  # scale the data is on whatever the units are. The steepness candidates are
  # the ones the dose-finding literature uses.
  theta_grid <- if (given_theta) {
    as.numeric(ed50)
  } else {
    # ED50 is a scale parameter, so the grid is spaced on the log scale and
    # reaches well below the smallest observed dose and above the largest.
    # A grid that only spans the observed doses pins ED50 to its own edge
    # whenever the data wants a curve that saturates immediately, and returns
    # that edge as though it were an estimate.
    lo <- min(obs) / 10
    hi <- max(obs) * 2
    exp(seq(log(lo), log(hi), length.out = 25L))
  }
  h_grid <- if (hyperbolic) {
    1
  } else if (!is.null(steepness)) {
    as.numeric(steepness)
  } else {
    c(0.5, 0.75, 1, 1.5, 2, 3)
  }

  out <- make(stats::median(theta_grid), h_grid[ceiling(length(h_grid) / 2)])
  out$profile    <- TRUE
  out$candidates <- expand.grid(ed50 = theta_grid, steepness = h_grid,
                                KEEP.OUT.ATTRS = FALSE)
  out$respecify  <- function(row) make(row$ed50, row$steepness)
  out
}

#' Fit the shape once per candidate and keep the best
#'
#' Chosen on the log likelihood, which for a fixed number of coefficients
#' orders the candidates the same way the deviance does. A candidate the
#' engine cannot fit is skipped rather than failing the call: a grid is
#' expected to contain values the data will not support.
#'
#' @keywords internal
#' @noRd
dr_profile_shape <- function(shape, frame, sm, intake, proc, method,
                             covariance, method.smd) {
  if (!isTRUE(shape$profile)) return(shape)
  cand <- shape$candidates
  best <- NULL
  best_ll <- -Inf
  for (i in seq_len(nrow(cand))) {
    s <- try(shape$respecify(cand[i, , drop = FALSE]), silent = TRUE)
    if (inherits(s, "try-error")) next
    f <- try(suppressWarnings(
      dr_fit(frame, sm = sm, intake = intake, proc = proc, method = method,
             covariance = covariance, method.smd = method.smd, shape = s)),
      silent = TRUE)
    if (inherits(f, "try-error")) next
    ll <- dr_try_num(stats::logLik(f))
    if (is.finite(ll) && ll > best_ll) {
      best_ll <- ll
      best    <- s
    }
  }
  if (is.null(best)) {
    stop("easydr: no ", shape$name, " candidate could be fitted to this ",
         "table. The Emax family needs doses above zero and enough levels ",
         "to estimate where the curve bends.", call. = FALSE)
  }
  best$profiled_over <- nrow(cand)
  # A profile that stops at the end of its own grid has not found a maximum,
  # it has run out of room. The value is still returned, because it is the
  # best of what was tried, but it is flagged so that nothing downstream
  # reports it as an estimate of where the curve bends.
  edge <- range(cand$ed50)
  if (isTRUE(all.equal(best$ed50, edge[1L])) ||
      isTRUE(all.equal(best$ed50, edge[2L]))) {
    best$ed50_at_boundary <- TRUE
    best$label <- paste0(best$label, " (at the edge of the search range)")
  }
  best
}
