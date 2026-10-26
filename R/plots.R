# The four plots, and the one thing all of them need: an estimate per study.
#
# A dose-response fit has no single effect per study the way a pairwise
# meta-analysis does. What it has, after the two-stage first stage, is one
# coefficient vector per study with its covariance. Every plot here is built
# from those: a forest of them, a funnel of them, and a bubble of them against
# a study characteristic. A one-stage fit never produces them, so these plots
# say so instead of inventing something to draw.

#' One estimate per study, on the scale the model was fitted on
#'
#' For a one-term shape the estimate is the study's own trend per unit of dose,
#' which is what the coefficient already is. For anything with more than one
#' term there is no single trend, so the estimate is the predicted effect at a
#' stated dose against the fit's reference dose, and the dose is reported
#' rather than assumed.
#'
#' The weight is the random-effects inverse-variance weight for the same
#' contrast: the within-study variance plus the between-study variance the fit
#' estimated, both projected onto the contrast. A common-effect fit has a zero
#' between-study term, which is the point of it.
#'
#' @keywords internal
#' @noRd
dr_study_effects <- function(x, at = NULL, what = "this plot") {
  if (!inherits(x, "easydr_result")) {
    stop("easydr: ", what, " takes a fitted result from easydr().",
         call. = FALSE)
  }
  bi <- x$fit$bi
  Si <- x$fit$Si
  if (is.null(bi) || is.null(Si) || !is.matrix(bi)) {
    stop("easydr: ", what, " needs one estimate per study, and only the ",
         "two-stage first stage produces them. This fit is one-stage, which ",
         "estimates every study inside a single model. Refit with ",
         "proc = \"2stage\" if you want it.", call. = FALSE)
  }

  frame   <- x$frame
  studies <- unique(as.character(frame$studlab))
  ids     <- vapply(studies, function(s)
    frame$id[as.character(frame$studlab) == s][1L], numeric(1))
  xref <- x$reference$curve
  unit <- if (is.null(x$dose_unit)) "dose unit" else x$dose_unit

  one_term <- x$shape$terms == 1L
  if (one_term && is.null(at)) {
    cvec    <- matrix(1, nrow = 1L, ncol = 1L)
    label   <- dr_label(x$sm, x$dose_unit)
    at_used <- NA_real_
  } else {
    if (is.null(at)) {
      # The first dose the fit already reports, so the plot and the printed
      # table are talking about the same thing.
      cand    <- x$predicted$dose[!x$predicted$is_ref]
      at_used <- if (length(cand) > 0L) cand[1L] else max(frame$dose)
    } else {
      if (!is.numeric(at) || length(at) != 1L || !is.finite(at)) {
        stop("easydr: `at` must be a single dose, such as at = 10.",
             call. = FALSE)
      }
      at_used <- at
    }
    cvec  <- dr_contrast(x$shape, at_used, xref)
    label <- paste0(x$sm, " at ", at_used, " ", unit, " against ", xref)
  }

  Psi <- x$Psi
  if (is.null(Psi)) Psi <- matrix(0, ncol(bi), ncol(bi))
  Psi  <- as.matrix(Psi)
  tau2 <- as.numeric(cvec %*% Psi %*% t(cvec))

  rows <- lapply(seq_along(studies), function(i) {
    id <- ids[i]
    if (is.na(id) || id > nrow(bi)) {
      return(data.frame(study = studies[i], levels = NA_integer_,
                        est = NA_real_, se = NA_real_, var_w = NA_real_,
                        stringsAsFactors = FALSE))
    }
    b   <- as.numeric(bi[id, ])
    V   <- as.matrix(Si[[id]])
    est <- as.numeric(cvec %*% b)
    vw  <- as.numeric(cvec %*% V %*% t(cvec))
    data.frame(
      study  = studies[i],
      levels = length(unique(frame$dose[as.character(frame$studlab) ==
                                          studies[i]])),
      est    = est,
      se     = sqrt(vw),
      var_w  = vw,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)

  zcrit <- stats::qnorm(1 - (1 - x$provenance$level) / 2)
  tab$lb <- tab$est - zcrit * tab$se
  tab$ub <- tab$est + zcrit * tab$se
  w <- 1 / (tab$var_w + max(tau2, 0))
  tab$weight <- 100 * w / sum(w, na.rm = TRUE)
  rownames(tab) <- NULL

  # The pooled value for the same contrast, so the diamond and the squares are
  # answering one question.
  pooled_est <- as.numeric(cvec %*% as.numeric(coef(x$fit)))
  pooled_var <- as.numeric(cvec %*% as.matrix(vcov(x$fit)) %*% t(cvec))
  pooled <- list(est = pooled_est, se = sqrt(pooled_var),
                 lb = pooled_est - zcrit * sqrt(pooled_var),
                 ub = pooled_est + zcrit * sqrt(pooled_var))

  list(table = tab, pooled = pooled, label = label, at = at_used,
       ratio = x$ratio, sm = x$sm, unit = unit, level = x$provenance$level,
       one_term = one_term && is.na(at_used))
}

#' Leave a margin round a range, multiplicatively on a ratio scale
#'
#' A log axis needs a multiplicative margin: adding a constant to a ratio makes
#' the gap at the bottom much larger than the gap at the top.
#'
#' @keywords internal
#' @noRd
dr_pad_range <- function(r, ratio, frac = 0.06) {
  if (!all(is.finite(r))) return(r)
  if (isTRUE(ratio)) {
    lr <- log(r)
    d  <- max(diff(lr), 1e-8) * frac
    return(exp(c(lr[1L] - d, lr[2L] + d)))
  }
  d <- max(diff(r), 1e-8) * frac
  c(r[1L] - d, r[2L] + d)
}

#' Spread the points of one column evenly, the same way every time
#'
#' Random jitter would move the picture between two runs of the same script,
#' which is the last thing a figure in a review should do.
#'
#' @keywords internal
#' @noRd
dr_spread <- function(group, width = 0.26) {
  out <- numeric(length(group))
  for (g in unique(group)) {
    i <- which(group == g)
    out[i] <- if (length(i) == 1L) 0 else seq(-width, width,
                                              length.out = length(i))
  }
  out
}

#' Axis ticks on the reported scale, for an axis drawn on the model scale
#' @keywords internal
#' @noRd
dr_plot_ticks <- function(lo, hi, ratio) {
  if (!ratio) {
    v <- pretty(c(lo, hi))
    return(list(at = v, labels = format(v)))
  }
  v <- pretty(c(exp(lo), exp(hi)))
  v <- v[v > 0]
  list(at = log(v), labels = format(v))
}

#' Forest plot of the study trends and the pooled one
#'
#' @description
#' One row per study, laid out the way a `meta` forest plot is: the study label
#' on the left, a square sized by its weight with its confidence interval, the
#' value with its interval on the right, and the pooled estimate as a diamond
#' underneath.
#'
#' @section What each row is:
#' A dose-response study contributes a trend, not an effect, so that is what a
#' row shows: the study's own slope per unit of dose, from the two-stage first
#' stage. For a shape with more than one term there is no single slope, so the
#' row becomes the predicted effect at one dose against the fit's reference
#' dose, and the heading names the dose. Pass `at` to choose it, for any shape.
#'
#' @section Why some fits cannot be drawn this way:
#' A one-stage fit estimates every study inside a single model and never
#' produces a per-study coefficient vector. There is nothing to put on the rows,
#' so this refuses rather than drawing the pooled line eight times.
#'
#' @param x a fitted result from [easydr()].
#' @param at a dose at which to compare studies, instead of their trends.
#'   Required in effect for curves, where no single trend exists; the default
#'   is the first dose the fit already reports.
#' @param main plot title. `NULL` writes one from the fit.
#' @param xlab label for the effect axis. `NULL` writes one from `sm`.
#' @param xlim limits for the effect axis, on the reported scale, such as
#'   `c(0.9, 1.5)` for a ratio.
#' @param order `"data"` keeps the order of the table, `"effect"` sorts by the
#'   estimate, `"weight"` sorts by weight.
#' @param digits digits for the printed estimates beside the rows.
#' @param col.square,col.diamond,col.null colours for the study squares, the
#'   pooled diamond and the null line.
#' @param ... further arguments passed to [graphics::plot.default()].
#'
#' @return The plotted table, invisibly: one row per study with `est`, `se`,
#'   `lb`, `ub` and `weight`, plus the `pooled` estimate. Returned so the
#'   numbers behind the picture can be checked.
#'
#' @seealso [easydr_curve()] for the fitted dose-response itself.
#'
#' @examples
#' fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
#'               dose.unit = "g/day")
#' easydr_forest(fit)
#'
#' # A curve has no single trend, so the rows become one dose.
#' rcs <- easydr(alcohol_cancer, sm = "IRR", shape = "rcs",
#'               dose.unit = "g/day")
#' easydr_forest(rcs, at = 25)
#'
#' @export
easydr_forest <- function(x, at = NULL, main = NULL, xlab = NULL, xlim = NULL,
                          order = c("data", "effect", "weight"), digits = 3,
                          col.square = "#1F2A44", col.diamond = "#C0392B",
                          col.null = "#8A939B", ...) {
  order <- match.arg(order)
  e   <- dr_study_effects(x, at, what = "a forest plot")
  tab <- e$table
  tab <- switch(order,
                data   = tab,
                effect = tab[order(tab$est), , drop = FALSE],
                weight = tab[order(-tab$weight), , drop = FALSE])
  k <- nrow(tab)

  null_v <- 0
  lo <- min(c(tab$lb, e$pooled$lb, null_v), na.rm = TRUE)
  hi <- max(c(tab$ub, e$pooled$ub, null_v), na.rm = TRUE)
  if (!is.null(xlim)) {
    v  <- if (e$ratio) log(xlim) else xlim
    lo <- min(v)
    hi <- max(v)
  }
  pad <- 0.04 * (hi - lo)
  lo  <- lo - pad
  hi  <- hi + pad
  mapx <- function(v) (v - lo) / (hi - lo)

  back <- if (e$ratio) exp else identity
  fmt  <- function(v) format(round(back(v), digits), nsmall = digits)
  txt  <- ifelse(is.na(tab$est), "not estimable",
                 paste0(fmt(tab$est), " (", fmt(tab$lb), " to ",
                        fmt(tab$ub), ")"))
  pooled_txt <- paste0(fmt(e$pooled$est), " (", fmt(e$pooled$lb), " to ",
                       fmt(e$pooled$ub), ")")

  if (is.null(main)) {
    main <- paste0("Dose-response trend per study (", x$sm_label, ")")
  }
  if (is.null(xlab)) xlab <- e$label

  op <- par(mar = c(4.5, 0.5, 3.5, 0.5))
  on.exit(par(op), add = TRUE)
  plot(NA, xlim = c(-0.95, 1.95), ylim = c(-2.1, k + 1.9),
       axes = FALSE, xlab = "", ylab = "", main = main, ...)

  ys <- rev(seq_len(k))
  text(-0.92, k + 1.1, "Study", adj = 0, font = 2, cex = 0.85)
  text(1.05, k + 1.1, paste0(x$sm, " (", format(100 * e$level), "% CI)"),
       adj = 0, font = 2, cex = 0.85)
  text(1.93, k + 1.1, "Weight", adj = 1, font = 2, cex = 0.85)
  segments(-0.95, k + 0.7, 1.95, k + 0.7, col = col.null)

  text(-0.92, ys, tab$study, adj = 0, cex = 0.85)
  text(1.05, ys, txt, adj = 0, cex = 0.85)
  text(1.93, ys, paste0(format(round(tab$weight, 1), nsmall = 1), "%"),
       adj = 1, cex = 0.85)

  segments(mapx(null_v), -1.5, mapx(null_v), k + 0.7, col = col.null,
           lty = 2)
  ok <- !is.na(tab$est)
  if (any(ok)) {
    segments(mapx(tab$lb[ok]), ys[ok], mapx(tab$ub[ok]), ys[ok],
             col = col.square)
    # Area proportional to weight, which is what a reader reads off a square.
    side <- 0.55 * sqrt(tab$weight / max(tab$weight, na.rm = TRUE))
    for (j in which(ok)) {
      polygon(mapx(tab$est[j]) + c(-1, 1, 1, -1) * side[j] * 0.045,
              ys[j] + c(-1, -1, 1, 1) * side[j] * 0.5,
              col = col.square, border = NA)
    }
  }

  segments(-0.95, -0.35, 1.95, -0.35, col = col.null)
  polygon(mapx(c(e$pooled$lb, e$pooled$est, e$pooled$ub, e$pooled$est)),
          -1.1 + c(0, 0.35, 0, -0.35), col = col.diamond, border = NA)
  text(-0.92, -1.1, "Pooled", adj = 0, font = 2, cex = 0.85)
  text(1.05, -1.1, pooled_txt, adj = 0, font = 2, cex = 0.85)
  text(1.93, -1.1, "100.0%", adj = 1, font = 2, cex = 0.85)

  ticks <- dr_plot_ticks(lo, hi, e$ratio)
  keep  <- ticks$at >= lo & ticks$at <= hi
  axis(1, at = mapx(ticks$at[keep]), labels = ticks$labels[keep],
       cex.axis = 0.85, line = -0.5)
  mtext(xlab, side = 1, line = 1.9, at = 0.5, cex = 0.9)

  invisible(list(table = tab, pooled = e$pooled, label = e$label, at = e$at))
}

#' Plot the fitted dose-response curve
#'
#' @description
#' The pooled dose-response against dose, with its confidence band, optionally
#' its prediction band, and the doses the fit reports marked and labelled. The
#' numbers were all computed by [easydr()]; this chooses what appears and how.
#'
#' @section What the two bands mean:
#' The confidence band describes where the pooled curve lies. The prediction
#' band describes where a new study's curve would fall, so it is wider, and a
#' common-effect fit has none at all because it assumes there is no
#' between-study variation for one to describe. Showing the prediction band and
#' calling it uncertainty about the pooled curve is the usual mistake.
#'
#' Doses beyond the observed range are drawn, because [easydr()] will predict
#' there if asked, but they are marked: the model has no data to support that
#' part of the curve.
#'
#' @param x a fitted result from [easydr()].
#' @param main plot title. `NULL` writes one from the fit, naming the shape.
#' @param xlab,ylab axis labels. `NULL` writes them from the dose unit and
#'   `sm`.
#' @param xlim,ylim axis limits on the reported scale. `NULL` fits them to
#'   whatever is drawn.
#' @param round digits for the labels beside the marked doses. `NULL` prints no
#'   labels.
#' @param prediction `TRUE` adds the prediction band when the fit has one.
#'   `FALSE`, the default, draws the confidence band alone.
#' @param points `TRUE`, the default, marks the doses the fit reports.
#' @param reference `TRUE`, the default, draws the null line and marks the
#'   reference dose.
#' @param col,col.ci,col.pi colours for the curve, the confidence band and the
#'   prediction band.
#' @param ... further arguments passed to [graphics::plot.default()].
#'
#' @return The plotted curve table, invisibly: `dose`, `est`, `ci.lb`,
#'   `ci.ub`, `pi.lb`, `pi.ub` and `extrapolated`.
#'
#' @seealso [easydr_forest()] for what each study contributed.
#'
#' @examples
#' fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
#'               dose.unit = "g/day", at = c(10, 25, 50))
#' easydr_curve(fit, main = "Alcohol and colorectal cancer", round = 2)
#'
#' # A spline, with the prediction band shown and the axis pinned.
#' rcs <- easydr(alcohol_cancer, sm = "IRR", shape = "rcs",
#'               dose.unit = "g/day", at = c(10, 25, 50))
#' easydr_curve(rcs, prediction = TRUE, ylim = c(0.8, 2))
#'
#' @export
easydr_curve <- function(x, main = NULL, xlab = NULL, ylab = NULL,
                         xlim = NULL, ylim = NULL, round = NULL,
                         prediction = FALSE, points = TRUE, reference = TRUE,
                         col = "#1F2A44", col.ci = "#C5C9CE",
                         col.pi = "#EDEFF2", ...) {
  if (!inherits(x, "easydr_result")) {
    stop("easydr: easydr_curve() takes a fitted result from easydr().",
         call. = FALSE)
  }
  tab  <- x$curve
  unit <- if (is.null(x$dose_unit)) "dose" else x$dose_unit
  show_pi <- isTRUE(prediction) && isTRUE(x$prediction$estimable) &&
    !all(is.na(tab$pi.lb))
  if (isTRUE(prediction) && !show_pi) {
    warning("easydr: this fit has no prediction interval (",
            x$prediction$reason, "), so only the confidence band is drawn.",
            call. = FALSE)
  }

  if (is.null(xlab)) xlab <- paste0("Dose (", unit, ")")
  if (is.null(ylab)) ylab <- x$sm_label
  if (is.null(main)) main <- paste0("Dose-response: ", x$shape$label)
  if (is.null(xlim)) xlim <- range(tab$dose)
  if (is.null(ylim)) {
    v    <- c(tab$ci.lb, tab$ci.ub, if (show_pi) c(tab$pi.lb, tab$pi.ub))
    ylim <- range(c(v, if (x$ratio) 1 else 0), na.rm = TRUE)
    ylim <- dr_pad_range(ylim, x$ratio)
  }

  op <- par(mar = c(4.8, 4.5, 3.8, 1))
  on.exit(par(op), add = TRUE)
  plot(NA, xlim = xlim, ylim = ylim, log = if (x$ratio) "y" else "",
       xlab = xlab, ylab = ylab, main = main, las = 1, ...)

  if (show_pi) {
    polygon(c(tab$dose, rev(tab$dose)), c(tab$pi.lb, rev(tab$pi.ub)),
            col = col.pi, border = NA)
  }
  polygon(c(tab$dose, rev(tab$dose)), c(tab$ci.lb, rev(tab$ci.ub)),
          col = col.ci, border = NA)
  if (isTRUE(reference)) {
    abline(h = if (x$ratio) 1 else 0, col = "#8A939B", lty = 2)
    abline(v = x$reference$curve, col = "#8A939B", lty = 3)
  }
  lines(tab$dose, tab$est, col = col, lwd = 2.2)

  if (any(tab$extrapolated)) {
    ex <- tab[tab$extrapolated, , drop = FALSE]
    points(ex$dose, ex$est, pch = 4, col = "#C0392B", cex = 0.7)
    mtext("x marks doses beyond the observed range", side = 1, line = 3.4,
          adj = 1, cex = 0.75, col = "#C0392B")
  }

  marks <- x$predicted[!x$predicted$is_ref, , drop = FALSE]
  if (isTRUE(points) && nrow(marks) > 0L) {
    points(marks$dose, marks$est, pch = 21, bg = "white", col = col,
           cex = 1.1, lwd = 1.6)
    if (!is.null(round)) {
      text(marks$dose, marks$est,
           format(base::round(marks$est, round), nsmall = round),
           pos = 3, offset = 0.6, cex = 0.8, col = col)
    }
  }

  bands <- paste0(format(100 * x$provenance$level), "% CI")
  if (show_pi) {
    bands <- paste0(bands, " and ", format(100 * x$provenance$level), "% PI")
  }
  mtext(paste0("Reference ", x$reference$curve, " ", unit, ".  Band: ", bands),
        side = 3, line = 0.3, adj = 0, cex = 0.78, col = "#5B6475")

  invisible(tab)
}

#' Egger's regression test for small-study effects
#'
#' Regresses each study's standard normal deviate on its precision. A non-zero
#' intercept is the asymmetry: small studies sitting away from the pooled value.
#'
#' The test is for the trend estimates, one per study, which the two-stage
#' first stage produces. It is underpowered below about ten studies, and that
#' is reported alongside rather than left to be remembered.
#'
#' @keywords internal
#' @noRd
dr_egger <- function(est, se, k_warn = 10L) {
  ok  <- is.finite(est) & is.finite(se) & se > 0
  est <- est[ok]
  se  <- se[ok]
  k   <- length(est)
  if (k < 3L) {
    return(list(estimable = FALSE, k = k,
                reason = paste0("Egger's test needs at least 3 studies with a ",
                                "standard error, and ", k, " ",
                                if (k == 1L) "is" else "are", " available")))
  }
  # Every study the same size means there is no small-study contrast, and the
  # precision column is constant: the regression is rank deficient and has no
  # intercept to read. That is a real answer about the data, not a failure.
  if (length(unique(se)) < 2L) {
    return(list(estimable = FALSE, k = k,
                reason = paste0("every study has the same standard error, so ",
                                "there is no contrast between large and small ",
                                "studies for the regression to pick up")))
  }
  m  <- stats::lm(I(est / se) ~ I(1 / se))
  cf <- summary(m)$coefficients
  if (nrow(cf) < 2L || anyNA(cf[, 1L])) {
    return(list(estimable = FALSE, k = k,
                reason = "the precisions do not vary enough to fit the line"))
  }
  list(estimable = TRUE, k = k,
       intercept = cf[1L, 1L], se = cf[1L, 2L],
       t = cf[1L, 3L], df = k - 2L, p = cf[1L, 4L],
       slope = cf[2L, 1L],
       underpowered = k < k_warn, k_warn = k_warn)
}

#' Funnel plot with a test for small-study effects
#'
#' @description
#' Each study's dose-response trend against its standard error, with the pooled
#' value and the pseudo confidence contours, plus Egger's regression test for
#' asymmetry.
#'
#' @section Why Egger and not Deeks:
#' Deeks' funnel plot asymmetry test is for diagnostic test accuracy, where the
#' effective sample size rather than the standard error is the sensible
#' horizontal axis. It does not apply here. A dose-response review has one trend
#' estimate per study with a standard error, which is exactly the setting
#' Egger's regression was written for, so that is the test: regress each study's
#' standard normal deviate on its precision and read the intercept.
#'
#' @section What asymmetry is and is not:
#' A funnel plot asks whether small studies sit away from the pooled value. It
#' cannot say why. Publication bias produces that pattern, and so do real
#' differences between small and large studies, poorer methods in smaller ones,
#' and chance. Below about ten studies the test has little power and a null
#' result means almost nothing, which the output says.
#'
#' A curve has no single trend, so for a shape with more than one term the
#' comparison is the predicted effect at one dose, named in the output.
#'
#' @param x a fitted result from [easydr()].
#' @param at a dose at which to compare studies. See [easydr_forest()].
#' @param main plot title. `NULL` writes one from the fit.
#' @param xlab,ylab axis labels. `NULL` writes them from `sm`.
#' @param contours `TRUE`, the default, draws the pseudo confidence contours
#'   around the pooled value.
#' @param labels `TRUE` writes the study labels beside the points.
#' @param digits digits for the test reported under the plot.
#' @param ... further arguments passed to [graphics::plot.default()].
#'
#' @return Invisibly, a list with the plotted `table` and `egger`: the
#'   intercept with its standard error, `t`, `df` and `p`, the number of
#'   studies, and whether that number makes the test underpowered.
#'
#' @references
#' Egger M, Davey Smith G, Schneider M, Minder C (1997). Bias in meta-analysis
#' detected by a simple, graphical test. BMJ 315, 629-634.
#'
#' Sterne JAC, Sutton AJ, Ioannidis JPA, et al. (2011). Recommendations for
#' examining and interpreting funnel plot asymmetry in meta-analyses of
#' randomised controlled trials. BMJ 343, d4002.
#'
#' @examples
#' fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
#'               dose.unit = "g/day")
#' f <- easydr_funnel(fit)
#' f$egger
#'
#' @export
easydr_funnel <- function(x, at = NULL, main = NULL, xlab = NULL, ylab = NULL,
                          contours = TRUE, labels = FALSE, digits = 3, ...) {
  e   <- dr_study_effects(x, at, what = "a funnel plot")
  tab <- e$table
  ok  <- !is.na(tab$est) & !is.na(tab$se) & tab$se > 0
  if (!any(ok)) {
    stop("easydr: no study has an estimate with a standard error, so there is ",
         "nothing to plot.", call. = FALSE)
  }
  eg <- dr_egger(tab$est, tab$se)

  back <- if (e$ratio) exp else identity
  if (is.null(main)) main <- "Funnel plot of the study trends"
  if (is.null(xlab)) xlab <- e$label
  if (is.null(ylab)) ylab <- "Standard error of the trend"

  se_max <- max(tab$se[ok]) * 1.12
  halfw  <- 1.96 * se_max
  xr <- range(c(tab$est[ok], e$pooled$est - halfw, e$pooled$est + halfw))

  op <- par(mar = c(7.2, 4.5, 3.5, 1))
  on.exit(par(op), add = TRUE)
  plot(NA, xlim = if (e$ratio) back(xr) else xr, ylim = c(se_max, 0),
       log = if (e$ratio) "x" else "", xlab = xlab, ylab = ylab,
       main = main, las = 1, ...)

  if (isTRUE(contours)) {
    for (z in c(1.96, 2.58)) {
      lines(back(c(e$pooled$est, e$pooled$est - z * se_max)), c(0, se_max),
            col = "#C5C9CE", lty = 2)
      lines(back(c(e$pooled$est, e$pooled$est + z * se_max)), c(0, se_max),
            col = "#C5C9CE", lty = 2)
    }
  }
  abline(v = back(e$pooled$est), col = "#C0392B")
  points(back(tab$est[ok]), tab$se[ok], pch = 21, bg = "#1F2A44",
         col = "white", cex = 1.3)
  if (isTRUE(labels)) {
    u    <- par("usr")
    mid  <- if (e$ratio) 10^mean(u[1:2]) else mean(u[1:2])
    side <- ifelse(back(tab$est[ok]) > mid, 2L, 4L)
    text(back(tab$est[ok]), tab$se[ok], tab$study[ok], pos = side, cex = 0.7,
         col = "#5B6475")
  }

  note <- if (isTRUE(eg$estimable)) {
    paste0("Egger's test: intercept ",
           format(round(eg$intercept, digits), nsmall = digits),
           " (SE ", format(round(eg$se, digits), nsmall = digits),
           "), t = ", format(round(eg$t, 2), nsmall = 2),
           " on ", eg$df, " df, p = ", format.pval(eg$p, digits = 3),
           ", ", eg$k, " studies")
  } else {
    paste0("Egger's test: not estimable (", eg$reason, ")")
  }
  mtext(note, side = 1, line = 4.1, adj = 0, cex = 0.78)
  if (isTRUE(eg$estimable) && isTRUE(eg$underpowered)) {
    mtext(paste0("With fewer than ", eg$k_warn, " studies this test has ",
                 "little power; a null result means very little."),
          side = 1, line = 5.1, adj = 0, cex = 0.78, col = "#C0392B")
  }
  mtext("Asymmetry has several causes. Publication bias is only one.",
        side = 1, line = 6.1, adj = 0, cex = 0.75, col = "#5B6475")

  invisible(list(table = tab, pooled = e$pooled, egger = eg, label = e$label,
                 at = e$at))
}

#' Bubble plot of a dose-response meta-regression
#'
#' @description
#' Each study's own trend against the moderator, with the bubble area
#' proportional to the weight the study carries, and the line the
#' meta-regression fitted. A categorical moderator gets one column of bubbles
#' per level with the fitted level means drawn across them.
#'
#' @section What the line is:
#' The line comes from the meta-regression in `x`, not from an ordinary
#' regression through the bubbles. The bubbles are the observed study trends
#' from the unadjusted fit, which is why a study can sit off the line: that
#' distance is the residual heterogeneity the printed output quantifies.
#'
#' @section Why only one-term shapes:
#' A bubble plot needs one number per study on the vertical axis. A curve has no
#' single trend, so there is nothing to put there, and projecting it onto one
#' dose would quietly answer a different question from the regression the line
#' comes from. For curves, read [easydr_metareg()]'s table instead.
#'
#' @param x an object from [easydr_metareg()].
#' @param mod which moderator to draw, when the regression has more than one.
#'   `NULL` uses the first.
#' @param main plot title. `NULL` writes one naming the moderator.
#' @param xlab,ylab axis labels. `NULL` writes them from the moderator and
#'   `sm`.
#' @param ylim axis limits on the reported scale.
#' @param labels `TRUE` writes the study labels beside the bubbles.
#' @param col,col.line colours for the bubbles and the fitted line.
#' @param ... further arguments passed to [graphics::plot.default()].
#'
#' @return Invisibly, a list with the plotted `table` (study, moderator value,
#'   trend, standard error, weight) and the `fitted` line or level means.
#'
#' @seealso [easydr_metareg()], which fits it and reports the test.
#'
#' @examples
#' fit <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
#' mr <- easydr_metareg(fit, ~ design)
#' easydr_bubble(mr)
#'
#' @export
easydr_bubble <- function(x, mod = NULL, main = NULL, xlab = NULL,
                          ylab = NULL, ylim = NULL, labels = FALSE,
                          col = "#4E7CA8", col.line = "#C0392B", ...) {
  if (!inherits(x, "easydr_moderator")) {
    stop("easydr: easydr_bubble() takes the result of easydr_metareg().",
         call. = FALSE)
  }
  inh <- x$inherited
  if (inh$shape$terms != 1L) {
    stop("easydr: a bubble plot needs one trend per study on the vertical ",
         "axis, and a ", inh$shape$label, " has none: no single number ",
         "describes it. Read the table from easydr_metareg() instead, or ",
         "refit with shape = \"linear\" if a trend is what you want to ",
         "regress.", call. = FALSE)
  }
  if (is.null(mod)) mod <- x$moderators[1L]
  if (!mod %in% x$moderators) {
    stop("easydr: ", dQuote(mod), " is not a moderator in this regression. ",
         "It has: ", paste(x$moderators, collapse = ", "), ".", call. = FALSE)
  }

  # The bubbles are the unadjusted study trends, refitted through the settings
  # the regression inherited, so nothing about the model is re-derived here.
  base_fit <- dr_refit(inh, inh$studies)
  e   <- dr_study_effects(base_fit, NULL, what = "a bubble plot")
  tab <- e$table

  st      <- x$study_table
  tab$mod <- st[[mod]][match(tab$study, st$.study)]
  tab <- tab[!is.na(tab$est) & !is.na(tab$mod), , drop = FALSE]
  if (nrow(tab) < 2L) {
    stop("easydr: fewer than two studies have both a trend and a moderator ",
         "value, so there is nothing to draw.", call. = FALSE)
  }

  back       <- if (x$ratio) exp else identity
  continuous <- identical(unname(x$types[[mod]]), "continuous")
  if (is.null(ylab)) ylab <- dr_label(x$sm, x$dose_unit)
  if (is.null(xlab)) xlab <- mod
  if (is.null(main)) main <- paste0("Dose-response trend against ", mod)
  # The bubbles carry no interval, so the axis is scaled to the estimates and
  # the fitted line. Scaling it to the confidence limits would leave the
  # bubbles in a band at the middle of an empty plot.
  if (is.null(ylim)) {
    ylim <- dr_pad_range(range(back(tab$est), na.rm = TRUE), x$ratio, 0.18)
  }

  # The fitted line is the regression's, evaluated on the moderator. For a
  # one-term shape the coefficients are the trend at the reference level plus
  # the change in trend, so the model matrix gives the fitted trend directly.
  b <- coef(x$fit)
  predict_at <- function(values) {
    nd <- st[rep(1L, length(values)), , drop = FALSE]
    nd[[mod]] <- values
    mm <- stats::model.matrix(stats::reformulate(x$moderators), data = nd)
    as.numeric(mm %*% b[seq_len(ncol(mm))])
  }

  op <- par(mar = c(4.8, 4.8, 3.8, 1))
  on.exit(par(op), add = TRUE)

  if (continuous) {
    xr <- range(tab$mod)
    plot(NA, xlim = xr + c(-0.04, 0.04) * diff(xr), ylim = ylim,
         log = if (x$ratio) "y" else "", xlab = xlab, ylab = ylab,
         main = main, las = 1, ...)
    abline(h = if (x$ratio) 1 else 0, col = "#8A939B", lty = 2)
    grid_v <- seq(xr[1L], xr[2L], length.out = 60)
    fitted <- data.frame(mod = grid_v, trend = predict_at(grid_v))
    lines(fitted$mod, back(fitted$trend), col = col.line, lwd = 2)
    px <- tab$mod
  } else {
    lv <- unique(as.character(tab$mod))
    px <- match(as.character(tab$mod), lv)
    # Studies in one level would otherwise be drawn on the same x. Spread them
    # evenly across the column instead of jittering at random, so the picture
    # is the same every time it is drawn.
    px <- px + dr_spread(px)
    plot(NA, xlim = c(0.6, length(lv) + 0.4), ylim = ylim,
         log = if (x$ratio) "y" else "", xlab = xlab, ylab = ylab,
         main = main, las = 1, axes = FALSE, ...)
    axis(1, at = seq_along(lv), labels = lv, cex.axis = 0.9)
    axis(2, las = 1)
    box()
    abline(h = if (x$ratio) 1 else 0, col = "#8A939B", lty = 2)
    fitted <- data.frame(mod = lv, trend = predict_at(lv),
                         stringsAsFactors = FALSE)
    segments(seq_along(lv) - 0.3, back(fitted$trend),
             seq_along(lv) + 0.3, back(fitted$trend),
             col = col.line, lwd = 2.4)
  }

  cex <- 2.4 * sqrt(tab$weight / max(tab$weight))
  points(px, back(tab$est), pch = 21, bg = adjustcolor(col, 0.55),
         col = col, cex = cex, lwd = 1.4)
  if (isTRUE(labels)) {
    # A label on the right of a point near the right edge runs off the plot,
    # so points past the midpoint are labelled on their left instead.
    mid  <- mean(par("usr")[1:2])
    side <- ifelse(px > mid, 2L, 4L)
    text(px, back(tab$est), tab$study, pos = side, cex = 0.7,
         col = "#5B6475", offset = 0.3 + 0.08 * cex)
  }

  note <- if (isTRUE(x$test$estimable)) {
    paste0("Between-moderator test: X2 = ", round(x$test$chi2, 2),
           " (df = ", x$test$df, "), p = ", format.pval(x$test$p, digits = 3))
  } else {
    "Between-moderator test: not estimable"
  }
  mtext(paste0(note, ".  Bubble area is the weight the study carries."),
        side = 3, line = 0.3, adj = 0, cex = 0.78, col = "#5B6475")

  invisible(list(table = tab, fitted = fitted, moderator = mod,
                 test = x$test))
}
