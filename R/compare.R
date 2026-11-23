# Two comparison workflows, after Crippa's own scripts.
#
# fp_splines_example.txt fits the candidate shapes on one table, ranks them by
# AIC, and overlays the fitted curves on a single panel.
#
# knots_location.txt takes the 10th, 25th, 50th, 75th and 90th percentiles of
# the dose distribution, forms every three-knot combination of them, refits a
# restricted cubic spline for each, and overlays those.
#
# Both answer the same worry from two directions: how much of the curve is the
# data, and how much is the choice of model.

#' Refit a result under a different shape, keeping everything else
#' @keywords internal
#' @noRd
dr_refit_shape <- function(f, shape, knots = NULL) {
  # The powers and the fixed Emax parameters belong to the shape that was
  # fitted, not to the one being tried.
  dr_refit_with(f, list(shape = shape, knots = knots, powers = NULL,
                        ed50 = NULL, steepness = NULL))
}

#' Colours and line types for overlaid curves, keyed by label
#'
#' The panel and the table under it draw the same curves, so both take their
#' colour and line type from this one lookup, by name. The table is sorted by
#' AIC and the panel by the order the shapes were given, and matching by
#' position would give a row another curve's line.
#'
#' @keywords internal
#' @noRd
dr_overlay_style <- function(labels) {
  pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
           "#56B4E9", "#444444", "#8C564B", "#7F3C8D", "#B8860B")
  lty <- c("solid", "dashed", "dotdash", "longdash", "twodash", "dotted")
  n <- length(labels)
  list(colour   = stats::setNames(rep_len(pal, n), labels),
       linetype = stats::setNames(rep_len(lty, n), labels))
}

#' Overlay several fitted curves on one panel
#'
#' The panels of `mult.curves` answer "what does each shape look like". This
#' answers "how far apart are they", which needs them on the same axes and in
#' the same panel, and is the right panel of Crippa and Orsini's figure 3.
#' The table under the panel is the legend, so the panel draws none unless
#' asked to (`legend = TRUE`, for a figure drawn without its table).
#'
#' @keywords internal
#' @noRd
dr_overlay_panel <- function(fits, labels, ratio, sm_label, unit,
                             xlab = NULL, ylab = NULL,
                             legend.title = "Model", style = NULL,
                             legend = FALSE) {
  rows <- do.call(rbind, lapply(seq_along(fits), function(i) {
    d <- fits[[i]]$curve
    data.frame(dose = d$dose, est = d$est,
               model = factor(labels[i], levels = labels),
               stringsAsFactors = FALSE)
  }))
  if (is.null(xlab)) xlab <- paste0("Dose (", unit, ")")
  if (is.null(ylab)) ylab <- sm_label
  if (is.null(style)) style <- dr_overlay_style(labels)

  p <- ggplot2::ggplot(rows, ggplot2::aes(x = .data$dose, y = .data$est,
                                          colour = .data$model,
                                          linetype = .data$model)) +
    ggplot2::geom_hline(yintercept = if (isTRUE(ratio)) 1 else 0,
                        linetype = "dashed", colour = "grey55",
                        linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::scale_colour_manual(values = style$colour, breaks = labels,
                                 name = legend.title) +
    ggplot2::scale_linetype_manual(values = style$linetype, breaks = labels,
                                   name = legend.title) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey92"),
                   axis.text = ggplot2::element_text(colour = "black"),
                   legend.position = if (isTRUE(legend)) "right" else "none",
                   legend.key.width = grid::unit(26, "pt"),
                   legend.text = ggplot2::element_text(size = 8),
                   legend.title = ggplot2::element_text(size = 8),
                   plot.margin = ggplot2::margin(6, 10, 4, 6))
  if (isTRUE(ratio)) p <- p + ggplot2::scale_y_log10()
  p
}

#' The table under an overlay, which is also its legend
#'
#' Each row carries the line of its curve beside the name. The note goes
#' under the table, inside the table's own columns, so it lines up with the
#' names and nothing sits between the panel and the table but the table's
#' thin top rule.
#'
#' @keywords internal
#' @noRd
dr_overlay_table <- function(df, style, best = integer(0), note = NULL,
                             base_size = 12) {
  nm <- as.character(df[[1L]])
  tg <- dr_stats_grob(df, base_size = base_size,
                      key = list(colour = unname(style$colour[nm]),
                                 linetype = unname(style$linetype[nm])),
                      bold_rows = best)
  if (is.null(note)) return(tg)
  # Wrapped to the table's own width, so a long note never runs past it.
  note <- strwrap(paste(note, collapse = " "),
                  width = max(40L, floor(dr_grob_width_in(tg) / 0.07)))
  tg <- gtable::gtable_add_rows(tg, grid::unit(13 * length(note) + 10,
                                               "points"))
  gtable::gtable_add_grob(
    tg, grid::textGrob(paste(note, collapse = "\n"), x = grid::unit(5, "pt"),
                       y = grid::unit(1, "npc") - grid::unit(6, "pt"),
                       hjust = 0, vjust = 1,
                       gp = grid::gpar(fontsize = 9, col = "grey30",
                                       lineheight = 1.2)),
    t = nrow(tg), l = 1, r = ncol(tg), name = "note", clip = "off")
}

#' Draw an overlay with its table underneath, and save it if asked
#' @keywords internal
#' @noRd
dr_overlay_draw <- function(p, tg, main, foot, spec) {
  # A width asked for is kept unless the table under the panel needs more:
  # a table cut at both edges shows nothing.
  if (!is.null(spec) && !is.null(spec$width) && !is.null(tg)) {
    spec$width <- max(spec$width, dr_grob_width_in(tg) + 0.4)
  }
  title_in <- if (is.null(main)) 0 else 1.9 * 12 / 72
  w_in <- max(7.6, dr_grob_width_in(tg) + 0.5)
  panel_in <- min(max(0.46 * w_in, 3.2), 4.6)
  if (is.null(tg)) {
    # No table: the note is the footer, wrapped to the figure's width and
    # measured the way every other figure's footer is.
    foot <- strwrap(paste(foot, collapse = " "),
                    width = max(40L, floor(w_in / 0.068)))
    size <- dr_fig_size(w_in, panel_in + title_in +
                          dr_footer_in(foot, NULL) + 0.1)
    con <- dr_begin_file(spec, size[["width"]], size[["height"]])
    on.exit(dr_end_file(con), add = TRUE)
    return(invisible(dr_draw_stack(p, NULL, title = main, footer = foot,
                                   title_hjust = 0.5)))
  }
  # The table sits flush under the panel at its own height. Padding here
  # would open a band of white between the two; the few points kept are
  # what holds the note off the bottom edge of a saved file.
  th <- dr_grob_height(tg, pad = 6)
  size <- dr_fig_size(w_in, panel_in + title_in + dr_height_in(th))
  con <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)
  invisible(dr_draw_stack(p, NULL, title = main, footer = tg, footer_h = th,
                          title_hjust = 0.5))
}

#' Fit the candidate shapes and compare them
#'
#' Fits each shape in `shapes` to the table the result came from, ranks them
#' by AIC, and overlays the fitted curves on one panel with the comparison
#' table underneath. This is the workflow of Crippa's `fp_splines_example`
#' script, at
#' <https://alecri.github.io/downloads/codes/fp_splines_example.txt>.
#'
#' The shapes being compared are not nested inside one another: a spline is
#' not a quadratic with a term added, and an Emax curve is neither. No Wald
#' test holds between them, which is why the comparison is by information
#' criterion. Both are computed by maximum likelihood, because a REML
#' likelihood depends on the design matrix and two shapes would not be
#' comparable, and both charge for the parameters an Emax profile searched
#' over.
#'
#' @param x a fitted result from [easydr()].
#' @param shapes the shapes to fit and compare.
#' @param plot `TRUE` draws the overlay.
#' @param main a title, or `FALSE` for none.
#' @param digits digits for the reported effects.
#' @param result `TRUE` puts the comparison table under the figure. `FALSE`
#'   leaves the panel on its own, which is what a slide wants.
#' @param wald,deviance,r2,r2adj,aic,delta.aic,bic `TRUE` (`T`) or `FALSE`
#'   (`F`): show or hide that column in the table under the figure. `wald` is
#'   the p value for non-linearity, `deviance` the deviance test p value,
#'   `r2` the R-squared, `r2adj` the adjusted R-squared, `aic` the AIC itself,
#'   `delta.aic` its distance from the best shape, and `bic` the BIC. By
#'   default the figure shows the Wald test, the deviance p value, the
#'   adjusted R-squared and the Delta AIC.
#' @param deviance.value `TRUE` adds the deviance itself beside its p value.
#'   `FALSE`, the default, shows the p value only.
#'   The shapes are ranked by AIC whichever are shown, and the returned
#'   `table` always holds every column.
#' @param filetype,path,file.name,width,height where to save the figure, as
#'   in [easydr_curve()].
#'
#' @return An object of class `easydr_comparison`, holding `table`, ordered by
#'   AIC, and `fits`, the fitted models by shape.
#'
#' @examples
#' \donttest{
#' data(aripiprazole)
#' easydr.object <- easydr(aripiprazole, sm = "MD", shape = "linear",
#'                         dose.unit = "mg/day")
#' easydr_compare_shapes(easydr.object)
#' }
#'
#' @export
easydr_compare_shapes <- function(x,
                                  shapes = c("linear", "quadratic", "rcs",
                                             "fp", "hyperbolic_emax"),
                                  plot = TRUE, main = NULL, digits = 3,
                                  result = TRUE,
                                  wald = TRUE, deviance = TRUE,
                                  deviance.value = FALSE, r2 = FALSE,
                                  r2adj = TRUE, aic = FALSE,
                                  delta.aic = TRUE, bic = FALSE,
                                  filetype = c("png", "pdf", "tiff"),
                                  path = NULL, file.name = NULL,
                                  width = NULL, height = NULL) {
  dr_require_single(x, "easydr_compare_shapes()")
  if (!is.character(shapes) || length(shapes) < 2L) {
    stop("easydr: `shapes` names at least two shapes to compare.",
         call. = FALSE)
  }
  fits <- lapply(shapes, function(s) {
    if (identical(s, x$shape$name)) x else dr_refit_shape(x, s)
  })
  names(fits) <- shapes
  failed <- shapes[vapply(fits, is.null, logical(1))]
  fits   <- fits[!vapply(fits, is.null, logical(1))]
  if (length(fits) < 2L) {
    stop("easydr: fewer than two of those shapes could be fitted to this ",
         "table, so there is nothing to compare.", call. = FALSE)
  }

  m   <- dr_multi(fits, names(fits), match.call())
  tab <- dr_multi_table(m, digits)
  tab <- tab[order(tab$AIC, na.last = TRUE), , drop = FALSE]
  rownames(tab) <- NULL

  if (isTRUE(plot)) {
    unit <- if (is.null(x$dose_unit)) "dose" else x$dose_unit
    if (is.null(main)) main <- "Candidate shapes"
    if (identical(main, FALSE)) main <- NULL
    labs  <- dr_shape_title(names(fits))
    style <- dr_overlay_style(labs)
    p  <- dr_overlay_panel(fits, labs, x$ratio, x$sm_label, unit,
                           ylab = dr_sm_axis(x$sm),
                           legend.title = "Shape", style = style,
                           legend = !isTRUE(result))
    # The table is the point of this figure, so it is set larger than the
    # summary tables under a single curve. It is also the legend, and its
    # first row, the lowest AIC, is set in bold.
    shown <- dr_test_columns(tab, wald, deviance, r2, r2adj, aic, delta.aic,
                             bic, deviance.value)
    shown$Shape <- dr_shape_title(shown$Shape)
    foot <- paste0("Ranked by AIC. Reference ", x$reference$curve, " ",
                   unit, ".")
    best <- if (is.finite(tab$AIC[1L])) 1L else integer(0)
    tg <- if (isTRUE(result)) {
      dr_overlay_table(shown, style, best = best, note = foot)
    } else NULL
    spec <- dr_file_spec(file.name, filetype, path, width, height,
                         what = "easydr_compare_shapes()")
    dr_overlay_draw(p, tg, main, foot, spec)
  }

  structure(list(table = tab, fits = fits, failed = failed,
                 sm_label = x$sm_label, unit = x$dose_unit),
            class = "easydr_comparison")
}

#' Print a shape comparison
#'
#' @param x an object from [easydr_compare_shapes()].
#' @param ... ignored.
#' @return `x`, invisibly. Called for the printed output.
#' @export
print.easydr_comparison <- function(x, ...) {
  cat("Candidate shapes, ranked by AIC (", x$sm_label, ")\n\n", sep = "")
  tab <- x$table
  tab$Shape <- dr_shape_title(tab$Shape)
  print(tab, row.names = FALSE, right = FALSE)
  if (length(x$failed)) {
    cat("\nCould not be fitted to this table: ",
        paste(x$failed, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Refit a spline across knot locations
#'
#' Takes the percentiles of the dose distribution, forms every combination of
#' `nknots` of them, refits a restricted cubic spline at each, and overlays
#' the curves. This is the workflow of Crippa's `knots_location` script, at
#' <https://alecri.github.io/downloads/codes/knots_location.txt>, whose
#' default of the 10th, 25th, 50th, 75th and 90th percentiles taken three at a
#' time gives ten curves.
#'
#' Knot placement is a decision the protocol should carry, and this shows what
#' the decision was worth. Curves that lie on top of one another say the
#' result is the data. Curves that fan apart say the reader is looking at the
#' knots.
#'
#' @param x a fitted result from [easydr()].
#' @param percentiles the percentiles of the dose distribution to draw knots
#'   from.
#' @param nknots how many knots in each combination.
#' @param sets a list of knot vectors, in dose units, to use instead of the
#'   percentile combinations, such as the sets a published paper drew. When
#'   given, `percentiles` and `nknots` are ignored.
#' @param plot `TRUE` draws the overlay.
#' @param main a title, or `FALSE` for none.
#' @param result `TRUE` puts the table of knot sets under the figure. `FALSE`
#'   leaves the curves on their own, which is what a slide wants.
#' @param filetype,path,file.name,width,height where to save the figure, as
#'   in [easydr_curve()].
#'
#' @return An object of class `easydr_knots`, holding `table`, one row per
#'   knot set, and `fits`.
#'
#' @examples
#' \donttest{
#' data(alcohol_cvd)
#' easydr.object <- easydr(alcohol_cvd, sm = "RR", shape = "rcs",
#'                         dose.unit = "g/day")
#' easydr_knot_sensitivity(easydr.object)
#' }
#'
#' @export
easydr_knot_sensitivity <- function(x,
                                    percentiles = c(0.1, 0.25, 0.5, 0.75,
                                                    0.9),
                                    nknots = 3L, sets = NULL, plot = TRUE,
                                    main = NULL, result = TRUE,
                                    filetype = c("png", "pdf", "tiff"),
                                    path = NULL, file.name = NULL,
                                    width = NULL, height = NULL) {
  dr_require_single(x, "easydr_knot_sensitivity()")
  given <- !is.null(sets)
  if (given) {
    if (!is.list(sets) || length(sets) < 2L ||
        !all(vapply(sets, function(k) is.numeric(k) && length(k) >= 3L,
                    logical(1)))) {
      stop("easydr: `sets` is a list of at least two knot vectors, each ",
           "with three or more doses, such as ",
           "list(c(0, 10, 30), c(0.5, 10, 30)).", call. = FALSE)
    }
    sets <- lapply(sets, function(k) sort(as.numeric(k)))
  } else {
  if (!is.numeric(percentiles) || any(percentiles <= 0) ||
      any(percentiles >= 1) || length(percentiles) < nknots) {
    stop("easydr: `percentiles` are positions in the dose distribution ",
         "between 0 and 1, and there must be at least `nknots` of them.",
         call. = FALSE)
  }
  pts <- round(stats::quantile(x$frame$dose, percentiles, names = FALSE), 2)
  pts <- sort(unique(pts))
  if (length(pts) < nknots) {
    stop("easydr: those percentiles give only ", length(pts), " distinct ",
         "doses, which is fewer than the ", nknots, " knots asked for. ",
         "A short dose range cannot carry this sensitivity analysis.",
         call. = FALSE)
  }
  sets <- utils::combn(pts, nknots, simplify = FALSE)
  if (length(sets) < 2L) {
    stop("easydr: those percentiles give only one set of ", nknots,
         " knots, so there is nothing to compare it against. Name more ",
         "percentiles, such as the default c(0.1, 0.25, 0.5, 0.75, 0.9), ",
         "which gives ten sets of three.", call. = FALSE)
  }
  }
  labs <- vapply(sets, function(k) paste(k, collapse = ", "), character(1))

  fits <- lapply(sets, function(k) dr_refit_shape(x, "rcs", knots = k))
  ok   <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2L) {
    stop("easydr: only ", sum(ok), " of the ", length(sets), " knot sets ",
         "could be fitted, so there is no sensitivity to show. A spline ",
         "needs more distinct doses than it has knots.", call. = FALSE)
  }
  fits <- fits[ok]
  labs <- labs[ok]
  names(fits) <- labs

  at  <- dr_result_at(x)
  tab <- do.call(rbind, lapply(seq_along(fits), function(i) {
    f   <- fits[[i]]
    col <- dr_result_column(f, at, 3)
    data.frame(Knots = labs[i], Effect = col[1L], CI = col[2L],
               Deviance = if (is.null(f$gof$deviance)) NA_real_ else
                 round(f$gof$deviance, 2),
               AIC = round(dr_ic(f, "AIC", dr_refit_ml(f)), 1),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  names(tab)[2L] <- dr_pooled_at(fits[[1L]], at)$label
  rownames(tab)  <- NULL

  if (isTRUE(plot)) {
    unit <- if (is.null(x$dose_unit)) "dose" else x$dose_unit
    if (is.null(main)) main <- "The same spline, at every knot placement"
    if (identical(main, FALSE)) main <- NULL
    style <- dr_overlay_style(labs)
    p  <- dr_overlay_panel(fits, labs, x$ratio, x$sm_label, unit,
                           ylab = dr_sm_axis(x$sm),
                           legend.title = "Knots", style = style,
                           legend = !isTRUE(result))
    what <- if (given) {
      paste0(length(fits), " of the ", length(sets), " knot sets given. ")
    } else {
      paste0(length(fits), " of the ", length(sets), " ", nknots,
             "-knot combinations of the ",
             paste0(format(100 * percentiles), collapse = ", "),
             "th percentiles. ")
    }
    foot <- paste0(what, "Curves that lie together say the result ",
                   "is the data; curves that fan apart say it is the knots.")
    # The table is the point of this figure, so it is set larger than the
    # summary tables under a single curve. It is also the legend, and the
    # knot set with the lowest AIC is set in bold.
    best <- if (any(is.finite(tab$AIC))) which.min(tab$AIC) else integer(0)
    tg <- if (isTRUE(result)) {
      dr_overlay_table(tab, style, best = best, note = foot)
    } else NULL
    spec <- dr_file_spec(file.name, filetype, path, width, height,
                         what = "easydr_knot_sensitivity()")
    dr_overlay_draw(p, tg, main, foot, spec)
  }

  structure(list(table = tab, fits = fits, knots = sets[ok],
                 sm_label = x$sm_label, unit = x$dose_unit),
            class = "easydr_knots")
}

#' Print a knot sensitivity analysis
#'
#' @param x an object from [easydr_knot_sensitivity()].
#' @param ... ignored.
#' @return `x`, invisibly. Called for the printed output.
#' @export
print.easydr_knots <- function(x, ...) {
  cat("Knot placement sensitivity (", x$sm_label, ")\n\n", sep = "")
  print(x$table, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' Keep only the test columns asked for
#'
#' The comparison table holds every test; a figure shows the ones switched on.
#'
#' @keywords internal
#' @noRd
dr_test_columns <- function(tab, wald = TRUE, deviance = TRUE, r2 = FALSE,
                            r2adj = FALSE, aic = FALSE, delta.aic = TRUE,
                            bic = FALSE, deviance.value = TRUE) {
  drop <- c(if (!isTRUE(wald)) "p nonlinearity",
            if (!isTRUE(deviance)) c("Deviance", "p deviance"),
            if (!isTRUE(deviance.value)) "Deviance",
            if (!isTRUE(r2)) "R-squared",
            if (!isTRUE(r2adj)) "R-squared adj",
            if (!isTRUE(aic)) "AIC",
            if (!isTRUE(delta.aic)) "Delta AIC",
            if (!isTRUE(bic)) "BIC")
  tab[, setdiff(names(tab), drop), drop = FALSE]
}

#' The name a reader knows a shape by
#'
#' The `shape` argument takes short codes; a table or legend a reader sees
#' spells them out.
#'
#' @keywords internal
#' @noRd
dr_shape_title <- function(x) {
  titles <- c(linear = "Linear", quadratic = "Quadratic", cubic = "Cubic",
              log = "Logarithmic", rcs = "Restricted cubic spline",
              fp = "Fractional polynomial", emax = "Emax",
              hyperbolic_emax = "Hyperbolic Emax")
  out <- unname(titles[x])
  ifelse(is.na(out), x, out)
}

#' The test lines under a single curve
#'
#' One line per test switched on, each from the fit's own numbers: the Wald
#' test from dosresmeta's non-linearity test, the deviance and R-squared from
#' dosresmeta::gof(), and the information criteria from the maximum
#' likelihood refit the comparison uses.
#'
#' @keywords internal
#' @noRd
dr_test_notes <- function(x, wald = TRUE, deviance = FALSE, r2 = FALSE,
                          r2adj = FALSE, aic = FALSE, bic = FALSE) {
  out <- character(0)
  if (isTRUE(wald)) out <- c(out, dr_linearity_note(x, lead = ""))
  g <- x$gof
  has_gof <- !is.null(g) && !identical(g$available, FALSE)
  if (isTRUE(deviance)) {
    out <- c(out, if (has_gof && !is.null(g$deviance)) {
      paste0("Deviance: D = ", formatC(g$deviance, format = "f", digits = 2),
             " (df = ", g$df, "), ", dr_p_lab(g$p))
    } else "Deviance: not available for this fit")
  }
  if (isTRUE(r2) || isTRUE(r2adj)) {
    r <- c(if (isTRUE(r2) && has_gof && !is.null(g$r2))
             paste0("R-squared = ", formatC(g$r2, format = "f", digits = 3)),
           if (isTRUE(r2adj) && has_gof && !is.null(g$r2adj))
             paste0("adjusted R-squared = ",
                    formatC(g$r2adj, format = "f", digits = 3)))
    out <- c(out, if (length(r)) paste(r, collapse = ", ") else
      "R-squared: not available for this fit")
  }
  if (isTRUE(aic) || isTRUE(bic)) {
    ml <- try(dr_refit_ml(x), silent = TRUE)
    if (inherits(ml, "try-error")) ml <- NULL
    ic <- c(if (isTRUE(aic)) paste0("AIC = ", formatC(dr_ic(x, "AIC", ml), format = "f", digits = 1)),
            if (isTRUE(bic)) paste0("BIC = ", formatC(dr_ic(x, "BIC", ml), format = "f", digits = 1)))
    out <- c(out, paste(ic, collapse = ", "))
  }
  out
}
