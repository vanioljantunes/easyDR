# Several shapes fitted to one table, held together.
#
# Comparing shapes is a normal part of a dose-response review: the protocol
# names a primary model and a secondary one, and the paper reports both. Doing
# that by hand means a list of fits, a loop and a layout, which is what the
# book used to show. This holds the fits instead, keeps them in the order they
# were named, and lets easydr_curve() draw them as panels.
#
# It is deliberately not an easydr_result. A forest, a funnel, a leave-one-out
# and a table each answer a question about one model, and quietly answering it
# for whichever shape happened to be named first would be worse than refusing.

#' Hold one fit per shape
#' @keywords internal
#' @noRd
dr_multi <- function(fits, shapes, call) {
  first <- fits[[1L]]
  structure(
    list(fits      = fits,
         shapes    = shapes,
         sm        = first$sm,
         sm_label  = first$sm_label,
         ratio     = first$ratio,
         dose_unit = first$dose_unit,
         reference = first$reference,
         counts    = first$counts,
         call      = call),
    class = "easydr_multi")
}

#' Refuse a multi fit, and say how to reach one of its fits
#'
#' The message names the shapes it holds and shows the expression that gets
#' one of them, because "takes a fitted result" on its own leaves the reader
#' with an object and no obvious way to open it.
#'
#' @keywords internal
#' @noRd
dr_require_single <- function(x, what) {
  if (inherits(x, "easydr_multi")) {
    stop("easydr: ", what, " reports one model, and this holds ",
         length(x$fits), " (", paste(x$shapes, collapse = ", "), "). ",
         "Pick the one you mean, such as x$fits$", x$shapes[1L],
         ". To see them side by side use ",
         "easydr_curve(x, mult.curves = TRUE).", call. = FALSE)
  }
  if (!inherits(x, "easydr_result")) {
    stop("easydr: ", what, " takes a fitted result from easydr().",
         call. = FALSE)
  }
  invisible(x)
}

#' The same fit, estimated by maximum likelihood
#'
#' Restricted maximum likelihood is the right default for estimating the
#' between-study variance, and the wrong basis for comparing models whose
#' fixed effects differ. A REML likelihood is computed on contrasts that
#' depend on the design matrix, so two shapes have REML likelihoods of
#' different quantities and the difference between them means nothing. The
#' comparison therefore refits under ML.
#'
#' The shape is handed back already resolved, so a profiled Emax is not
#' profiled a second time and the refit is the same curve.
#'
#' @keywords internal
#' @noRd
dr_refit_ml <- function(f) {
  if (identical(f$provenance$method, "ml")) return(f)
  inh <- try(dr_inherited(f), silent = TRUE)
  if (inherits(inh, "try-error")) return(NULL)
  args <- inh$args
  args$method.tau <- "ML"
  args$shape <- f$shape
  out <- try(suppressWarnings(
    do.call(inh$entry, c(list(data = inh$data), args))), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

#' AIC or BIC, charged for every parameter the fit actually spent
#'
#' The engine counts the coefficients it estimated. A profiled shape also
#' spent the parameters the profile searched over, and they were estimated
#' from the same data, so leaving them uncounted hands the Emax family a
#' discount no other shape gets and makes the comparison meaningless.
#'
#' @keywords internal
#' @noRd
dr_ic <- function(f, which = c("AIC", "BIC"), ml = NULL) {
  which <- match.arg(which)
  use   <- if (is.null(ml)) f else ml
  base  <- dr_try_num(if (which == "AIC") stats::AIC(use$fit) else
                        stats::BIC(use$fit))
  extra <- f$shape$profiled_params
  if (is.null(extra) || !is.finite(base)) return(base)
  penalty <- if (which == "AIC") 2 else log(max(f$counts$categories, 2))
  base + penalty * as.numeric(extra)
}

#' A number, or NA when the engine will not give one
#' @keywords internal
#' @noRd
dr_try_num <- function(expr) {
  v <- try(suppressWarnings(as.numeric(expr)), silent = TRUE)
  if (inherits(v, "try-error") || length(v) != 1L || !is.finite(v)) {
    return(NA_real_)
  }
  v
}

#' One row per shape: what each model costs and what it buys
#'
#' The columns are the ones a protocol compares shapes on. Terms is how many
#' coefficients the shape spends. Departure from linearity is the Wald test on
#' the non-linear terms, which a one-term shape has nothing to report for.
#' Deviance is the fit against the saturated model. AIC and BIC are the fit
#' with the terms charged for, which is what lets shapes that are not nested
#' inside one another be compared at all.
#'
#' @keywords internal
#' @noRd
dr_multi_table <- function(x, digits = 3) {
  rows <- lapply(x$shapes, function(s) {
    f  <- x$fits[[s]]
    nl <- f$nonlinearity
    # One refit per shape, reused by both information criteria.
    ml <- dr_refit_ml(f)
    data.frame(
      Shape = s,
      Terms = f$shape$terms,
      `p nonlinearity` = if (isTRUE(nl$estimable)) dr_p(nl$p) else "-",
      Deviance = if (is.null(f$gof$deviance)) NA_real_ else
        round(f$gof$deviance, 2),
      `R-squared` = if (is.null(f$gof$r2)) NA_real_ else round(f$gof$r2, 3),
      AIC = round(dr_ic(f, "AIC", ml), 1),
      BIC = round(dr_ic(f, "BIC", ml), 1),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Print several shapes fitted to the same table
#'
#' @param x an object from [easydr()] called with more than one `shape`.
#' @param digits significant digits for reported effects.
#' @param ... ignored.
#'
#' @return `x`, invisibly. Called for the printed output.
#'
#' @examples
#' data(alcohol_cancer)
#' easydr.object <- easydr(alcohol_cancer, sm = "RR",
#'                         shape = c("linear", "quadratic", "rcs"),
#'                         dose.unit = "g/day")
#' easydr.object
#'
#' # One of them, to go on with.
#' easydr.object$fits$rcs
#'
#' @export
print.easydr_multi <- function(x, digits = 3, ...) {
  cat(length(x$fits), " shapes fitted to the same table (",
      x$sm_label, ")\n", sep = "")
  unit <- if (is.null(x$dose_unit)) "dose unit" else x$dose_unit
  cat(x$counts$studies, " studies, ", x$counts$categories,
      " dose categories, dose range ", x$counts$dose_min, " to ",
      x$counts$dose_max, " ", unit, "\n\n", sep = "")

  print(dr_multi_table(x, digits), row.names = FALSE, right = FALSE)

  cat("\nDeviance scores the fitted curve against a saturated model and",
      "\njudges the shape alone. AIC and BIC also pay for the",
      "\nbetween-study variance, and are computed by maximum",
      "\nlikelihood, because a REML likelihood depends on the design",
      "\nand two shapes would not be comparable. Both charge for every",
      "\nparameter the shape spent, including the ones an Emax profile",
      "\nsearched over.\n")
  cat("\nThe two can disagree: a shape can track the doses closely and",
      "\nstill leave the studies far apart. Say in the protocol which",
      "\nof them decides, before seeing either.\n")
  cat("\nOne fit: x$fits$", x$shapes[1L], "\n", sep = "")
  cat("Side by side: easydr_curve(x, mult.curves = TRUE)\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Several curves in one figure
# ---------------------------------------------------------------------------

#' How many rows and columns for n panels
#'
#' Fixed by the author rather than derived: two and three go side by side,
#' four is a square, five is three over two, six is two rows of three. Beyond
#' six there is no layout, and there are only six shapes to fit anyway.
#'
#' @keywords internal
#' @noRd
dr_panel_grid <- function(n) {
  if (n < 1L) {
    stop("easydr: there is nothing to draw.", call. = FALSE)
  }
  if (n > 6L) {
    stop("easydr: ", n, " panels is more than this layout covers. Draw at ",
         "most six, naming them with curve.panels.", call. = FALSE)
  }
  switch(as.character(n),
         "1" = list(nrow = 1L, ncol = 1L),
         "2" = list(nrow = 1L, ncol = 2L),
         "3" = list(nrow = 1L, ncol = 3L),
         "4" = list(nrow = 2L, ncol = 2L),
         # Three on the first row and two on the second, so the last cell is
         # left empty rather than stretching one panel across it.
         "5" = list(nrow = 2L, ncol = 3L,
                    layout = rbind(c(1, 2, 3), c(4, 5, NA))),
         "6" = list(nrow = 2L, ncol = 3L))
}

#' The dose and effect range every panel should share
#'
#' Panels drawn on their own scales cannot be compared, which is the only
#' reason to put them side by side. The range is the union over the shapes
#' being drawn, padded the way a single curve pads its own.
#'
#' @keywords internal
#' @noRd
dr_panel_limits <- function(fits, ratio) {
  curves <- lapply(fits, function(f) f$curve)
  xs <- range(unlist(lapply(curves, function(d) range(d$dose))))
  ys <- range(unlist(lapply(curves, function(d)
    range(c(d$est, d$ci.lb, d$ci.ub), na.rm = TRUE))), na.rm = TRUE)
  list(xlim = xs, ylim = dr_pad_range(ys, ratio))
}

#' One column per shape, under the panels
#'
#' Every column reports the same quantity at the same dose, or the columns
#' would not be comparable. A one-term shape left to itself reports its slope
#' and a spline reports an effect at a dose, so a common dose is chosen and
#' every shape is asked for that.
#'
#' @keywords internal
#' @noRd
dr_panel_table <- function(fits, shapes, at = NULL, digits = 3) {
  first <- fits[[1L]]
  if (is.null(at)) {
    cand <- first$predicted$dose[!first$predicted$is_ref]
    at <- if (length(cand) > 0L) cand[1L] else max(first$frame$dose)
  }
  quant <- c(dr_pooled_at(first, at)$label,
             paste0(format(100 * first$provenance$level), "% CI"),
             "p value", "I2", "Q (df, p)", "Studies")
  out <- data.frame(quant, check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- " "
  for (s in shapes) out[[s]] <- dr_result_column(fits[[s]], at, digits)
  out
}

#' Draw one panel per shape, on shared axes, with one table underneath
#'
#' @keywords internal
#' @noRd
dr_curve_multi <- function(x, mult.curves = FALSE, curve.panels = NULL,
                           main = NULL, xlim = NULL, ylim = NULL,
                           result = TRUE, result.at = NULL,
                           digits = 3, filetype = c("png", "pdf", "tiff"),
                           path = NULL, file.name = NULL, width = NULL,
                           height = NULL, ...) {
  if (!inherits(x, "easydr_multi")) {
    stop("easydr: mult.curves draws one panel per shape, and this fit holds ",
         "one. Fit several, such as shape = c(\"linear\", \"rcs\"), then ",
         "draw them together.", call. = FALSE)
  }
  if (!isTRUE(mult.curves)) {
    stop("easydr: this holds ", length(x$fits), " shapes (",
         paste(x$shapes, collapse = ", "), "). Set mult.curves = TRUE to ",
         "draw them side by side, or pick one with x$fits$", x$shapes[1L],
         ".", call. = FALSE)
  }

  panels <- if (is.null(curve.panels)) x$shapes else curve.panels
  if (!is.character(panels) || !length(panels)) {
    stop("easydr: `curve.panels` names the shapes to draw, in the order to ",
         "draw them, such as curve.panels = c(\"rcs\", \"linear\").",
         call. = FALSE)
  }
  missing <- setdiff(panels, x$shapes)
  if (length(missing)) {
    stop("easydr: `curve.panels` names ", paste(missing, collapse = ", "),
         ", which ", if (length(missing) == 1L) "is" else "are",
         " not fitted here. This holds ", paste(x$shapes, collapse = ", "),
         ".", call. = FALSE)
  }
  grid_spec <- dr_panel_grid(length(panels))
  fits <- x$fits[panels]
  # Shared limits unless the caller named their own, which they may want when
  # a figure has to match one already published.
  lims <- dr_panel_limits(fits, x$ratio)
  if (!is.null(xlim)) lims$xlim <- xlim
  if (!is.null(ylim)) lims$ylim <- ylim

  # Each panel is titled with its own shape and carries no footer of its own:
  # the figure has one footer, under the whole thing.
  plots <- lapply(panels, function(s) {
    built <- dr_curve_panel(fits[[s]], main = FALSE,
                            xlim = lims$xlim, ylim = lims$ylim, ...)
    built$plot +
      ggplot2::ggtitle(s) +
      ggplot2::theme(plot.title = ggplot2::element_text(
        size = 10, face = "bold", hjust = 0.5))
  })

  body <- if (is.null(grid_spec$layout)) {
    gridExtra::arrangeGrob(grobs = plots, nrow = grid_spec$nrow,
                           ncol = grid_spec$ncol)
  } else {
    gridExtra::arrangeGrob(grobs = plots, layout_matrix = grid_spec$layout)
  }

  if (is.null(main)) {
    main <- paste0(x$sm_label, " by shape, same studies")
  }
  if (identical(main, FALSE)) main <- NULL

  unit <- if (is.null(x$dose_unit)) "dose" else x$dose_unit
  foot <- paste0("Reference ", x$reference$curve, " ", unit,
                 ".  Bands: ", format(100 * fits[[1L]]$provenance$level),
                 "% CI.  All panels share one dose and one effect axis, so ",
                 "the shapes can be compared by eye.")

  tg <- if (isTRUE(result)) {
    dr_stats_grob(dr_panel_table(fits, panels, result.at, digits))
  } else {
    NULL
  }
  fs <- dr_footer_stack(foot, tg)

  # Width grows with the columns and height with the rows, and the table and
  # title take their own measured space on top of that.
  spec <- dr_file_spec(file.name, filetype, path, width, height,
                       what = "easydr_curve(mult.curves = TRUE)")
  title_in <- if (is.null(main)) 0 else 1.9 * 12 / 72
  foot_in  <- dr_height_in(fs$height)
  w_in     <- max(3.6 * grid_spec$ncol, dr_grob_width_in(tg) + 0.5)
  size <- dr_fig_size(w_in, 2.9 * grid_spec$nrow + title_in + foot_in + 0.2,
                      max.w = 18, max.h = 24)
  con  <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)

  dr_draw_stack(body, NULL, title = main, footer = fs$grob,
                footer_h = fs$height, title_hjust = 0.5)
  invisible(list(shapes = panels, table = dr_multi_table(x, digits),
                 layout = grid_spec))
}
