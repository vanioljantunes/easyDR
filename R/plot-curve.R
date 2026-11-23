#' Plot the fitted dose-response curve
#'
#' @description
#' The pooled dose-response against dose, with its confidence band, optionally
#' its prediction band, and the doses the fit reports marked and labelled.
#' Every number was computed by [easydr()]; this chooses what appears and how.
#'
#' @section What the two bands mean:
#' The confidence band describes where the pooled curve lies. The prediction
#' band describes where a new study's curve would fall, which makes it wider.
#' A common-effect fit has no prediction band at all, because it assumes there
#' is no between-study variation for one to describe. Showing the prediction
#' band and calling it uncertainty about the pooled curve is the usual
#' mistake.
#'
#' Doses beyond the observed range are drawn, because [easydr()] will predict
#' there if asked, but they are marked: the model has no data to support that
#' part of the curve.
#'
#' @section Reading the axis:
#' A ratio measure goes on a log axis, so that a halving and a doubling sit
#' the same distance from the null. The dashed horizontal is that null, at 1
#' for a ratio and 0 for a difference. The dotted vertical marks the reference
#' dose every effect is measured against.
#'
#' @param x a fitted result from [easydr()].
#' @param main plot title. `NULL` writes one naming the shape; `FALSE` draws
#'   none.
#' @param xlab,ylab axis labels. `NULL` writes them from the dose unit and the
#'   summary measure, such as "Risk ratio (RR)".
#' @param xlim,ylim axis limits on the reported scale, for example
#'   `ylim = c(0.8, 2)`. `NULL` fits them to whatever is drawn.
#' @param round digits for the labels printed beside the marked doses. `NULL`
#'   marks the doses without labelling them.
#' @param prediction `TRUE` adds the prediction band when the fit has one.
#'   `FALSE`, the default, draws the confidence band alone.
#' @param points `TRUE`, the default, marks the doses the fit reports, which
#'   are the ones named in `at`.
#' @param studies `TRUE`, the default, draws the evidence behind the curve:
#'   one point per reported dose category, at that dose on the x axis and at
#'   the effect the study reported for it against its own reference dose on
#'   the y axis. Each point's area follows the participants behind it: the
#'   category's own and those of its study's reference category, with a
#'   floor so the smallest stays visible. The reference
#'   categories are not drawn, because every one of them sits on the null by
#'   construction. The curve and its marked doses are drawn over the points,
#'   so where they meet the curve stays readable. With several panels
#'   (`mult.curves = TRUE`) the default is `FALSE`, because the points crowd
#'   the small panels.
#' @param ref.shape `TRUE`, the default, gives the study points one shape per
#'   reference dose (circle, square, triangle, diamond, inverted triangle,
#'   star), with a "Study reference dose (unit)" key in one row under the
#'   panel listing each value. A study's effects are measured against its own reference, so
#'   a point is read against that dose, not the curve's. `FALSE` draws every
#'   point as a circle. It is ignored when subgroup curves are drawn, which use
#'   shape for the levels, with more than six reference doses, which run out
#'   of shapes, and when every study shares one reference dose, where a key
#'   would say nothing.
#' @param studlab `FALSE`, the default, leaves the study names off. `TRUE`
#'   writes each study's name once, beside its highest-dose point, so a reader
#'   can tell which study a run of points belongs to. It needs
#'   `studies = TRUE`.
#' @param reference `TRUE`, the default, draws the null line and marks the
#'   reference dose.
#' @param rug `TRUE` ticks the doses the fit reports along the bottom: the
#'   ones `at` named in [easydr()], and the reference. It marks where the
#'   reported effects sit on the axis, not where the data is. The default,
#'   `NULL`, draws it unless knot lines are drawn instead.
#' @param knot.lines `TRUE`, the default, draws a dashed vertical line at each
#'   knot of a restricted cubic spline, in place of the rug. Ignored for
#'   shapes without knots.
#' @param grid `FALSE`, the default, draws the panel without grid lines.
#'   `TRUE` adds the major grid.
#' @param subtitle `TRUE`, the default, writes three lines under the title,
#'   items separated by " · ": the shape, the number of studies and of
#'   participants and person-years; the headline effect in bold with its
#'   confidence interval and p value; and the model's performance (I², the
#'   Q test, the linearity and deviance tests as p values, and the adjusted
#'   R², see `wald`, `deviance` and `r2adj`). For an increment of dose the
#'   headline is the effect per that increment: the slope of a linear shape,
#'   and for any other shape the linear trend of the same studies, refitted
#'   as a straight line the way [easydr_loo()] does, labelled "(linear
#'   trend)". For a dose contrast, `result.at = c(a, b)`, it is the curve's
#'   own effect at `a` against `b`. With subgroups drawn it reports the
#'   overall fit.
#' @param results.table the doses to report in the table under the figure,
#'   each against the reference, such as `results.table = c(2, 4, 6)`. They
#'   are also marked on the curve and labelled on the dose axis. A dose above
#'   the highest one observed extends the curve, dashed, with a fainter band.
#'   `NULL`, the default, uses the doses `at` named in [easydr()].
#' @param result `TRUE`, the default, puts a summary table under the figure:
#'   the effect with its confidence interval, its p value, the heterogeneity
#'   and the number of studies. The figure then carries its own numbers. With
#'   subgroups drawn, the table gains one column per level and ends with the
#'   overall fit.
#' @param wald,deviance,r2,r2adj,aic,bic `TRUE` (`T`) or `FALSE` (`F`): add
#'   or leave out a line under the panel for that test. `wald` is the
#'   linearity test (on by default), `deviance` the deviance against a
#'   saturated model with its p value (on by default), `r2` the R-squared,
#'   `r2adj` (on by default, in the subtitle beside the deviance) the
#'   adjusted R-squared, `aic` and `bic` the information criteria of this fit,
#'   from its one-stage maximum likelihood refit.
#' @param target `TRUE` puts the target doses (ED50, ED80, ED max) in a box in
#'   the corner of the panel, for an Emax fit, which computes them. Ignored
#'   for other shapes.
#' @param result.at what the table reports. `NULL` (the default) uses the
#'   `result.at` given to [easydr()], which is per unit of dose unless set
#'   there; one number (an increment of dose) or `c(numerator, denominator)`
#'   overrides it for this plot.
#' @param table `"horizontal"` (the default) puts the quantities across the
#'   top of the summary table, `"vertical"` puts them down its first column.
#' @param subgroup `TRUE` draws one curve per subgroup level, `FALSE` draws the
#'   overall fit alone. The default, `NULL`, draws the levels when the fit
#'   carries them. Each level is told apart by colour and by line type, and its
#'   marked doses by shape. The band stays the overall one. `TRUE` on a fit
#'   with no subgroup is an error rather than a silent overall curve; see the
#'   `subgroup` argument of [easydr()].
#' @param digits digits for the numbers in that table.
#' @param col,fill.ci,fill.pi colours for the curve, the confidence band and
#'   the prediction band.
#' @param footer `TRUE`, the default, writes the line under the panel that
#'   explains the band and the points. `FALSE` leaves it out, for small
#'   tiles where it cannot be read.
#' @param legend optional free text drawn under the figure, one line per
#'   element.
#'
#' @return The plotted curve table, invisibly: `dose`, `est`, `ci.lb`,
#'   `ci.ub`, `pi.lb`, `pi.ub` and `extrapolated`.
#'
#' @seealso [easydr_forest()] for what each study contributed.
#'
#' @examples
#' easydr.object <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
#'                         dose.unit = "g/day", at = c(10, 25, 50))
#'
#' # The defaults: confidence band, reported doses marked, null line, and a
#' # rug of the exposure categories the studies actually reported.
#' easydr_curve(easydr.object)
#'
#' # Titled, labelled to two decimals, axis pinned, prediction band shown.
#' easydr_curve(easydr.object,
#'              main  = "Alcohol and colorectal cancer",
#'              xlab  = "Alcohol intake (g/day)",
#'              ylab  = "Rate ratio",
#'              ylim  = c(0.9, 1.8),
#'              round = 2,
#'              prediction = TRUE)
#'
#' @param filetype the format a saved figure is written in: `"png"` (the
#'   default), `"pdf"` or `"tiff"`. Only used when `file.name` is given.
#' @param path the directory to write the figure into. `NULL` is the working
#'   directory.
#' @param file.name the name of the file to write. `NULL`, the default, draws
#'   on the current device and writes nothing. The extension is added when it
#'   is missing.
#' @param width,height the size of the saved figure in inches. `NULL` lets
#'   easyDR measure the figure and choose a size that fits it.
#' @param mult.curves `TRUE` draws one panel per shape when `x` holds several,
#'   on shared axes so the shapes can be compared by eye. Two and three
#'   panels go side by side, four is a square, five is three over two, and
#'   six is two rows of three.
#' @param curve.panels the shapes to draw, in the order to draw them, such as
#'   `curve.panels = c("rcs", "linear")`. `NULL` draws every shape the fit
#'   holds, in the order they were named to [easydr()].
#' @export
easydr_curve <- function(x, main = NULL, xlab = NULL, ylab = NULL,
                         xlim = NULL, ylim = NULL, round = NULL,
                         prediction = FALSE, points = TRUE, studies = TRUE,
                         ref.shape = TRUE, studlab = FALSE, reference = TRUE,
                         rug = NULL, knot.lines = TRUE, grid = FALSE,
                         subtitle = TRUE, footer = TRUE,
                         results.table = NULL,
                         result = TRUE, result.at = NULL,
                         table = "horizontal", target = TRUE,
                         wald = TRUE, deviance = TRUE, r2 = FALSE,
                         r2adj = TRUE, aic = FALSE, bic = FALSE,
                         subgroup = NULL, digits = 3, col = "#1F2A44",
                         fill.ci = "#8C9AB3", fill.pi = "grey90",
                         legend = NULL,
                         filetype = c("png", "pdf", "tiff"),
                         path = NULL, file.name = NULL,
                         width = NULL, height = NULL,
                         mult.curves = FALSE,
                         curve.panels = NULL) {
  # Several shapes go to the panel drawer, which builds one of these plots
  # per shape and lays them out. A single fit asked for panels is a mistake
  # worth naming rather than quietly drawing one panel.
  result.at <- dr_result_at(x, result.at)
  if (inherits(x, "easydr_multi") || isTRUE(mult.curves)) {
    # Small panels side by side have no room for the study points, so they
    # are left out unless asked for.
    if (missing(studies)) studies <- FALSE
    return(dr_curve_multi(
      x, mult.curves = mult.curves, curve.panels = curve.panels,
      main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
      round = round, prediction = prediction, points = points,
      studies = studies, ref.shape = ref.shape, studlab = studlab,
      reference = reference,
      rug = rug, knot.lines = knot.lines, grid = grid,
      result = result, result.at = result.at, digits = digits,
      col = col, fill.ci = fill.ci, fill.pi = fill.pi, legend = legend,
      filetype = filetype, path = path, file.name = file.name,
      width = width, height = height))
  }
  dr_require_single(x, "easydr_curve()")
  xref <- x$reference$curve
  rt <- if (is.null(results.table)) {
    x$predicted$dose[!x$predicted$is_ref]
  } else {
    if (!is.numeric(results.table) || anyNA(results.table)) {
      stop("easydr: `results.table` takes the doses to report, such as ",
           "results.table = c(2, 4, 6).", call. = FALSE)
    }
    sort(unique(as.numeric(results.table)))
  }
  rt <- rt[rt != xref]

  built   <- dr_curve_panel(
    x, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
    round = round, prediction = prediction, points = points,
    studies = studies, ref.shape = ref.shape, studlab = studlab,
    reference = reference, mark.doses = rt,
    axis_in = if (is.null(width)) 7.2 else width,
    rug = rug, knot.lines = knot.lines, grid = grid,
    subgroup = subgroup, col = col, fill.ci = fill.ci,
    fill.pi = fill.pi, legend = legend, target = target,
    tests = list(wald = wald, deviance = deviance, r2 = r2, r2adj = r2adj,
                 aic = aic, bic = bic))
  # The key sits right of the panel, which keeps its place on the left.
  p       <- dr_centre_legend(built$plot, centre = isTRUE(built$use_sub))
  marg    <- attr(p, "panel_margins")
  foot    <- if (isTRUE(footer)) built$foot else NULL
  tab     <- built$tab
  use_sub <- built$use_sub
  main    <- built$main

  spec <- dr_file_spec(file.name, filetype, path, width, height,
                       what = "easydr_curve()")
  # The summary sits under the figure, below the footer, the way the funnel
  # puts Egger's numbers under its panel. Inside the panel it covered the
  # curve on a fit that rises into its own corner, and it could hold only one
  # column, which left nowhere for a subgroup to go.
  # The effects at chosen doses sit under the figure; the effect per
  # increment is in the subtitle above it.
  tg <- if (isTRUE(result) && isTRUE(use_sub)) {
    # By subgroup the table reports what the levels are compared on: the
    # effect per increment, then heterogeneity and fit, one column per level.
    dr_curve_table_grob(dr_subgroup_summary(x, result.at, digits),
                        title = paste0("Results by ", x$subgroup$var),
                        note = "")
  } else if (isTRUE(result) && length(rt)) {
    dt <- dr_dose_table(x, rt, digits, use_sub)
    if (identical(match.arg(table, c("horizontal", "vertical")), "vertical")) {
      m <- as.data.frame(t(as.matrix(dt[, -1L, drop = FALSE])),
                         stringsAsFactors = FALSE)
      names(m) <- paste(names(dt)[1L], dt[[1L]])
      dt <- cbind(data.frame(` ` = names(dt)[-1L], check.names = FALSE,
                             stringsAsFactors = FALSE), m)
    }
    if (identical(match.arg(table, c("horizontal", "vertical")), "vertical")) {
      g <- dr_stats_grob(dt)
      keep <- g$layout$name != "top-rule"
      g$grobs  <- g$grobs[keep]
      g$layout <- g$layout[keep, , drop = FALSE]
      g
    } else {
      ci_col <- paste0(format(100 * x$provenance$level), "% CI")
      if (identical(names(dt)[2L], x$sm)) names(dt)[2L] <- dr_curve_ylab(x$sm)
      if (ci_col %in% names(dt)) {
        v <- dt[[ci_col]]
        if (!any(grepl("-", v))) dt[[ci_col]] <- gsub(" to ", "\u2013", v)
        names(dt)[names(dt) == ci_col] <- paste0(format(100 * x$provenance$level),
                                                 "% confidence interval")
      }
      unit_r <- if (is.null(x$dose_unit)) "" else paste0(" ", x$dose_unit)
      dr_curve_table_grob(
        dt, title = paste0("Predicted ", dr_sm_label(x$sm), "s by dose"),
        note = paste0("Compared with ", dr_num_lab(xref), unit_r))
    }
  } else {
    NULL
  }
  w_in <- max(7.2, if (is.null(tg)) 0 else if (!is.null(attr(tg, "natural_in")))
    attr(tg, "natural_in") else dr_grob_width_in(tg) + 0.5)

  # The footer lines start where the panel box starts and wrap at its right
  # edge, so the text reads as belonging to the plot above it.
  dev_w <- if (!is.null(spec) && !is.null(spec$width)) {
    spec$width
  } else if (is.null(spec) && !is.null(grDevices::dev.list())) {
    grDevices::dev.size("in")[1L]
  } else {
    dr_fig_size(w_in, 1)[["width"]]
  }
  foot_x <- grid::unit(4, "pt")
  foot_r <- grid::unit(4, "pt")
  panel_cm <- dev_w * 2.54 * 0.8
  if (length(marg) == 2L && all(is.finite(marg))) {
    foot_x   <- grid::unit(marg[["left"]], "cm")
    foot_r   <- grid::unit(marg[["right"]], "cm")
    panel_cm <- dev_w * 2.54 - sum(marg)
    if (length(foot)) foot <- dr_wrap_text(foot, dev_w * 2.54 - 1.2, 7.5)
  }
  # Caption, a thin rule, then the table: one block under the panel, so the
  # numbers read as part of the figure rather than as a second one.
  fs <- dr_curve_footer(foot, tg, text_x = foot_x, right = foot_r,
                        centre = TRUE, size = 7.5)

  # The panel is the figure here, and the band under it is whatever the
  # footer needs. The width follows the summary table when there is one,
  # because a subgroup table with four columns is wider than a plot of a
  # single curve wants to be.
  foot_in  <- if (is.null(fs)) 0 else dr_height_in(fs$height)
  sub_g    <- if (isTRUE(subtitle)) {
    # What the curve rests on, small and grey; then, after a gap, the result
    # and the tests, the values in bold.
    st <- function(size, col, gap) list(size = size, col = col,
                                        strong_col = "black", sep_col = "grey70",
                                        gap = gap)
    dr_rich_grob(dr_rich_lines(
      dr_curve_subtitle_tokens(x, result.at, digits, wald = wald,
                               deviance = deviance, subgroup = use_sub,
                               r2adj = r2adj),
      list(st(7.5, "grey40", 0), st(8.5, "grey20", 6), st(8.5, "grey20", 1)),
      dev_w * 2.54 - 1))
  } else {
    NULL
  }
  title_in <- (if (is.null(main)) 0 else 1.9 * 12 / 72) +
    (if (is.null(sub_g)) 0 else dr_height_in(attr(sub_g, "height")))
  size <- dr_fig_size(w_in,
                      min(max(0.50 * w_in, 3.1), 4.4) + foot_in + title_in +
                        0.12)
  con  <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)

  dr_draw_stack(p, NULL, title = main,
                footer = if (is.null(fs)) NULL else fs$grob,
                footer_h = if (is.null(fs)) NULL else fs$height,
                title_hjust = 0.5, subtitle = sub_g)
  invisible(tab)
}

#' The block under a curve: caption, a thin rule, and the results table
#'
#' The caption is small and grey, because it explains marks rather than
#' reporting results. The rule spans the panel, not the device, so the table
#' reads as belonging to the plot above it, and the table follows with only a
#' few points of space.
#'
#' @keywords internal
#' @noRd
dr_curve_footer <- function(text, table_grob, size = 8.5,
                            text_x = grid::unit(4, "pt"),
                            right = grid::unit(4, "pt"), centre = FALSE) {
  # Centred under the figure, the caption and the table span its width.
  if (isTRUE(centre)) {
    text_x <- grid::unit(0.22, "in")
    right  <- grid::unit(0.22, "in")
  }
  if (!length(text) && is.null(table_grob)) return(NULL)
  grobs   <- list()
  heights <- list()
  if (length(text)) {
    n <- length(strsplit(paste(text, collapse = "\n"), "\n", fixed = TRUE)[[1]])
    grobs <- c(grobs, list(grid::textGrob(
      paste(text, collapse = "\n"),
      x = if (isTRUE(centre)) grid::unit(0.5, "npc") else text_x,
      y = grid::unit(1, "npc") - grid::unit(3, "pt"),
      hjust = if (isTRUE(centre)) 0.5 else 0, vjust = 1,
      gp = grid::gpar(fontsize = size, lineheight = 1.2, col = "grey35"))))
    heights <- c(heights, list(grid::unit(
      1.2 * size * n + if (is.null(table_grob)) 10 else 14, "points")))
  }
  if (!is.null(table_grob)) {
    rule <- grid::segmentsGrob(
      x0 = text_x, x1 = grid::unit(1, "npc") - right,
      y0 = grid::unit(0.5, "npc"), y1 = grid::unit(0.5, "npc"),
      gp = grid::gpar(col = "grey80", lwd = 0.6))
    grobs   <- c(grobs, list(rule, table_grob))
    heights <- c(heights, list(grid::unit(8, "points"),
                               dr_grob_height(table_grob, pad = 10)))
  }
  hv <- do.call(grid::unit.c, heights)
  list(grob = gridExtra::arrangeGrob(grobs = grobs, ncol = 1, heights = hv),
       height = sum(hv))
}

#' The axis title for a summary measure, in sentence case
#'
#' "Risk ratio (RR)", the way a journal writes it, rather than the title case
#' the other figures use.
#'
#' @keywords internal
#' @noRd
dr_curve_ylab <- function(sm) {
  lab <- dr_sm_label(sm)
  if (identical(lab, sm)) return(sm)
  paste0(toupper(substring(lab, 1, 1)), substring(lab, 2), " (", sm, ")")
}

#' The effect line the curve's subtitle leads with
#'
#' A straight line has one slope, so an increment of dose reports it. A curve
#' has none, so an increment asked of a curve is answered with the linear
#' trend of the same studies, refitted the way easydr_loo() does it, and the
#' label says so. A dose contrast is the curve's own prediction at the first
#' dose against the second. When the fit cannot be refitted (a result made
#' before refitting was supported) the curve's own contrast is reported.
#'
#' @keywords internal
#' @noRd
dr_curve_headline <- function(x, at = NULL, digits = 3) {
  spec <- dr_result_at(x, at)
  unit <- if (is.null(x$dose_unit)) "dose unit" else x$dose_unit
  src  <- x
  lab  <- NULL
  if (dr_wants_trend(x, spec)) {
    trend <- tryCatch({
      inh <- dr_inherited(x)
      res <- NULL
      utils::capture.output(res <- suppressMessages(suppressWarnings(
        dr_refit_trend(inh, inh$studies))))
      res
    }, error = function(e) NULL)
    if (!is.null(trend)) {
      src <- trend
      lab <- paste0(dr_label(dr_unit_one(unit, spec$value), spec$value),
                    " (linear trend)")
    }
  }
  p <- dr_pooled_at(src, spec)
  if (is.null(lab)) {
    lab <- if (identical(spec$type, "pair")) {
      paste0("At ", format(spec$value[1L]), " ",
             dr_unit_one(unit, spec$value[1L]), " vs ", format(spec$value[2L]))
    } else {
      dr_label(dr_unit_one(unit, spec$value), spec$value)
    }
  }
  back <- if (x$ratio) exp else identity
  fmt  <- function(v) format(round(back(v), digits), nsmall = digits)
  lvl  <- format(100 * x$provenance$level)
  lo <- fmt(p$lb)
  hi <- fmt(p$ub)
  # An en dash between the limits, as journals print them, unless a limit is
  # negative, where "-0.2-0.1" would misread.
  ci <- if (grepl("^-", lo) || grepl("^-", hi)) {
    paste0(lo, " to ", hi)
  } else {
    paste0(lo, "\u2013", hi)
  }
  list(label = lab,
       est = paste0(x$sm, " ", fmt(p$est)),
       ci = paste0("(", lvl, "% CI ", ci, "),"),
       p = dr_p_lab(p$p),
       effect = paste0(x$sm, " ", fmt(p$est), " (", lvl, "% CI ", ci, "), ",
                       dr_p_lab(p$p)),
       pooled = p, trend = !identical(src, x))
}

#' A unit written for one of it: "1 cup/day", not "1 cups/day"
#'
#' Only the first word before a slash loses a plural "s", and only when the
#' amount is exactly one.
#'
#' @keywords internal
#' @noRd
dr_unit_one <- function(unit, value) {
  if (length(value) && isTRUE(all.equal(value[1L], 1))) {
    sub("^([A-Za-z]+[^s])s(/|$)", "\\1\\2", unit)
  } else {
    unit
  }
}

#' The subtitle under a curve's title, as tokens
#'
#' Three lines: what the curve is and what it rests on, the headline effect
#' in bold, and the model's performance (heterogeneity and the tests, each as
#' its p value). Items are separated by " · " and kept whole when a line
#' wraps, so a test never splits from its p value.
#'
#' @keywords internal
#' @noRd
dr_curve_subtitle_tokens <- function(x, at = NULL, digits = 3, wald = TRUE,
                                     deviance = TRUE, subgroup = FALSE,
                                     r2adj = FALSE) {
  pv  <- dr_p_lab
  # Items are parted by a light dot with room on both sides. An item is a
  # plain label and, where it has one, a bold value.
  dot <- list(text = "\u2002\u00b7\u2002", role = "sep", space = FALSE)
  sep <- function(items) {
    out <- list()
    for (i in seq_along(items)) {
      it <- items[[i]]
      if (!is.list(it)) it <- list(it)
      if (i > 1L) out <- c(out, list(dot))
      out <- c(out, list(list(text = it[[1L]], space = FALSE)))
      if (length(it) > 1L) {
        out <- c(out, list(list(text = it[[2L]], strong = TRUE, space = TRUE)))
      }
    }
    out
  }
  # What the curve is and what it rests on. Counts of people and of
  # person-time are different things, so a fit mixing incidence-rate studies
  # with the others reports both rather than adding them up.
  fr  <- x$frame
  big <- function(v) if (v >= 1e6) {
    paste0(formatC(v / 1e6, format = "fg", digits = 3), " million")
  } else {
    format(round(v), big.mark = ",", scientific = FALSE, trim = TRUE)
  }
  is_pt <- fr$type %in% "ir"
  nm <- sub(",.*$", "", x$shape$label)
  nm <- sub(" in dose$", "", nm)
  nm <- paste0(toupper(substring(nm, 1, 1)), substring(nm, 2))
  info <- c(paste(nm, "dose\u2013response meta-analysis"),
            paste(x$counts$studies, "studies"),
            if (any(!is_pt) && all(is.finite(fr$n[!is_pt])))
              paste(big(sum(fr$n[!is_pt])), "participants"),
            if (any(is_pt) && all(is.finite(fr$n[is_pt])))
              paste(big(sum(fr$n[is_pt])), "person-years"))

  hl <- dr_curve_headline(x, at, digits)

  het  <- x$heterogeneity
  perf <- if (isTRUE(het$estimable)) {
    list(list("I\u00b2 =", sprintf("%.1f%%", het$i2)), list("Q test", pv(het$p)))
  } else {
    list("I\u00b2 not estimable")
  }
  if (isTRUE(wald) && !identical(x$shape$name, "linear")) {
    nl <- x$nonlinearity
    perf <- c(perf, list(if (!is.null(nl) && isTRUE(nl$estimable)) {
      list("Non-linearity", pv(nl$p))
    } else "Non-linearity test not available"))
  }
  g  <- x$gof
  ok <- !is.null(g) && !identical(g$available, FALSE)
  if (isTRUE(deviance)) {
    perf <- c(perf, list(if (ok && !is.null(g$p) && is.finite(g$p)) {
      list("Deviance", pv(g$p))
    } else "Deviance test not available"))
  }
  # Adjusted R-squared (Discacciati, 2017): how much of the variation in the
  # reported effects the curve explains, charged for its number of terms.
  if (isTRUE(r2adj)) {
    perf <- c(perf, list(if (ok && !is.null(g$r2adj) && is.finite(g$r2adj)) {
      list("Adjusted R\u00b2 =", formatC(g$r2adj, format = "f", digits = 2))
    } else "Adjusted R\u00b2 not available"))
  }
  # A curve drawn by subgroup carries the test for subgroup interaction
  # (Berlin, 1993), the comparison its several lines invite.
  if (isTRUE(subgroup)) {
    t <- x$subgroup$test
    perf <- c(perf, list(if (!is.null(t) && isTRUE(t$estimable)) {
      list("Subgroup interaction", pv(t$p))
    } else "Subgroup interaction test not available"))
  }

  list(sep(as.list(info)),
       list(list(text = paste0(hl$label, ":"), space = FALSE),
            list(text = hl$est, strong = TRUE, space = TRUE),
            list(text = hl$ci, space = TRUE),
            list(text = hl$p, strong = TRUE, space = TRUE)),
       sep(perf))
}

#' Build the curve panel, without drawing it
#'
#' Split out of `easydr_curve()` so that several shapes can be laid out side
#' by side: the panel drawer needs the ggplot object, and the drawing step
#' needs to happen once for the whole figure rather than once per panel.
#'
#' @keywords internal
#' @noRd
dr_curve_panel <- function(x, main = NULL, xlab = NULL, ylab = NULL,
                           xlim = NULL, ylim = NULL, round = NULL,
                           prediction = FALSE, points = TRUE, studies = TRUE,
                           ref.shape = TRUE, studlab = FALSE,
                           reference = TRUE, rug = NULL,
                           knot.lines = TRUE, grid = FALSE, mark.doses = NULL,
                           axis_in = 7.2,
                           subgroup = NULL, col = "#1F2A44",
                           fill.ci = "#8C9AB3", fill.pi = "grey90",
                           legend = NULL, target = TRUE,
                           tests = list(wald = TRUE)) {
  use_sub <- dr_use_subgroups(x, subgroup, "easydr_curve()")
  tab  <- x$curve
  xref <- x$reference$curve
  if (is.null(mark.doses)) {
    mark.doses <- x$predicted$dose[!x$predicted$is_ref]
  }
  # A dose asked for above the highest observed one carries the curve on to
  # it, so the reported effect has a curve to sit on.
  top <- max(tab$dose)
  if (length(mark.doses) && max(mark.doses) > top) {
    ext <- dr_predicted(x$fit, sort(unique(c(xref, seq(top, max(mark.doses),
                                                         length.out = 25)))),
                        xref, x$provenance$level, x$sm,
                        x$provenance$single_study, x$shape,
                        range(x$frame$dose))$table
    ext <- ext[ext$dose > top, , drop = FALSE]
    tab <- rbind(tab, ext[, names(tab)])
  }
  unit <- if (is.null(x$dose_unit)) "dose" else x$dose_unit
  show_pi <- isTRUE(prediction) && isTRUE(x$prediction$estimable) &&
    !all(is.na(tab$pi.lb))
  if (isTRUE(prediction) && !show_pi) {
    warning("easydr: this fit has no prediction interval (",
            x$prediction$reason, "), so only the confidence band is drawn.",
            call. = FALSE)
  }

  if (is.null(xlab)) xlab <- paste0("Dose (", unit, ")")
  if (is.null(ylab)) ylab <- dr_curve_ylab(x$sm)
  if (is.null(main)) main <- paste0("Dose-response: ", x$shape$label)
  if (identical(main, FALSE)) main <- NULL
  # The studies behind the curve. A ratio measure already carries each category
  # as a log contrast against its study reference, so exponentiating it is the
  # reported effect. A continuous measure carries the group mean instead, and
  # the effect is that mean minus the study's own reference mean. The reference
  # categories themselves are left out: they sit on the null by construction and
  # would draw a row of points that says nothing. Each point keeps the dose its
  # study measured against, which is what its shape shows.
  obs <- NULL
  if (isTRUE(studies)) {
    obs <- x$frame
    is_ref   <- obs$is_ref
    ref      <- vapply(split(obs$y[is_ref], obs$studlab[is_ref]),
                       function(z) z[1L], numeric(1))
    ref_dose <- vapply(split(obs$dose[is_ref], obs$studlab[is_ref]),
                       function(z) z[1L], numeric(1))
    ref_n    <- vapply(split(obs$n[is_ref], obs$studlab[is_ref]),
                       function(z) z[1L], numeric(1))
    obs <- obs[!is_ref, , drop = FALSE]
    obs$est <- if (x$ratio) exp(obs$y) else obs$y - ref[obs$studlab]
    # A standardised mean difference puts each point on the SMD scale the curve
    # is drawn on, standardised by dosresmeta the same way the fit was.
    if (identical(x$sm, "SMD")) {
      fr  <- x$frame
      msd <- x$provenance$method.smd
      if (is.null(msd)) msd <- "cohens"
      smd <- unlist(unname(lapply(split(seq_len(nrow(fr)), fr$id), function(ix) {
        ix <- ix[order(!fr$is_ref[ix])]
        z  <- dosresmeta::covar.smd(y = fr$y[ix], sd = fr$cases[ix],
                                    n = fr$n[ix], measure = "smd",
                                    method = msd)$y
        stats::setNames(as.numeric(z), ix)
      })))
      keep <- which(!fr$is_ref)
      obs$est <- unname(smd[as.character(keep)])
    }
    obs$ref_dose <- unname(ref_dose[as.character(obs$studlab)])
    # A contrast rests on two groups, so its size is both: the category and
    # the reference it was measured against. With no counts at all, every
    # point gets the same size.
    obs$size_n <- obs$n + unname(ref_n[as.character(obs$studlab)])
    if (!any(is.finite(obs$size_n))) obs$size_n <- 1
    obs$size_n[!is.finite(obs$size_n)] <- stats::median(obs$size_n, na.rm = TRUE)
  } else if (isTRUE(studlab)) {
    warning("easydr: studlab needs studies = TRUE, so no labels were drawn.",
            call. = FALSE)
  }

  if (is.null(xlim)) xlim <- range(tab$dose)
  if (is.null(ylim)) {
    # The bubbles are data, so the axis makes room for them rather than
    # clipping a reported effect out of the figure.
    v    <- c(tab$ci.lb, tab$ci.ub, if (show_pi) c(tab$pi.lb, tab$pi.ub),
              obs$est)
    ylim <- dr_pad_range(range(c(v, if (x$ratio) 1 else 0), na.rm = TRUE),
                         x$ratio)
  }
  null_v <- if (x$ratio) 1 else 0

  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$dose))
  if (show_pi) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$pi.lb, ymax = .data$pi.ub), fill = fill.pi)
  }
  # Inside the observed doses the band is solid; past them it fades, and the
  # curve turns dashed, because there is no data there to support it. The
  # boundary point belongs to both pieces so they meet.
  obs_rng <- range(x$frame$dose)
  inside  <- tab$dose >= obs_rng[1L] - 1e-9 & tab$dose <= obs_rng[2L] + 1e-9
  tab_in  <- tab[inside, , drop = FALSE]
  edge    <- tab_in[which.max(tab_in$dose), , drop = FALSE]
  tab_out <- tab[!inside & tab$dose > obs_rng[2L], , drop = FALSE]
  if (nrow(tab_out)) tab_out <- rbind(edge, tab_out)
  # The band is light and slightly transparent, so the study points and the
  # knot lines behind it stay visible.
  p <- p + ggplot2::geom_ribbon(
    data = tab_in,
    ggplot2::aes(ymin = .data$ci.lb, ymax = .data$ci.ub), fill = fill.ci,
    alpha = 0.4)
  if (nrow(tab_out)) {
    p <- p + ggplot2::geom_ribbon(
      data = tab_out,
      ggplot2::aes(ymin = .data$ci.lb, ymax = .data$ci.ub), fill = fill.ci,
      alpha = 0.18)
  }

  # A spline's knots are where its pieces join, so they are drawn as lines
  # across the panel; the rug they replace only marked the reported doses.
  knots <- if (identical(x$shape$name, "rcs")) x$shape$knots else NULL
  show_knots <- isTRUE(knot.lines) && length(knots) > 0L
  if (is.null(rug)) rug <- !show_knots
  # A light dashed line at every labelled dose (the knots, the doses in the
  # table, the highest dose) lets a reader carry a value down to the axis.
  # The breaks are the ones the axis uses, worked out the same way below.
  if (show_knots || length(mark.doses)) {
    cand_v <- unique(signif(c(range(x$frame$dose), knots, mark.doses), 6))
    span_v <- diff(range(cand_v))
    gap_v  <- if (span_v > 0) span_v * 0.32 / max(axis_in - 1, 1) else 0
    br_v <- numeric(0)
    for (b in cand_v) if (all(abs(b - br_v) >= gap_v)) br_v <- c(br_v, b)
    br_v <- c(br_v[br_v > min(x$frame$dose)], if (show_knots) knots)
    p <- p + ggplot2::geom_vline(xintercept = unique(br_v), linetype = "dashed",
                                 colour = "grey78", linewidth = 0.3)
  }

  if (isTRUE(reference)) {
    p <- p + ggplot2::geom_hline(yintercept = null_v, linetype = "dashed",
                                 colour = "grey50", linewidth = 0.35)
    # The reference dose is marked only when it sits inside the panel: at the
    # lowest dose it would only double the axis edge.
    if (x$reference$curve > xlim[1L] + 1e-9) {
      p <- p + ggplot2::geom_vline(xintercept = x$reference$curve,
                                   linetype = "dotted", colour = "grey65",
                                   linewidth = 0.3)
    }
  }
  # Shapes for the reference doses, in the order a reader names them. A
  # subgroup figure already spends shape on its levels, and a seventh
  # reference would have to reuse one, so both fall back to circles.
  # The five fillable shapes come first, so they share one fill and one
  # outline; the star, which has no fill, is the sixth.
  ref_shapes <- c(21L, 22L, 24L, 23L, 25L, 8L)
  key_title <- paste0("Study reference dose (", unit, ")")
  ref_levels <- if (is.null(obs)) numeric(0) else sort(unique(obs$ref_dose))
  # One reference dose for every study tells the reader nothing a key could
  # add, so the key appears only when the references differ.
  shape_key  <- isTRUE(ref.shape) && !is.null(obs) &&
    length(ref_levels) > 1L && length(ref_levels) <= length(ref_shapes)
  if (isTRUE(ref.shape) && !is.null(obs) &&
      length(ref_levels) > length(ref_shapes)) {
    warning("easydr: the studies use ", length(ref_levels), " reference ",
            "doses and there are ", length(ref_shapes), " shapes, so every ",
            "study point is drawn as a circle.", call. = FALSE)
  }
  if (!is.null(obs)) {
    pt_fill <- grDevices::adjustcolor(col, alpha.f = 0.18)
    pt_line <- grDevices::adjustcolor(col, alpha.f = 0.55)
    # By subgroup, each study's circles take the colour of its level's curve,
    # so a reader can tell which points a curve was fitted to. A study whose
    # level could not be fitted keeps the overall colour.
    if (use_sub) {
      drawn_lv <- x$subgroup$levels[
        vapply(x$subgroup$fits, Negate(is.null), logical(1))]
      lv_col <- dr_subgroup_style(drawn_lv)$colour
      pt_lv  <- unname(x$subgroup$by_study[as.character(obs$studlab)])
      base   <- ifelse(pt_lv %in% drawn_lv, lv_col[pt_lv], col)
      pt_fill <- grDevices::adjustcolor(base, alpha.f = 0.18)
      pt_line <- grDevices::adjustcolor(base, alpha.f = 0.55)
    }
    if (shape_key) {
      labs_ref <- dr_num_lab(ref_levels)
      obs$ref_lab <- factor(dr_num_lab(obs$ref_dose), levels = labs_ref)
      p <- p +
        ggplot2::geom_point(
          data = obs, inherit.aes = FALSE,
          ggplot2::aes(x = .data$dose, y = .data$est, shape = .data$ref_lab,
                       size = .data$size_n),
          colour = pt_line, fill = pt_fill, stroke = 0.5) +
        ggplot2::scale_shape_manual(
          values = stats::setNames(ref_shapes[seq_along(labs_ref)], labs_ref),
          name = key_title, drop = FALSE,
          guide = ggplot2::guide_legend(
            title.position = "left", nrow = 1L, order = 2L,
            override.aes = list(size = 2.8, stroke = 0.5,
                                colour = grDevices::adjustcolor(col, alpha.f = 0.55),
                                fill = grDevices::adjustcolor(col, alpha.f = 0.18))))
    } else {
      p <- p + ggplot2::geom_point(
        data = obs, inherit.aes = FALSE,
        ggplot2::aes(x = .data$dose, y = .data$est, size = .data$size_n),
        shape = 21L, colour = pt_line, fill = pt_fill, stroke = 0.5)
    }
    p <- p + ggplot2::scale_size(range = c(1, 7), guide = "none")
    if (isTRUE(studlab)) {
      # One label per study, at its highest dose. A label on the right half of
      # the panel is written to the left of its bubble, so the last study does
      # not have its name cut off by the panel edge, and consecutive labels
      # alternate above and below, which is what keeps two studies that end at
      # the same dose from printing on top of each other.
      lab <- do.call(rbind, lapply(split(obs, obs$studlab), function(d)
        d[which.max(d$dose), , drop = FALSE]))
      mid <- mean(range(obs$dose))
      lab$hj <- ifelse(lab$dose > mid, 1.1, -0.1)
      lab$vj <- ifelse(seq_len(nrow(lab)) %% 2L == 0L, 1.6, -0.8)
      p <- p + ggplot2::geom_text(
        data = lab, inherit.aes = FALSE,
        ggplot2::aes(x = .data$dose, y = .data$est, label = .data$studlab,
                     hjust = .data$hj, vjust = .data$vj),
        size = 2.6, colour = "grey25")
    }
  }

  p <- p + ggplot2::geom_line(data = tab_in, ggplot2::aes(y = .data$est),
                              colour = col, linewidth = 1.3,
                              lineend = "round")
  if (nrow(tab_out)) {
    p <- p + ggplot2::geom_line(data = tab_out, ggplot2::aes(y = .data$est),
                                colour = col, linewidth = 1.1,
                                linetype = "dashed")
  }

  # One curve per subgroup level, over the overall band. Each level is told
  # apart twice over, by colour and by line type, because colour alone is lost
  # in print and to a colour-blind reader. The band stays the overall one: a
  # ribbon per level turns the panel to mud at three levels and hides the
  # curves it is meant to qualify.
  if (use_sub) {
    drawn <- x$subgroup$levels[
      vapply(x$subgroup$fits, Negate(is.null), logical(1))]
    sty <- dr_subgroup_style(drawn)
    sub_curve <- do.call(rbind, lapply(drawn, function(lv) {
      d <- x$subgroup$fits[[lv]]$curve
      data.frame(dose = d$dose, est = d$est, level = lv,
                 stringsAsFactors = FALSE)
    }))
    sub_marks <- do.call(rbind, lapply(drawn, function(lv) {
      d <- x$subgroup$fits[[lv]]$predicted
      d <- d[!d$is_ref, , drop = FALSE]
      if (nrow(d) == 0L) return(NULL)
      data.frame(dose = d$dose, est = d$est, level = lv,
                 stringsAsFactors = FALSE)
    }))
    sub_curve$level <- factor(sub_curve$level, levels = drawn)
    p <- p + ggplot2::geom_line(
      data = sub_curve, inherit.aes = FALSE,
      ggplot2::aes(x = .data$dose, y = .data$est, colour = .data$level,
                   linetype = .data$level), linewidth = 0.75)
    if (!is.null(sub_marks) && isTRUE(points)) {
      sub_marks$level <- factor(sub_marks$level, levels = drawn)
      p <- p + ggplot2::geom_point(
        data = sub_marks, inherit.aes = FALSE,
        ggplot2::aes(x = .data$dose, y = .data$est, colour = .data$level),
        shape = 21L, fill = "white", size = 2.2, stroke = 0.8)
    }
    p <- p +
      ggplot2::scale_colour_manual(values = sty$colour, name = x$subgroup$var,
                                   drop = FALSE,
                                   guide = ggplot2::guide_legend(order = 1L)) +
      ggplot2::scale_linetype_manual(values = sty$linetype,
                                     name = x$subgroup$var, drop = FALSE,
                                     guide = ggplot2::guide_legend(order = 1L))
  }

  if (isTRUE(rug)) {
    # The doses the fit reports, which are the ones `at` named, plus the
    # reference. A tick per exposure category put forty-eight of them under an
    # eight-study figure and read as a band rather than as marks.
    obs <- data.frame(dose = sort(unique(x$predicted$dose)))
    p <- p + ggplot2::geom_rug(data = obs, inherit.aes = FALSE,
                               ggplot2::aes(x = .data$dose), sides = "b",
                               colour = "grey55", linewidth = 0.4,
                               length = grid::unit(0.025, "npc"))
  }
  marks <- if (length(mark.doses)) {
    m <- dr_predicted(x$fit, sort(unique(c(xref, mark.doses))), xref,
                      x$provenance$level, x$sm, x$provenance$single_study,
                      x$shape, range(x$frame$dose))$table
    m[!m$is_ref, , drop = FALSE]
  } else {
    x$predicted[0L, , drop = FALSE]
  }
  if (isTRUE(points) && nrow(marks) > 0L) {
    p <- p + ggplot2::geom_point(data = marks, inherit.aes = FALSE,
                                 ggplot2::aes(x = .data$dose, y = .data$est),
                                 shape = 21L, fill = "white", colour = col,
                                 size = 2.2, stroke = 0.8)
    # With levels drawn, a number beside the overall curve sits on top of a
    # level's marker and reads as that level's value. The table under the
    # figure carries every number anyway.
    if (!is.null(round) && !use_sub) {
      marks$lab <- format(base::round(marks$est, round), nsmall = round)
      p <- p + ggplot2::geom_text(
        data = marks, inherit.aes = FALSE,
        ggplot2::aes(x = .data$dose, y = .data$est, label = .data$lab),
        vjust = -1.1, size = 3, colour = col)
    }
  }

  p <- p +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(
                     colour = "grey60", fill = NA, linewidth = 0.4),
                   axis.ticks = ggplot2::element_line(colour = "grey60",
                                                      linewidth = 0.3),
                   axis.text = ggplot2::element_text(colour = "grey20",
                                                     size = 9),
                   axis.title = ggplot2::element_text(colour = "grey10",
                                                      size = 10),
                   panel.grid.major = if (isTRUE(grid))
                     ggplot2::element_line(colour = "grey94",
                                           linewidth = 0.3) else
                       ggplot2::element_blank(),
                   # The keys run in one row above the panel, just under
                   # the tests, so they take no width from it.
                   legend.position = if (use_sub || shape_key) "top" else "none",
                   legend.direction = "horizontal",
                   legend.box = "horizontal",
                   legend.box.just = "left",
                   legend.box.spacing = grid::unit(6, "pt"),
                   legend.spacing.x = grid::unit(18, "pt"),
                   legend.justification = "left",
                   legend.location = "plot",
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.key.width = grid::unit(if (use_sub) 22 else 12, "pt"),
                   legend.key.height = grid::unit(if (use_sub) 10 else 12, "pt"),
                   legend.key.spacing.x = grid::unit(if (use_sub) 5.5 else 8,
                                                     "pt"),
                   legend.text = ggplot2::element_text(
                     size = 8, colour = "grey20",
                     margin = ggplot2::margin(l = 2)),
                   legend.title = ggplot2::element_text(
                     size = 8, colour = "grey20", hjust = 0, vjust = 0.5,
                     face = "bold", lineheight = 1.1,
                     margin = ggplot2::margin(r = 8)),
                   legend.key.spacing.y = grid::unit(3, "pt"),
                   plot.margin = ggplot2::margin(6, 10, 2, 6))
  p <- if (x$ratio) {
    # Plain numbers on the log axis: 1000, not 1e+03.
    p + ggplot2::scale_y_log10(limits = ylim,
                               labels = function(b) dr_num_lab(signif(b, 3)))
  } else {
    p + ggplot2::scale_y_continuous(limits = ylim)
  }
  # The dose axis is labelled where the figure has something to say: the
  # lowest and highest observed doses, the knots, and the doses reported in
  # the table.
  # Labels closer than about a third of an inch run into each other, so the
  # breaks are kept in that order of importance and a later one is dropped
  # when it would sit too close to one already kept. Its knot line still
  # shows where it is.
  cand <- unique(signif(c(range(x$frame$dose), knots, mark.doses), 6))
  span <- diff(range(cand))
  gap  <- if (span > 0) span * 0.32 / max(axis_in - 1, 1) else 0
  br <- numeric(0)
  for (b in cand) if (all(abs(b - br) >= gap)) br <- c(br, b)
  br <- sort(br)
  if (length(br) > 1L) {
    p <- p + ggplot2::scale_x_continuous(breaks = br,
                                         labels = dr_num_lab(signif(br, 3)))
  }
  p <- p + ggplot2::coord_cartesian(xlim = xlim)

  lvl <- format(100 * x$provenance$level)
  bands <- if (show_pi) {
    paste0("Shaded bands: ", lvl, "% confidence interval (darker) and ", lvl,
           "% prediction interval (lighter).")
  } else {
    paste0("Shaded band: ", lvl, "% confidence interval.")
  }
  # Under the figure only what explains its marks; the result and the tests
  # are in the subtitle.
  extra <- tests[setdiff(names(tests), c("wald", "deviance", "r2adj"))]
  # Band and symbols read as one legend, so they share a row.
  knot_txt <- if (show_knots) {
    kl <- dr_num_lab(signif(knots, 3))
    kl <- if (length(kl) > 1L) {
      paste0(paste(kl[-length(kl)], collapse = ", "), " and ", kl[length(kl)])
    } else kl
    paste0(" Knots at ", kl, " ", unit, ".")
  }
  foot <- c(paste0(bands, knot_txt,
                   if (isTRUE(studies))
                     " Symbols: individual study estimates relative to their own reference dose.",
                   if (use_sub) " The dark line and the band are the overall fit.",
                   if (nrow(tab_out))
                     " Dashed: beyond the highest observed dose."),
            if (length(extra))
              do.call(dr_test_notes, c(list(x, wald = FALSE), extra)),
            legend)

  # The target doses sit in the corner the curve leaves empty: bottom right
  # when it rises, top right when it falls.
  if (isTRUE(target) && !is.null(x$target) && !isTRUE(use_sub)) {
    tg <- x$target
    lab <- paste0(c(names(tg$ed), "ED max"), "  ",
                  # Four significant figures, as the shape line under the
                  # figure prints the ED50, so the two read the same.
                  formatC(c(tg$ed, tg$ed_max), format = "fg", digits = 4),
                  " ", unit, collapse = "\n")
    ref <- if (isTRUE(x$ratio)) 1 else 0
    rising <- utils::tail(tab$est, 1) >= ref
    p <- p + ggplot2::annotate(
      "label", x = Inf, y = if (rising) -Inf else Inf, label = lab,
      hjust = 1.04, vjust = if (rising) -0.12 else 1.12,
      size = 3.3, family = "mono", fill = "white", colour = "#1F2A44",
      label.r = grid::unit(0, "pt"),
      label.padding = grid::unit(5, "pt"))
  }

  list(plot = p, foot = foot, tab = tab, use_sub = use_sub,
       main = main, side_info = shape_key)
}


#' The linearity test, as one clause of a figure caption
#'
#' Every curve with a non-linear term carries the Wald test of those terms
#' against the straight line, so a reader of the figure alone sees whether
#' the bend is supported. A linear shape has nothing to test and says so.
#'
#' @keywords internal
#' @noRd
dr_linearity_note <- function(x, lead = ".  ") {
  nl <- x$nonlinearity
  if (identical(x$shape$name, "linear")) {
    return(paste0(lead, "Non-linearity test: not applicable to a linear shape"))
  }
  if (is.null(nl) || !isTRUE(nl$estimable)) {
    return(paste0(lead, "Non-linearity test: not available for this shape"))
  }
  paste0(lead, "Non-linearity test: \u03c7\u00b2 = ", formatC(nl$statistic, format = "f", digits = 2),
         " (df = ", nl$df, "), ", dr_p_lab(nl$p))
}


#' A dose written the way a reader would: no trailing zeros, no padding
#'
#' @keywords internal
#' @noRd
dr_num_lab <- function(v) {
  vapply(v, function(z) format(z, trim = TRUE, drop0trailing = TRUE),
         character(1), USE.NAMES = FALSE)
}


#' Keep the panel centred when a key sits to its right
#'
#' A key on the right takes its width from the panel, which then sits left of
#' the summary table under it. The same width is added as empty space on the
#' left, so the panel is centred on the figure. Inside the key, the shapes and
#' their values are centred under the title rather than left under it.
#'
#' @keywords internal
#' @noRd
dr_centre_legend <- function(p, centre = TRUE) {
  g  <- ggplot2::ggplotGrob(p)
  side <- function(g, idx) tryCatch(
    sum(dr_measure(grid::convertWidth(g$widths[idx], "cm", valueOnly = TRUE))),
    error = function(e) NA_real_)
  margins <- function(g) {
    pc <- unique(g$layout$l[g$layout$name == "panel"])[1L]
    n  <- length(g$widths)
    if (is.na(pc)) return(c(left = NA_real_, right = NA_real_))
    c(left = side(g, seq_len(pc - 1L)), right = side(g, (pc + 1L):n))
  }
  at <- which(g$layout$name == "guide-box-right")
  if (!isTRUE(centre) || !length(at) || inherits(g$grobs[[at]], "zeroGrob")) {
    attr(g, "panel_margins") <- margins(g)
    return(g)
  }
  box <- g$grobs[[at]]
  gi  <- which(box$layout$name == "guides")
  if (length(gi) == 1L) {
    gd <- box$grobs[[gi]]
    lab <- grepl("^label", gd$layout$name)
    n   <- length(gd$widths)
    lc  <- max(gd$layout$r[lab])
    if (lc + 1L <= n - 1L) {
      cm <- tryCatch(grid::convertWidth(gd$widths[(lc + 1L):(n - 1L)], "cm",
                                        valueOnly = TRUE),
                     error = function(e) NULL)
      if (!is.null(cm) && sum(cm) > 0) {
        half <- grid::unit(sum(cm) / 2, "cm")
        gd$widths[2L] <- gd$widths[2L] + half
        gd$widths[(lc + 1L):(n - 1L)] <- grid::unit(0, "cm")
        gd$widths[n - 1L] <- half
        ti <- which(gd$layout$name == "title")
        gd$layout$l[ti] <- 2L
        gd$layout$r[ti] <- n - 1L
        box$grobs[[gi]] <- gd
        g$grobs[[at]] <- box
      }
    }
  }
  # The panel box, not the whole plot, is what the eye centres on: pad the
  # left by what the right side (key, spacing, margin) has beyond the left
  # side (axis title, tick labels, margin).
  m   <- margins(g)
  pad <- m[["right"]] - m[["left"]]
  if (is.finite(pad) && pad > 0) {
    g <- gtable::gtable_add_cols(g, grid::unit(pad, "cm"), pos = 0L)
  }
  attr(g, "panel_margins") <- margins(g)
  g
}


#' Wrap subtitle groups, each with its own type size and colours
#'
#' A group is one line of the subtitle (what the curve rests on, the result,
#' the tests); it wraps only when wider than the panel. Each line keeps the
#' style of its group and the space above it.
#'
#' @keywords internal
#' @noRd
dr_rich_lines <- function(groups, styles, width_cm) {
  old <- grDevices::dev.cur()
  tf  <- tempfile(fileext = ".png")
  grDevices::png(tf, width = 7, height = 7, units = "in", res = 300)
  on.exit({
    grDevices::dev.off()
    if (old > 1L) grDevices::dev.set(old)
    unlink(tf)
  }, add = TRUE)
  out <- list()
  for (i in seq_along(groups)) {
    st <- styles[[i]]
    lines <- dr_subtitle_wrap1(groups[[i]], width_cm, size = st$size)
    for (j in seq_along(lines)) {
      runs <- lines[[j]]
      # A dot left at either end of a wrapped line parts nothing.
      while (length(runs) && identical(runs[[1L]]$role, "sep")) runs <- runs[-1L]
      while (length(runs) && identical(runs[[length(runs)]]$role, "sep")) {
        runs <- runs[-length(runs)]
      }
      if (length(runs)) runs[[1L]]$text <- sub("^ ", "", runs[[1L]]$text)
      out[[length(out) + 1L]] <- c(st, list(runs = runs,
                                            gap = if (j == 1L) st$gap else 0))
    }
  }
  out
}

#' Draw styled subtitle lines, centred, from plain, bold and separator pieces
#'
#' Each line is one plotmath expression drawn once per kind of piece, the
#' other pieces held as invisible phantoms, so the pieces meet exactly while
#' each kind keeps its own colour.
#'
#' @keywords internal
#' @noRd
dr_rich_grob <- function(lines, align = 0.5) {
  kind <- function(r) {
    if (identical(r$role, "sep")) "sep" else if (isTRUE(r$strong)) "strong" else "plain"
  }
  piece <- function(r, show) {
    e <- if (isTRUE(r$strong)) bquote(bold(.(r$text))) else r$text
    if (show) e else bquote(phantom(.(e)))
  }
  kids <- list()
  y <- 0
  for (ln in lines) {
    lh <- 1.4 * ln$size
    y  <- y + ln$gap
    yy <- grid::unit(1, "npc") - grid::unit(y + 0.5 * lh, "points")
    for (k in c("plain", "strong", "sep")) {
      if (!any(vapply(ln$runs, kind, character(1)) == k)) next
      parts <- lapply(ln$runs, function(r) piece(r, kind(r) == k))
      kids[[length(kids) + 1L]] <- grid::textGrob(
        as.expression(as.call(c(as.name("paste"), parts))),
        x = grid::unit(align, "npc"), y = yy, hjust = align,
        gp = grid::gpar(fontsize = ln$size,
                        col = switch(k, plain = ln$col, strong = ln$strong_col,
                                     sep = ln$sep_col)))
    }
    y <- y + lh
  }
  g <- grid::gTree(children = do.call(grid::gList, kids))
  attr(g, "height") <- grid::unit(y + 6, "points")
  g
}

#' The results table under a curve: a title row, then the table across the
#' figure
#'
#' The title names what the rows are and the note on the right what they are
#' compared with. The header row is shaded, the first column (the dose) is
#' bold, and the numbers are right-aligned under their headers.
#'
#' @keywords internal
#' @noRd
dr_curve_table_grob <- function(dt, title, note, base_size = 8.5,
                                margin = grid::unit(0.22, "in")) {
  dt[] <- lapply(dt, function(v) ifelse(is.na(v), "", as.character(v)))
  nr <- nrow(dt)
  nc <- ncol(dt)
  meas <- function(s, face) {
    if (!length(s)) return(0)
    max(vapply(s, function(z) dr_measure(grid::convertWidth(grid::grobWidth(
      grid::textGrob(z, gp = grid::gpar(fontsize = base_size, fontface = face))),
      "in", valueOnly = TRUE)), numeric(1)))
  }
  nat <- vapply(seq_len(nc), function(j) {
    max(meas(names(dt)[j], "plain"), meas(dt[[j]], if (j == 1L) "bold" else "plain"))
  }, numeric(1))
  rel <- nat + mean(nat)
  row_h <- 2.3 * base_size
  hdr_h <- 2.1 * base_size
  ttl_h <- 2.6 * base_size
  gt <- gtable::gtable(
    widths = grid::unit.c(margin, grid::unit(rel, "null"), margin),
    heights = grid::unit(c(ttl_h, hdr_h, rep(row_h, nr)), "points"))
  pad <- grid::unit(7, "pt")
  cell <- function(lab, j, face, col, size = base_size) {
    left <- j == 1L
    grid::textGrob(lab, x = if (left) pad else grid::unit(1, "npc") - pad,
                   hjust = if (left) 0 else 1,
                   gp = grid::gpar(fontsize = size, fontface = face, col = col))
  }
  line <- function(y, col) grid::segmentsGrob(
    x0 = grid::unit(0, "npc"), x1 = grid::unit(1, "npc"),
    y0 = grid::unit(y, "npc"), y1 = grid::unit(y, "npc"),
    gp = grid::gpar(col = col, lwd = 0.6))
  # Title row: the name on the left, what the rows are compared with on the
  # right, both on the table's edges.
  gt <- gtable::gtable_add_grob(gt, grid::textGrob(
    title, x = grid::unit(0, "npc"), hjust = 0,
    gp = grid::gpar(fontsize = base_size + 1.5, fontface = "bold", col = "grey10")),
    t = 1, l = 2, r = nc + 1, name = "title")
  if (length(note) && nzchar(note)) {
    gt <- gtable::gtable_add_grob(gt, grid::textGrob(
      note, x = grid::unit(1, "npc"), hjust = 1,
      gp = grid::gpar(fontsize = base_size - 0.5, col = "grey45")),
      t = 1, l = 2, r = nc + 1, name = "note")
  }
  gt <- gtable::gtable_add_grob(gt, grid::rectGrob(
    gp = grid::gpar(fill = "grey95", col = NA)), t = 2, l = 2, r = nc + 1,
    name = "header-fill")
  for (j in seq_len(nc)) {
    gt <- gtable::gtable_add_grob(gt, cell(names(dt)[j], j, "plain", "grey30"),
                                  t = 2, l = j + 1, name = paste0("hdr-", j))
  }
  for (i in seq_len(nr)) {
    gt <- gtable::gtable_add_grob(gt, line(0, "grey88"), t = i + 2, l = 2,
                                  r = nc + 1, name = paste0("rule-", i))
    for (j in seq_len(nc)) {
      gt <- gtable::gtable_add_grob(
        gt, cell(dt[[j]][i], j, if (j == 1L) "bold" else "plain", "grey10"),
        t = i + 2, l = j + 1, name = paste0("cell-", i, "-", j))
    }
  }
  attr(gt, "natural_in") <- sum(rel) / sum(rel) * (sum(nat) + nc * 0.4) +
    2 * grid::convertWidth(margin, "in", valueOnly = TRUE)
  gt
}


#' The table under a curve drawn by subgroup: one column per level and the
#' overall fit, one row per result
#'
#' The first row is the effect per increment (or the dose contrast asked
#' for) with its interval and p value; then heterogeneity (I2 and the Q
#' test), the deviance test and the adjusted R2.
#'
#' @keywords internal
#' @noRd
dr_subgroup_summary <- function(x, at = NULL, digits = 3) {
  fits <- c(x$subgroup$fits, list(Overall = x))
  fits <- fits[!vapply(fits, is.null, logical(1))]
  back <- if (x$ratio) exp else identity
  fmt  <- function(v) format(round(back(v), digits), nsmall = digits)
  hl0  <- dr_curve_headline(x, at, digits)
  eff <- function(f) {
    h <- dr_curve_headline(f, at, digits)
    q <- h$pooled
    lo <- fmt(q$lb)
    hi <- fmt(q$ub)
    ci <- if (grepl("^-", lo) || grepl("^-", hi)) paste0(lo, " to ", hi) else
      paste0(lo, "\u2013", hi)
    paste0(fmt(q$est), " (", ci, "; ", dr_p_lab(q$p), ")")
  }
  het <- function(f) {
    h <- f$heterogeneity
    if (isTRUE(h$estimable)) {
      sprintf("I\u00b2 = %.1f%%; Q test %s", h$i2, dr_p_lab(h$p))
    } else "I\u00b2 not estimable"
  }
  dev <- function(f) {
    g <- f$gof
    if (!is.null(g) && !identical(g$available, FALSE) && !is.null(g$p) &&
        is.finite(g$p)) {
      dr_p_lab(g$p)
    } else "not available"
  }
  r2 <- function(f) {
    g <- f$gof
    if (!is.null(g) && !identical(g$available, FALSE) && !is.null(g$r2adj) &&
        is.finite(g$r2adj)) {
      formatC(g$r2adj, format = "f", digits = 2)
    } else "not available"
  }
  out <- data.frame(` ` = c(sub(" [(]linear trend[)]$", "", hl0$label),
                            "Heterogeneity", "Deviance", "Adjusted R\u00b2"),
                    check.names = FALSE, stringsAsFactors = FALSE)
  for (lv in names(fits)) {
    f <- fits[[lv]]
    out[[lv]] <- c(eff(f), het(f), dev(f), r2(f))
  }
  out
}
