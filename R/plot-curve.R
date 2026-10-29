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
#'   summary measure.
#' @param xlim,ylim axis limits on the reported scale, for example
#'   `ylim = c(0.8, 2)`. `NULL` fits them to whatever is drawn.
#' @param round digits for the labels printed beside the marked doses. `NULL`
#'   marks the doses without labelling them.
#' @param prediction `TRUE` adds the prediction band when the fit has one.
#'   `FALSE`, the default, draws the confidence band alone.
#' @param points `TRUE`, the default, marks the doses the fit reports, which
#'   are the ones named in `at`.
#' @param studies `TRUE` draws the evidence behind the curve: one bubble per
#'   reported dose category, at the effect that study reported for it against
#'   its own reference, with the bubble area scaled to the participants that
#'   category carried. The reference categories are not drawn, because every
#'   one of them sits on the null by construction. `FALSE` is the default.
#' @param studlabel `TRUE` writes each study's label once, beside its
#'   highest-dose bubble, so a reader can tell which study a run of bubbles
#'   belongs to. It needs `studies = TRUE`.
#' @param reference `TRUE`, the default, draws the null line and marks the
#'   reference dose.
#' @param rug `TRUE`, the default, ticks the doses the fit reports along the
#'   bottom: the ones `at` named in [easydr()], and the reference. It marks
#'   where the reported effects sit on the axis, not where the data is.
#' @param result `TRUE`, the default, puts a summary table under the figure:
#'   the effect with its confidence interval, its p value, the heterogeneity
#'   and the number of studies. The figure then carries its own numbers. With
#'   subgroups drawn, the table gains one column per level and ends with the
#'   overall fit.
#' @param result.at the dose that table reports, for a shape with more than one
#'   term. `NULL` uses the first dose the fit reports. Ignored for a straight
#'   line, which has one slope.
#' @param subgroup `TRUE` draws one curve per subgroup level, `FALSE` draws the
#'   overall fit alone. The default, `NULL`, draws the levels when the fit
#'   carries them. Each level is told apart by colour and by line type, and its
#'   marked doses by shape. The band stays the overall one. `TRUE` on a fit
#'   with no subgroup is an error rather than a silent overall curve; see the
#'   `subgroup` argument of [easydr()].
#' @param digits digits for the numbers in that table.
#' @param col,fill.ci,fill.pi colours for the curve, the confidence band and
#'   the prediction band.
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
                         prediction = FALSE, points = TRUE, studies = FALSE,
                         studlabel = FALSE, reference = TRUE,
                         rug = TRUE, result = TRUE, result.at = NULL,
                         subgroup = NULL, digits = 3, col = "#1F2A44",
                         fill.ci = "grey75", fill.pi = "grey90",
                         legend = NULL,
                         filetype = c("png", "pdf", "tiff"),
                         path = NULL, file.name = NULL,
                         width = NULL, height = NULL,
                         mult.curves = FALSE,
                         curve.panels = NULL) {
  # Several shapes go to the panel drawer, which builds one of these plots
  # per shape and lays them out. A single fit asked for panels is a mistake
  # worth naming rather than quietly drawing one panel.
  if (inherits(x, "easydr_multi") || isTRUE(mult.curves)) {
    return(dr_curve_multi(
      x, mult.curves = mult.curves, curve.panels = curve.panels,
      main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
      round = round, prediction = prediction, points = points,
      studies = studies, studlabel = studlabel, reference = reference,
      rug = rug, result = result, result.at = result.at, digits = digits,
      col = col, fill.ci = fill.ci, fill.pi = fill.pi, legend = legend,
      filetype = filetype, path = path, file.name = file.name,
      width = width, height = height))
  }
  dr_require_single(x, "easydr_curve()")

  built   <- dr_curve_panel(
    x, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
    round = round, prediction = prediction, points = points,
    studies = studies, studlabel = studlabel, reference = reference,
    rug = rug, subgroup = subgroup, col = col, fill.ci = fill.ci,
    fill.pi = fill.pi, legend = legend)
  p       <- built$plot
  foot    <- built$foot
  tab     <- built$tab
  use_sub <- built$use_sub
  main    <- built$main

  # The summary sits under the figure, below the footer, the way the funnel
  # puts Egger's numbers under its panel. Inside the panel it covered the
  # curve on a fit that rises into its own corner, and it could hold only one
  # column, which left nowhere for a subgroup to go.
  if (isTRUE(result)) {
    tg <- dr_stats_grob(dr_summary_table(x, result.at, digits, use_sub))
    fs <- dr_footer_stack(foot, tg)
  } else {
    tg <- NULL
    fs <- NULL
  }

  # The panel is the figure here, and the band under it is whatever the
  # footer needs. The width follows the summary table when there is one,
  # because a subgroup table with four columns is wider than a plot of a
  # single curve wants to be.
  spec <- dr_file_spec(file.name, filetype, path, width, height,
                       what = "easydr_curve()")
  foot_in  <- if (is.null(fs)) {
    (1.35 * 9 * length(foot) + 8) / 72
  } else {
    dr_height_in(fs$height)
  }
  title_in <- if (is.null(main)) 0 else 1.9 * 12 / 72
  w_in     <- max(7.2, dr_grob_width_in(tg) + 0.5)
  size <- dr_fig_size(w_in,
                      min(max(0.50 * w_in, 3.1), 4.4) + foot_in + title_in +
                        0.12)
  con  <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)

  if (is.null(fs)) {
    dr_draw_stack(p, NULL, title = main, footer = foot, title_hjust = 0.5)
  } else {
    dr_draw_stack(p, NULL, title = main, footer = fs$grob,
                  footer_h = fs$height, title_hjust = 0.5)
  }
  invisible(tab)
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
                           prediction = FALSE, points = TRUE, studies = FALSE,
                           studlabel = FALSE, reference = TRUE, rug = TRUE,
                           subgroup = NULL, col = "#1F2A44",
                           fill.ci = "grey75", fill.pi = "grey90",
                           legend = NULL) {
  use_sub <- dr_use_subgroups(x, subgroup, "easydr_curve()")
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
  if (identical(main, FALSE)) main <- NULL
  # The studies behind the curve. A ratio measure already carries each category
  # as a log contrast against its study reference, so exponentiating it is the
  # reported effect. A continuous measure carries the group mean instead, and
  # the effect is that mean minus the study's own reference mean. The reference
  # categories themselves are left out: they sit on the null by construction and
  # would draw a row of bubbles that says nothing.
  obs <- NULL
  if (isTRUE(studies)) {
    obs <- x$frame
    ref <- vapply(split(obs$y[obs$is_ref], obs$studlab[obs$is_ref]),
                  function(z) z[1L], numeric(1))
    obs <- obs[!obs$is_ref, , drop = FALSE]
    obs$est <- if (x$ratio) exp(obs$y) else obs$y - ref[obs$studlab]
    if (!any(is.finite(obs$n))) obs$n <- 1
  } else if (isTRUE(studlabel)) {
    warning("easydr: studlabel needs studies = TRUE, so no labels were drawn.",
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
  p <- p + ggplot2::geom_ribbon(
    ggplot2::aes(ymin = .data$ci.lb, ymax = .data$ci.ub), fill = fill.ci)

  if (isTRUE(reference)) {
    p <- p +
      ggplot2::geom_hline(yintercept = null_v, linetype = "dashed",
                          colour = "grey35", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = x$reference$curve, linetype = "dotted",
                          colour = "grey45", linewidth = 0.4)
  }
  if (!is.null(obs)) {
    p <- p +
      ggplot2::geom_point(
        data = obs, inherit.aes = FALSE,
        ggplot2::aes(x = .data$dose, y = .data$est, size = .data$n),
        colour = col, fill = col, shape = 21L, alpha = 0.25, stroke = 0) +
      ggplot2::scale_size_area(max_size = 7, guide = "none")
    if (isTRUE(studlabel)) {
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

  p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$est), colour = col,
                              linewidth = 0.9)

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
        ggplot2::aes(x = .data$dose, y = .data$est, colour = .data$level,
                     shape = .data$level), fill = "white", size = 2.2,
        stroke = 0.8)
    }
    p <- p +
      ggplot2::scale_colour_manual(values = sty$colour, name = x$subgroup$var,
                                   drop = FALSE) +
      ggplot2::scale_linetype_manual(values = sty$linetype,
                                     name = x$subgroup$var, drop = FALSE) +
      ggplot2::scale_shape_manual(values = sty$shape, name = x$subgroup$var,
                                  drop = FALSE)
  }

  if (isTRUE(rug)) {
    # The doses the fit reports, which are the ones `at` named, plus the
    # reference. A tick per exposure category put forty-eight of them under an
    # eight-study figure and read as a band rather than as marks.
    obs <- data.frame(dose = sort(unique(x$predicted$dose)))
    p <- p + ggplot2::geom_rug(data = obs, inherit.aes = FALSE,
                               ggplot2::aes(x = .data$dose), sides = "b",
                               colour = "grey25", linewidth = 0.5,
                               length = grid::unit(0.035, "npc"))
  }
  if (any(tab$extrapolated)) {
    ex <- tab[tab$extrapolated, , drop = FALSE]
    p <- p + ggplot2::geom_point(data = ex, inherit.aes = FALSE,
                                 ggplot2::aes(x = .data$dose, y = .data$est),
                                 shape = 4L, colour = "#C0392B", size = 1.6)
  }

  marks <- x$predicted[!x$predicted$is_ref, , drop = FALSE]
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
                   legend.position = if (use_sub) "top" else "none",
                   legend.direction = "horizontal",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0, 0, 2, 0),
                   legend.key.width = grid::unit(22, "pt"),
                   legend.key.height = grid::unit(10, "pt"),
                   legend.text = ggplot2::element_text(size = 8),
                   legend.title = ggplot2::element_text(size = 8),
                   plot.margin = ggplot2::margin(6, 10, 2, 6))
  p <- if (x$ratio) {
    p + ggplot2::scale_y_log10(limits = ylim)
  } else {
    p + ggplot2::scale_y_continuous(limits = ylim)
  }
  p <- p + ggplot2::coord_cartesian(xlim = xlim)

  bands <- paste0(format(100 * x$provenance$level), "% CI")
  if (show_pi) {
    bands <- paste0(bands, " and ", format(100 * x$provenance$level), "% PI")
  }
  foot <- c(paste0("Reference ", x$reference$curve, " ", unit,
                   ".  Band: ", bands,
                   ".  Shape = ", x$shape$label,
                   if (isTRUE(studies))
                     ".  Bubbles: reported effects, sized by participants.",
                   if (use_sub)
                     ".  The dark line and the band are the overall fit.",
                   if (any(tab$extrapolated))
                     ".  Crosses mark doses beyond the observed range."),
            legend)

  list(plot = p, foot = foot, tab = tab, use_sub = use_sub,
       main = main)
}
