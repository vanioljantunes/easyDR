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
#' fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
#'               dose.unit = "g/day", at = c(10, 25, 50))
#'
#' # The defaults: confidence band, reported doses marked, null line, and a
#' # rug of the exposure categories the studies actually reported.
#' easydr_curve(fit)
#'
#' # Titled, labelled to two decimals, axis pinned, prediction band shown.
#' easydr_curve(fit,
#'              main  = "Alcohol and colorectal cancer",
#'              xlab  = "Alcohol intake (g/day)",
#'              ylab  = "Rate ratio",
#'              ylim  = c(0.9, 1.8),
#'              round = 2,
#'              prediction = TRUE)
#'
#' @export
easydr_curve <- function(x, main = NULL, xlab = NULL, ylab = NULL,
                         xlim = NULL, ylim = NULL, round = NULL,
                         prediction = FALSE, points = TRUE, reference = TRUE,
                         rug = TRUE, result = TRUE, result.at = NULL,
                         subgroup = NULL, digits = 3, col = "#1F2A44",
                         fill.ci = "grey75", fill.pi = "grey90",
                         legend = NULL) {
  if (!inherits(x, "easydr_result")) {
    stop("easydr: easydr_curve() takes a fitted result from easydr().",
         call. = FALSE)
  }
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
  if (is.null(xlim)) xlim <- range(tab$dose)
  if (is.null(ylim)) {
    v    <- c(tab$ci.lb, tab$ci.ub, if (show_pi) c(tab$pi.lb, tab$pi.ub))
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
                   if (use_sub)
                     ".  The dark line and the band are the overall fit.",
                   if (any(tab$extrapolated))
                     ".  Crosses mark doses beyond the observed range."),
            legend)

  # The summary sits under the figure, below the footer, the way the funnel
  # puts Egger's numbers under its panel. Inside the panel it covered the
  # curve on a fit that rises into its own corner, and it could hold only one
  # column, which left nowhere for a subgroup to go.
  if (isTRUE(result)) {
    tg <- dr_stats_grob(dr_summary_table(x, result.at, digits, use_sub))
    fs <- dr_footer_stack(foot, tg)
    dr_draw_stack(p, NULL, title = main, footer = fs$grob,
                  footer_h = fs$height, title_hjust = 0.5)
  } else {
    dr_draw_stack(p, NULL, title = main, footer = foot, title_hjust = 0.5)
  }
  invisible(tab)
}
