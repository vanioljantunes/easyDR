#' Plot the fitted dose-response curve
#'
#' @description
#' The pooled dose-response against dose, with its confidence band, optionally
#' its prediction band, and the doses the fit reports marked and labelled.
#' Every number was computed by [easydr()]; this chooses what appears and how.
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
#' @section Reading the axis:
#' A ratio measure is drawn on a log axis, because a halving and a doubling
#' should be the same distance from the null. The dashed horizontal is the
#' null, at 1 for a ratio and 0 for a difference, and the dotted vertical marks
#' the reference dose every effect is measured against.
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
#' @param rug `TRUE`, the default, ticks the observed dose categories along the
#'   bottom, so a gap in the data is visible under the curve drawn through it.
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
#' fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
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
                         rug = TRUE, col = "#1F2A44", fill.ci = "grey75",
                         fill.pi = "grey90", legend = NULL) {
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

  if (isTRUE(rug)) {
    obs <- data.frame(dose = sort(unique(x$frame$dose)))
    p <- p + ggplot2::geom_rug(data = obs, inherit.aes = FALSE,
                               ggplot2::aes(x = .data$dose), sides = "b",
                               colour = "grey35")
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
    if (!is.null(round)) {
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
                   if (any(tab$extrapolated))
                     ".  Crosses mark doses beyond the observed range."),
            legend)

  dr_draw_stack(p, NULL, title = main, footer = foot)
  invisible(tab)
}
