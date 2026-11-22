#' Funnel plot with a test for small-study effects
#'
#' @description
#' Each study's dose-response trend against its standard error, with the pooled
#' value, the pseudo confidence regions shaded the way `metafor::funnel()`
#' shades them, and Egger's regression test for asymmetry printed in a table
#' underneath.
#'
#' @section Why Egger and not Deeks:
#' Deeks' funnel plot asymmetry test belongs to diagnostic test accuracy,
#' where the effective sample size rather than the standard error is the
#' sensible horizontal axis. It does not apply here. A dose-response review
#' has one trend estimate per study with a standard error, which is the
#' setting Egger wrote his regression for: regress each study's standard
#' normal deviate on its precision and read the intercept.
#'
#' @section What asymmetry is and is not:
#' A funnel plot asks whether small studies sit away from the pooled value. It
#' cannot say why. Publication bias produces that pattern, and so do real
#' differences between small and large studies, poorer methods in the smaller
#' ones, and chance. Below about ten studies the test has little power and a
#' null result means almost nothing. The figure prints that warning.
#'
#' A curve has no single trend, so for a shape with more than one term the
#' comparison is the predicted effect at one dose, named on the axis.
#'
#' @param x a fitted result from [easydr()].
#' @param at what each study reports. `NULL` (the default) uses the
#'   `result.at` given to [easydr()]; one number (an increment of dose) or
#'   `c(numerator, denominator)` overrides it.
#' @param main plot title. `NULL` writes one; `FALSE` draws none.
#' @param table `"horizontal"` (the default) puts the statistics across the
#'   top of the table under the funnel, `"vertical"` puts them down its first
#'   column. Either way the table is headed "Egger's regression test".
#' @param xlab,ylab axis labels. `NULL` writes them from the fit.
#' @param levels the significance contours to shade, as percentages,
#'   outermost first.
#' @param shade fill colours for those regions, outermost first. The default
#'   greys get lighter outwards and stay pale, so the innermost region does
#'   not swamp the points.
#' @param labels `TRUE` writes the study labels beside the points.
#' @param digits digits for the table under the plot.
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
#' Viechtbauer W (2010). Conducting meta-analyses in R with the metafor
#' package. Journal of Statistical Software 36(3), 1-48.
#'
#' @examples
#' easydr.object <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
#'                         dose.unit = "g/day")
#' f <- easydr_funnel(easydr.object, labels = TRUE)
#' f$egger
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
#' @export
easydr_funnel <- function(x, at = NULL, main = NULL, xlab = NULL, ylab = NULL,
                          levels = c(99, 95, 90),
                          shade = c("grey96", "grey90", "grey83"),
                          labels = FALSE, digits = 4, table = "horizontal",
                          filetype = c("png", "pdf", "tiff"),
                          path = NULL, file.name = NULL,
                          width = NULL, height = NULL) {
  dr_require_single(x, "easydr_funnel()")
  e   <- dr_study_effects(x, at, what = "a funnel plot")
  tab <- e$table
  ok  <- !is.na(tab$est) & !is.na(tab$se) & tab$se > 0
  if (!any(ok)) {
    stop("easydr: no study has an estimate with a standard error, so there is ",
         "nothing to plot.", call. = FALSE)
  }
  if (length(levels) != length(shade)) {
    stop("easydr: `levels` and `shade` must have the same length; got ",
         length(levels), " and ", length(shade), ".", call. = FALSE)
  }
  eg <- dr_egger(tab$est, tab$se)
  # Classed so it prints as a test result. Unlisting and rounding it by hand
  # was the only way to read it before.
  class(eg) <- "easydr_egger"

  pts    <- tab[ok, , drop = FALSE]
  se_max <- max(pts$se) * 1.12
  ctr    <- e$pooled$est

  # The shaded wedges, widest first so the narrower ones draw over them. Each
  # is the region a study of that precision would fall in under the pooled
  # value, which is what metafor::funnel() shades.
  ord <- order(levels, decreasing = TRUE)
  wedges <- do.call(rbind, lapply(ord, function(i) {
    z <- stats::qnorm(1 - (1 - levels[i] / 100) / 2)
    data.frame(xx = c(ctr, ctr - z * se_max, ctr + z * se_max),
               yy = c(0, se_max, se_max),
               grp = factor(levels[i], levels = levels[ord]),
               stringsAsFactors = FALSE)
  }))
  fills <- stats::setNames(shade[ord], as.character(levels[ord]))
  # The legend reads from the narrowest contour outwards: 90%, 95%, 99%.
  brk <- as.character(sort(levels))

  if (is.null(main)) main <- "Funnel plot of the study trends"
  if (identical(main, FALSE)) main <- NULL
  if (is.null(xlab)) xlab <- e$label
  if (is.null(ylab)) ylab <- "Standard error"

  lo <- min(c(pts$est, ctr - stats::qnorm(0.995) * se_max))
  hi <- max(c(pts$est, ctr + stats::qnorm(0.995) * se_max))
  ticks <- dr_plot_ticks(lo, hi, e$ratio)
  keep  <- ticks$at >= lo & ticks$at <= hi

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = wedges,
                          ggplot2::aes(x = .data$xx, y = .data$yy,
                                       fill = .data$grp, group = .data$grp),
                          colour = "grey60", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = fills, breaks = brk,
                               labels = paste0(brk, "%"),
                               name = "Significance contours:") +
    ggplot2::geom_vline(xintercept = ctr, linetype = "dashed",
                        colour = "black", linewidth = 0.45) +
    ggplot2::geom_point(data = pts,
                        ggplot2::aes(x = .data$est, y = .data$se),
                        shape = 24L, fill = "white", colour = "black",
                        size = 2.1, stroke = 0.6)
  if (isTRUE(labels)) {
    # Above the marker, not beside it: a label beside a triangle sits on top
    # of it as soon as two studies have a similar precision.
    p <- p + ggplot2::geom_text(
      data = pts, ggplot2::aes(x = .data$est, y = .data$se,
                               label = .data$study),
      vjust = -0.9, size = 2.6, colour = "grey25")
  }
  p <- p +
    ggplot2::scale_x_continuous(breaks = ticks$at[keep],
                                labels = ticks$labels[keep]) +
    ggplot2::scale_y_reverse(limits = c(se_max, 0)) +
    ggplot2::coord_cartesian(xlim = c(lo, hi)) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(colour = "grey60", linewidth = 0.25))) +
    ggplot2::theme_bw(base_family = "sans") +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey93",
                                                            linewidth = 0.3),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top",
                   legend.direction = "horizontal",
                   legend.justification = "left",
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.box.spacing = grid::unit(4, "pt"),
                   legend.key.size = grid::unit(10, "pt"),
                   legend.text = ggplot2::element_text(
                     size = 8, margin = ggplot2::margin(l = 3, r = 8)),
                   legend.title = ggplot2::element_text(
                     size = 8, margin = ggplot2::margin(r = 6)),
                   plot.subtitle = ggplot2::element_text(size = 8,
                                                         lineheight = 1.3,
                                                         colour = "grey25"),
                   plot.margin = ggplot2::margin(6, 10, 2, 6))

  stats_df <- if (isTRUE(eg$estimable)) {
    data.frame(
      Statistic = c("Number of studies", "Egger intercept", "Standard error",
                    "t", "p value"),
      Value = c(as.character(eg$k),
                format(round(eg$intercept, digits), nsmall = digits),
                format(round(eg$se, digits), nsmall = digits),
                paste0(format(round(eg$t, 2), nsmall = 2), " on ", eg$df,
                       " df"),
                dr_p(eg$p)),
      stringsAsFactors = FALSE)
  } else {
    data.frame(Statistic = c("Number of studies", "Egger's test"),
               Value = c(as.character(eg$k), "not estimable"),
               stringsAsFactors = FALSE)
  }

  note <- character()
  if (isTRUE(eg$estimable) && isTRUE(eg$underpowered)) {
    note <- paste0("Fewer than ", eg$k_warn, " studies: the test has little ",
                   "power and a null result means very little.")
  } else if (!isTRUE(eg$estimable)) {
    note <- eg$reason
  }
  note <- c(note,
            "Asymmetry has several causes. Publication bias is only one.")

  # The caveats sit above the panel with the legend, and the table is the only
  # thing under it. Two stacked bands beneath one figure is one too many:
  # whichever grob measures short lets the other one land on top of it.
  p  <- p + ggplot2::labs(subtitle = paste(note, collapse = "\n"))
  tab_g <- dr_funnel_table(stats_df, table)
  tg    <- dr_funnel_footer(tab_g)

  # A funnel is read as a shape, so its panel stays close to square. The
  # subtitle carries the caveats and is drawn inside the panel, so it is
  # counted in: it squeezes the panel from the top.
  spec <- dr_file_spec(file.name, filetype, path, width, height,
                       what = "easydr_funnel()")
  title_in <- if (is.null(main)) 0 else 1.9 * 12 / 72
  # The subtitle and the shaded-region legend are both drawn above the panel
  # and both take their space from it. Counting only the subtitle pushed the
  # last row of the table off the bottom of the page.
  note_in  <- (1.35 * 11 * length(note) + 10) / 72 + 0.40
  w_in     <- max(7.0, dr_grob_width_in(tab_g) + 0.5)
  size <- dr_fig_size(w_in,
                      min(max(0.52 * w_in, 3.3), 4.6) + note_in + title_in +
                        dr_height_in(attr(tg, "height")) + 0.25)
  con  <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)

  g  <- dr_draw_stack(p, NULL, title = main, footer = tg,
                      footer_h = attr(tg, "height"))

  invisible(list(table = tab, pooled = e$pooled, egger = eg, label = e$label,
                 at = e$at, grob = g))
}

#' Egger's test as the table under the funnel
#'
#' The shared table style (bold labels over a thin rule, light row lines) with
#' the test's name as a plain bold caption in place of a grey header cell.
#'
#' @keywords internal
#' @noRd
dr_funnel_table <- function(stats_df, table = "horizontal", base_size = 9) {
  table <- match.arg(table, c("horizontal", "vertical"))
  if (identical(table, "horizontal")) {
    m <- as.data.frame(t(stats_df$Value), stringsAsFactors = FALSE)
    names(m) <- stats_df$Statistic
  } else {
    m <- stats_df
  }
  g <- dr_stats_grob(m, base_size)
  # The rule above the column labels would double the separator drawn across
  # the whole footer, so it goes.
  keep <- g$layout$name != "top-rule"
  g$grobs  <- g$grobs[keep]
  g$layout <- g$layout[keep, , drop = FALSE]
  g <- gtable::gtable_add_rows(g, grid::unit(2.1 * base_size, "points"),
                               pos = 0)
  gtable::gtable_add_grob(
    g, grid::textGrob("Egger\u2019s regression test", x = grid::unit(5, "pt"),
                      hjust = 0,
                      gp = grid::gpar(fontsize = base_size + 1,
                                      fontface = "bold")),
    t = 1, l = 1, r = ncol(g), name = "egger-title", clip = "off")
}

#' The footer band: a light rule across the figure, the table just under it
#' @keywords internal
#' @noRd
dr_funnel_footer <- function(tab_g, gap = 4, pad = 8) {
  th <- grid::convertHeight(sum(tab_g$heights), "points", valueOnly = TRUE)
  tab_g$vp <- grid::viewport(x = grid::unit(6, "pt"),
                             y = grid::unit(1, "npc") - grid::unit(gap, "pt"),
                             just = c("left", "top"), width = sum(tab_g$widths),
                             height = grid::unit(th, "points"))
  rule <- grid::segmentsGrob(
    x0 = grid::unit(6, "pt"), x1 = grid::unit(1, "npc") - grid::unit(10, "pt"),
    y0 = grid::unit(1, "npc"), y1 = grid::unit(1, "npc"),
    gp = grid::gpar(col = "grey75", lwd = 0.6))
  g <- grid::grobTree(rule, tab_g, name = "egger-footer")
  attr(g, "height") <- grid::unit(th + gap + pad, "points")
  g
}
