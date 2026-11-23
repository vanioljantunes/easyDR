# The picture a leave-one-out analysis should produce.
#
# easyDTA draws the summary curve once per omitted study, letters each one, and
# puts the numbers in a bordered table underneath so the figure reads on its
# own. The same picture answers the same question here: every refit is a whole
# dose-response curve, so the curves are what go on the plot, and a line that
# wanders away from the bold one names the study that was carrying the result.

#' Draw the leave-one-out curves and their table
#'
#' @keywords internal
#' @noRd
dr_loo_plot <- function(out, object, main = NULL, xlab = NULL, ylab = NULL,
                        ylim = NULL, digits = 3, position = "top",
                        filetype = c("png", "pdf", "tiff"),
                        path = NULL, file.name = NULL,
                        width = NULL, height = NULL) {
  # `out` carries the refits and the table; `object` is the fit they were all
  # taken from, which is what the bold curve and the final table row show.
  base  <- object
  tab   <- out$table
  fits  <- out$refits
  ratio <- isTRUE(out$ratio)
  unit  <- if (is.null(out$dose_unit)) "dose" else out$dose_unit

  ok <- tab$status == "fitted" & vapply(fits, function(f)
    inherits(f, "easydr_result"), logical(1))
  if (!any(ok)) {
    stop("easydr: no refit succeeded, so there is no curve to draw. The ",
         "printed table gives the reason against each row.", call. = FALSE)
  }

  # One letter per omitted study, in the order the table is already sorted by,
  # which is how far the result moved. A reader looks up the letter, not the
  # colour, so the figure survives being printed in black and white.
  keys <- if (sum(ok) <= 26L) LETTERS[seq_len(sum(ok))] else
    as.character(seq_len(sum(ok)))
  tab$key <- NA_character_
  tab$key[ok] <- keys

  curves <- do.call(rbind, lapply(which(ok), function(i) {
    cv <- fits[[i]]$curve
    data.frame(dose = cv$dose, est = cv$est, key = tab$key[i],
               stringsAsFactors = FALSE)
  }))
  pooled <- data.frame(dose = base$curve$dose, est = base$curve$est)

  if (is.null(ylab)) ylab <- dr_sm_axis(object$sm)
  if (is.null(xlab)) xlab <- paste0("Dose (", unit, ")")
  if (is.null(main)) main <- "Dose-response omitting each study"
  if (identical(main, FALSE)) main <- NULL
  if (is.null(ylim)) {
    ylim <- dr_pad_range(range(c(curves$est, pooled$est,
                                 if (ratio) 1 else 0), na.rm = TRUE), ratio)
  }

  # Each curve is labelled where it sits furthest from the pooled one, which
  # is the point a reader is looking at anyway.
  lab_at <- do.call(rbind, lapply(split(curves, curves$key), function(d) {
    gap <- abs(d$est - pooled$est[match(d$dose, pooled$dose)])
    d[which.max(gap), , drop = FALSE]
  }))

  p <- ggplot2::ggplot(curves, ggplot2::aes(x = .data$dose, y = .data$est)) +
    ggplot2::geom_hline(yintercept = if (ratio) 1 else 0, linetype = "dashed",
                        colour = "grey45", linewidth = 0.4) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$key), linewidth = 0.55) +
    ggplot2::geom_line(data = pooled, inherit.aes = FALSE,
                       ggplot2::aes(x = .data$dose, y = .data$est),
                       colour = "black", linewidth = 1.4) +
    ggplot2::geom_text(data = lab_at, inherit.aes = FALSE,
                       ggplot2::aes(x = .data$dose, y = .data$est,
                                    label = .data$key, colour = .data$key),
                       fontface = "bold", size = 3.1, vjust = -0.7,
                       show.legend = FALSE) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   plot.margin = ggplot2::margin(6, 10, 2, 6))
  p <- if (ratio) {
    p + ggplot2::scale_y_log10(limits = ylim)
  } else {
    p + ggplot2::scale_y_continuous(limits = ylim)
  }

  # The rows under the figure are drawn as a forest: one row per omission
  # with its interval on the effect axis, the fit with nothing omitted as a
  # diamond at the foot, and the statistics in columns either side.
  fmt <- function(v) ifelse(is.na(v), "-",
                            format(round(v, digits), nsmall = digits))
  pct <- function(v) ifelse(is.na(v), "-",
                            paste0(format(round(v, 1), nsmall = 1), "%"))
  r2f <- function(v) ifelse(is.na(v), "-", format(round(v, 2), nsmall = 2))
  shift <- ifelse(is.na(tab$shift), "-",
                  paste0(ifelse(tab$shift >= 0, "+", ""),
                         format(round(tab$shift, digits), nsmall = digits)))
  full_fit <- out$full$fit
  if (is.null(full_fit)) full_fit <- base
  full_i2 <- if (isTRUE(base$heterogeneity$estimable)) {
    pct(base$heterogeneity$i2)
  } else {
    "-"
  }

  n     <- nrow(tab)
  ypos  <- c(rev(seq_len(n)) + 1, 0.6)
  face  <- c(rep("plain", n), "bold")
  hdr_y <- n + 2
  ylim  <- c(0.2, hdr_y + 0.5)
  zebra <- dr_zebra(ypos[seq_len(n)])

  left_cols <- list(` ` = c(ifelse(is.na(tab$key), "", tab$key), ""),
                    Omitting = c(tab$label, "Overall"))
  right_cols <- list(
    x = c(ifelse(is.na(tab$lb), "-",
                 paste0(fmt(tab$estimate), " (", fmt(tab$lb), " to ",
                        fmt(tab$ub), ")")),
          paste0(fmt(out$full$est), " (", fmt(out$full$lb), " to ",
                 fmt(out$full$ub), ")")),
    p = c(dr_fmt_p(tab$p_effect), dr_fmt_p(if (is.null(out$full$p)) NA else out$full$p)))
  names(right_cols)[1L] <- paste0(object$sm, " (", format(100 * out$level), "% CI)")
  if (isTRUE(out$show_shift)) right_cols$Shift <- c(shift, "-")
  right_cols <- c(right_cols, list(
    "I\u00b2" = c(pct(tab$i2), full_i2),
    `p non-linearity` = c(dr_fmt_p(tab$p_nonlin),
                          dr_fmt_p(dr_sens_stat(full_fit$nonlinearity, "estimable", "p"))),
    `p deviance` = c(dr_fmt_p(tab$p_dev),
                     dr_fmt_p(dr_sens_stat(full_fit$gof, "available", "p"))),
    "Adjusted R\u00b2" = c(r2f(tab$r2adj),
                      r2f(dr_sens_stat(full_fit$gof, "available", "r2adj")))))

  p_left  <- dr_text_panel(left_cols, c(FALSE, FALSE), ypos = ypos,
                           face = face, header_y = hdr_y, ylim = ylim,
                           zebra = zebra)
  p_right <- dr_text_panel(right_cols, rep(TRUE, length(right_cols)),
                           ypos = ypos, face = face, header_y = hdr_y,
                           ylim = ylim, zebra = zebra,
                           center = rep(TRUE, length(right_cols)))

  tr  <- if (ratio) log else identity
  est <- tr(c(tab$estimate, out$full$est))
  lb  <- tr(c(tab$lb, out$full$lb))
  ub  <- tr(c(tab$ub, out$full$ub))
  null <- if (ratio) 0 else 0
  lo  <- min(c(lb, null), na.rm = TRUE); hi <- max(c(ub, null), na.rm = TRUE)
  pad <- 0.08 * (hi - lo); lo <- lo - pad; hi <- hi + pad
  rows <- data.frame(y = ypos[seq_len(n)], est = est[seq_len(n)],
                     lb = lb[seq_len(n)], ub = ub[seq_len(n)])
  rows <- rows[is.finite(rows$est), , drop = FALSE]
  dh  <- 0.32
  dia <- data.frame(x = c(lb[n + 1], est[n + 1], ub[n + 1], est[n + 1]),
                    y = 0.6 + c(0, dh, 0, -dh))
  ticks <- dr_plot_ticks(lo, hi, ratio)
  keep  <- ticks$at >= lo & ticks$at <= hi
  xlab_f <- out$label
  p_eff <- ggplot2::ggplot(rows, ggplot2::aes(x = .data$est, y = .data$y)) +
    zebra +
    ggplot2::annotate("segment", x = null, xend = null, y = 0.2,
                      yend = hdr_y - 0.5, colour = "grey25", linewidth = 0.4) +
    ggplot2::annotate("segment", x = est[n + 1], xend = est[n + 1], y = 0.2,
                      yend = hdr_y - 0.5, linetype = "dashed",
                      colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lb, xmax = .data$ub),
                           width = 0.26, orientation = "y", linewidth = 0.45) +
    ggplot2::geom_point(shape = 15L, size = 2.2) +
    ggplot2::geom_polygon(data = dia, inherit.aes = FALSE,
                          ggplot2::aes(x = .data$x, y = .data$y),
                          fill = "black") +
    ggplot2::scale_x_continuous(breaks = ticks$at[keep],
                                labels = ticks$labels[keep]) +
    ggplot2::coord_cartesian(xlim = c(lo, hi), ylim = ylim) +
    ggplot2::labs(x = xlab_f) +
    dr_gg_theme(axis = TRUE) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold", size = 9,
                                                        margin = ggplot2::margin(t = 4)))

  row_in  <- 0.24
  lw <- dr_panel_width(left_cols); rw <- dr_panel_width(right_cols)
  ew <- 2.6
  fg <- gridExtra::arrangeGrob(
    ggplot2::ggplotGrob(p_left), ggplot2::ggplotGrob(p_eff),
    ggplot2::ggplotGrob(p_right), ncol = 3,
    widths = grid::unit.c(grid::unit(lw, "in"), grid::unit(ew, "in"),
                          grid::unit(rw, "in")))
  fw_in <- lw + ew + rw
  fh_in <- (ylim[2] - ylim[1]) * row_in + 0.5

  spec <- dr_file_spec(file.name, filetype, path, width, height,
                       what = "easydr_loo()")
  title_in <- if (is.null(main)) 0 else 1.9 * 12 / 72
  side     <- position %in% c("left", "right")
  if (side) {
    w_in <- 6.2 + fw_in + 0.3
    size <- dr_fig_size(w_in, max(4.2, fh_in) + title_in + 0.12, max.w = 22)
  } else {
    w_in <- max(7.6, fw_in + 0.3)
    size <- dr_fig_size(w_in, 3.8 + title_in + fh_in + 0.2, max.w = 22)
  }
  con  <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)

  fh <- grid::unit(fh_in, "in")
  fwu <- grid::unit(fw_in, "in")
  g <- switch(position,
    top = dr_draw_stack(p, NULL, title = main, footer = fg, footer_h = fh),
    bottom = dr_draw_stack(gridExtra::arrangeGrob(
      grobs = list(fg, ggplot2::ggplotGrob(p)), ncol = 1,
      heights = grid::unit.c(fh + grid::unit(10, "pt"),
                             grid::unit(1, "null"))), NULL, title = main),
    left = dr_draw_stack(gridExtra::arrangeGrob(
      grobs = list(ggplot2::ggplotGrob(p), fg), ncol = 2,
      widths = grid::unit.c(grid::unit(1, "null"), fwu)), NULL, title = main),
    right = dr_draw_stack(gridExtra::arrangeGrob(
      grobs = list(fg, ggplot2::ggplotGrob(p)), ncol = 2,
      widths = grid::unit.c(fwu, grid::unit(1, "null"))), NULL, title = main))
  invisible(list(table = tab, curves = curves, grob = g))
}
