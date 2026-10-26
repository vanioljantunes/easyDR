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
                        ylim = NULL, digits = 4) {
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

  if (is.null(ylab)) ylab <- out$sm_label
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

  # The table under the figure. The bold last row is the fit with nothing
  # omitted, so every row above it reads as a departure from it.
  fmt <- function(v) ifelse(is.na(v), "-",
                            format(round(v, digits), nsmall = digits))
  shift <- ifelse(is.na(tab$shift), "-",
                  paste0(ifelse(tab$shift >= 0, "+", ""),
                         format(round(tab$shift, digits), nsmall = digits)))
  pct <- function(v) ifelse(is.na(v), "-",
                            paste0(format(round(v, 1), nsmall = 1), "%"))
  full_i2 <- if (isTRUE(base$heterogeneity$estimable)) {
    pct(base$heterogeneity$i2)
  } else {
    "-"
  }

  cell <- data.frame(
    ` `         = c(ifelse(is.na(tab$key), "", tab$key), ""),
    Omitting    = c(tab$label, "None: the full fit"),
    est         = c(fmt(tab$estimate), fmt(out$full$est)),
    `95% CI`    = c(ifelse(is.na(tab$lb), "-",
                           paste0(fmt(tab$lb), " to ", fmt(tab$ub))),
                    paste0(fmt(out$full$lb), " to ", fmt(out$full$ub))),
    Shift       = c(shift, "-"),
    `I-squared` = c(pct(tab$i2), full_i2),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(cell)[3L] <- if (identical(out$quantity, "slope")) {
    paste0(out$sm, " per 1 ", unit)
  } else {
    paste0(out$sm, " at ", out$dose, " ", unit)
  }

  tg <- dr_stats_grob(cell)
  # Embolden the final row, the one nothing was left out of.
  last <- nrow(cell) + 1L
  for (i in seq_len(ncol(cell))) {
    idx <- which(tg$layout$t == last & tg$layout$l == i &
                   grepl("core-fg", tg$layout$name))
    # tableGrob already stores the face as `font`; setting `fontface` beside
    # it makes grid refuse to draw the cell at all.
    if (length(idx)) {
      tg$grobs[[idx]]$gp$fontface <- NULL
      tg$grobs[[idx]]$gp$font <- 2L
    }
  }

  g <- dr_draw_stack(p, NULL, title = main, footer = tg,
                     footer_h = dr_grob_height(tg))
  invisible(list(table = tab, curves = curves, grob = g))
}
