# Drawn by easydr_metareg(), and internal for the same reason it is.

#' Bubble plot of a dose-response meta-regression
#' @keywords internal
#' @noRd
easydr_bubble <- function(x, mod = NULL, main = NULL, xlab = NULL,
                          ylab = NULL, ylim = NULL, ci = TRUE, labels = FALSE,
                          col = "#4E7CA8", col.line = "#C0392B",
                          legend = NULL,
                          filetype = c("png", "pdf", "tiff"),
                          path = NULL, file.name = NULL,
                          width = NULL, height = NULL) {
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
  tab <- dr_study_effects(base_fit, NULL, what = "a bubble plot")$table

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
  if (identical(main, FALSE)) main <- NULL

  # The fitted value and its standard error come from the regression. For a
  # one-term shape the coefficients are the trend at the reference level plus
  # the change in trend, so the model matrix gives both directly.
  b  <- coef(x$fit)
  V  <- as.matrix(vcov(x$fit))
  zc <- stats::qnorm(1 - (1 - x$level) / 2)
  predict_at <- function(values) {
    nd <- st[rep(1L, length(values)), , drop = FALSE]
    nd[[mod]] <- values
    mm <- stats::model.matrix(stats::reformulate(x$moderators), data = nd)
    k   <- seq_len(ncol(mm))
    est <- as.numeric(mm %*% b[k])
    se  <- sqrt(pmax(0, rowSums((mm %*% V[k, k, drop = FALSE]) * mm)))
    data.frame(trend = est, lb = est - zc * se, ub = est + zc * se)
  }

  tab$area <- 1.6 + 7.0 * sqrt(tab$weight / max(tab$weight))

  if (continuous) {
    grid_v <- seq(min(tab$mod), max(tab$mod), length.out = 80)
    fitted <- cbind(data.frame(mod = grid_v), predict_at(grid_v))
    if (is.null(ylim)) {
      ylim <- dr_pad_range(
        range(back(c(tab$est, fitted$trend,
                     if (isTRUE(ci)) c(fitted$lb, fitted$ub))), na.rm = TRUE),
        x$ratio, 0.10)
    }
    p <- ggplot2::ggplot()
    if (isTRUE(ci)) {
      p <- p + ggplot2::geom_ribbon(
        data = fitted, ggplot2::aes(x = .data$mod, ymin = back(.data$lb),
                                    ymax = back(.data$ub)), fill = "grey85")
    }
    p <- p +
      ggplot2::geom_hline(yintercept = if (x$ratio) 1 else 0,
                          linetype = "dashed", colour = "grey45",
                          linewidth = 0.4) +
      ggplot2::geom_line(data = fitted,
                         ggplot2::aes(x = .data$mod, y = back(.data$trend)),
                         colour = col.line, linewidth = 0.9) +
      ggplot2::geom_point(data = tab,
                          ggplot2::aes(x = .data$mod, y = back(.data$est),
                                       size = .data$area),
                          shape = 21L, fill = scales::alpha(col, 0.55),
                          colour = col, stroke = 0.7) +
      ggplot2::scale_size_identity()
    xpos <- tab$mod
  } else {
    lv        <- unique(as.character(tab$mod))
    fitted    <- cbind(data.frame(mod = lv), predict_at(lv))
    fitted$xi <- seq_along(lv)
    # Studies in one level would otherwise be drawn on the same x. Spread them
    # evenly across the column instead of jittering at random, so the picture
    # is the same every time it is drawn.
    xpos     <- match(as.character(tab$mod), lv)
    xpos     <- xpos + dr_spread(xpos)
    tab$xpos <- xpos
    if (is.null(ylim)) {
      # The interval on each level mean is drawn, so it has to fit: scaling to
      # the estimates alone clips the very thing ci = TRUE added.
      ylim <- dr_pad_range(
        range(back(c(tab$est, fitted$trend,
                     if (isTRUE(ci)) c(fitted$lb, fitted$ub))), na.rm = TRUE),
        x$ratio, 0.10)
    }
    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = if (x$ratio) 1 else 0,
                          linetype = "dashed", colour = "grey45",
                          linewidth = 0.4)
    if (isTRUE(ci)) {
      # Capped, and in the colour of the mean it belongs to: a bare grey
      # vertical in the middle of a column reads as an axis line.
      p <- p + ggplot2::geom_errorbar(
        data = fitted, ggplot2::aes(x = .data$xi, ymin = back(.data$lb),
                                    ymax = back(.data$ub)),
        width = 0.12, colour = col.line, linewidth = 0.6)
    }
    p <- p +
      ggplot2::geom_segment(data = fitted,
                            ggplot2::aes(x = .data$xi - 0.3,
                                         xend = .data$xi + 0.3,
                                         y = back(.data$trend),
                                         yend = back(.data$trend)),
                            colour = col.line, linewidth = 1.1) +
      ggplot2::geom_point(data = tab,
                          ggplot2::aes(x = .data$xpos, y = back(.data$est),
                                       size = .data$area),
                          shape = 21L, fill = scales::alpha(col, 0.55),
                          colour = col, stroke = 0.7) +
      ggplot2::scale_size_identity() +
      ggplot2::scale_x_continuous(breaks = seq_along(lv), labels = lv,
                                  limits = c(0.55, length(lv) + 0.45))
  }

  if (isTRUE(labels)) {
    # Above the bubble, not beside it, so a label cannot sit on the marker it
    # names however the bubbles are sized.
    lab_df <- data.frame(x = xpos, y = back(tab$est), label = tab$study,
                         nudge = tab$area)
    p <- p + ggplot2::geom_text(
      data = lab_df, ggplot2::aes(x = .data$x, y = .data$y,
                                  label = .data$label),
      vjust = -1.4, size = 2.6, colour = "grey25")
  }

  p <- p +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   plot.margin = ggplot2::margin(6, 10, 2, 6))
  p <- if (x$ratio) {
    p + ggplot2::scale_y_log10(limits = ylim)
  } else {
    p + ggplot2::scale_y_continuous(limits = ylim)
  }

  test <- if (isTRUE(x$test$estimable)) {
    paste0("Between-moderator test: X2 = ", round(x$test$chi2, 2),
           " on ", x$test$df, " df, ", dr_p_lab(x$test$p), ".")
  } else {
    "Between-moderator test: not estimable."
  }
  foot <- c(paste(test, "Bubble area is the weight the study carries."),
            "An association across studies, not a randomised comparison.",
            legend)

  spec <- dr_file_spec(file.name, filetype, path, width, height,
                       what = "easydr_bubble()")
  title_in <- if (is.null(main)) 0 else 1.9 * 12 / 72
  foot_in  <- (1.35 * 9 * length(foot) + 8) / 72
  size <- dr_fig_size(7.2, min(max(0.52 * 7.2, 3.2), 4.4) + title_in +
                        foot_in + 0.22)
  con  <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)

  g <- dr_draw_stack(p, NULL, title = main, footer = foot)

  invisible(list(table = tab, fitted = fitted, moderator = mod,
                 test = x$test, grob = g))
}
