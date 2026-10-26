# What every plot in the package needs, and the house style they share.
#
# A dose-response fit has no single effect per study the way a pairwise
# meta-analysis does. What it has, after the two-stage first stage, is one
# coefficient vector per study with its covariance. The forest, the funnel and
# the bubble are all built from those. A one-stage fit never produces them, so
# those plots say so instead of inventing something to draw.
#
# The layout follows easyDTA: ggplot panels composed side by side with
# gridExtra, zebra-striped study rows, an in-data bold header row so column
# titles line up across panels whatever the axis areas do, and a fixed height
# per row so a forest of six studies and one of thirty have the same spacing.

#' One estimate per study, on the scale the model was fitted on
#'
#' For a one-term shape the estimate is the study's own trend per unit of dose,
#' which is what the coefficient already is. For anything with more than one
#' term there is no single trend, so the estimate is the predicted effect at a
#' stated dose against the fit's reference dose, and the dose is reported
#' rather than assumed.
#'
#' The weight is the random-effects inverse-variance weight for the same
#' contrast: the within-study variance plus the between-study variance the fit
#' estimated, both projected onto the contrast. A common-effect fit has a zero
#' between-study term, which is the point of it.
#'
#' @keywords internal
#' @noRd
dr_study_effects <- function(x, at = NULL, what = "this plot") {
  if (!inherits(x, "easydr_result")) {
    stop("easydr: ", what, " takes a fitted result from easydr().",
         call. = FALSE)
  }
  bi <- x$fit$bi
  Si <- x$fit$Si
  if (is.null(bi) || is.null(Si) || !is.matrix(bi)) {
    stop("easydr: ", what, " needs one estimate per study, and only the ",
         "two-stage first stage produces them. This fit is one-stage, which ",
         "estimates every study inside a single model. Refit with ",
         "proc = \"2stage\" if you want it.", call. = FALSE)
  }

  frame   <- x$frame
  studies <- unique(as.character(frame$studlab))
  ids     <- vapply(studies, function(s)
    frame$id[as.character(frame$studlab) == s][1L], numeric(1))
  xref <- x$reference$curve
  unit <- if (is.null(x$dose_unit)) "dose unit" else x$dose_unit

  one_term <- x$shape$terms == 1L
  if (one_term && is.null(at)) {
    cvec    <- matrix(1, nrow = 1L, ncol = 1L)
    label   <- dr_label(x$sm, x$dose_unit)
    at_used <- NA_real_
  } else {
    if (is.null(at)) {
      # The first dose the fit already reports, so the plot and the printed
      # table are talking about the same thing.
      cand    <- x$predicted$dose[!x$predicted$is_ref]
      at_used <- if (length(cand) > 0L) cand[1L] else max(frame$dose)
    } else {
      if (!is.numeric(at) || length(at) != 1L || !is.finite(at)) {
        stop("easydr: `at` must be a single dose, such as at = 10.",
             call. = FALSE)
      }
      at_used <- at
    }
    cvec  <- dr_contrast(x$shape, at_used, xref)
    label <- paste0(x$sm, " at ", at_used, " ", unit, " against ", xref)
  }

  Psi <- x$Psi
  if (is.null(Psi)) Psi <- matrix(0, ncol(bi), ncol(bi))
  Psi  <- as.matrix(Psi)
  tau2 <- as.numeric(cvec %*% Psi %*% t(cvec))

  rows <- lapply(seq_along(studies), function(i) {
    id  <- ids[i]
    sub <- frame[as.character(frame$studlab) == studies[i], , drop = FALSE]
    base <- data.frame(
      study    = studies[i],
      levels   = length(unique(sub$dose)),
      dose_min = min(sub$dose),
      dose_max = max(sub$dose),
      stringsAsFactors = FALSE
    )
    if (is.na(id) || id > nrow(bi)) {
      return(cbind(base, est = NA_real_, se = NA_real_, var_w = NA_real_))
    }
    b  <- as.numeric(bi[id, ])
    V  <- as.matrix(Si[[id]])
    vw <- as.numeric(cvec %*% V %*% t(cvec))
    cbind(base, est = as.numeric(cvec %*% b), se = sqrt(vw), var_w = vw)
  })
  tab <- do.call(rbind, rows)

  zcrit <- stats::qnorm(1 - (1 - x$provenance$level) / 2)
  tab$lb <- tab$est - zcrit * tab$se
  tab$ub <- tab$est + zcrit * tab$se
  w <- 1 / (tab$var_w + max(tau2, 0))
  tab$weight <- 100 * w / sum(w, na.rm = TRUE)
  rownames(tab) <- NULL

  # The pooled value for the same contrast, so the diamond and the squares are
  # answering one question.
  pooled_est <- as.numeric(cvec %*% as.numeric(coef(x$fit)))
  pooled_var <- as.numeric(cvec %*% as.matrix(vcov(x$fit)) %*% t(cvec))
  pooled <- list(est = pooled_est, se = sqrt(pooled_var),
                 lb = pooled_est - zcrit * sqrt(pooled_var),
                 ub = pooled_est + zcrit * sqrt(pooled_var))

  list(table = tab, pooled = pooled, label = label, at = at_used,
       ratio = x$ratio, sm = x$sm, unit = unit, level = x$provenance$level,
       heterogeneity = x$heterogeneity,
       one_term = one_term && is.na(at_used))
}

#' The pooled effect for one contrast, from the fitted model alone
#'
#' `dr_study_effects()` needs the per-study coefficients, which a one-stage fit
#' never produces. The curve has to be able to report its own result whatever
#' the pooling stage was, so this reads `coef()` and `vcov()` and nothing else.
#'
#' A one-term shape with no dose named reports its slope. Anything else reports
#' the predicted effect at a dose, and says which dose that was.
#'
#' @keywords internal
#' @noRd
dr_pooled_at <- function(x, at = NULL) {
  one_term <- x$shape$terms == 1L
  xref <- x$reference$curve
  unit <- if (is.null(x$dose_unit)) "dose unit" else x$dose_unit

  if (one_term && is.null(at)) {
    cvec    <- matrix(1, nrow = 1L, ncol = 1L)
    label   <- dr_label(x$sm, x$dose_unit)
    at_used <- NA_real_
  } else {
    if (is.null(at)) {
      cand    <- x$predicted$dose[!x$predicted$is_ref]
      at_used <- if (length(cand) > 0L) cand[1L] else max(x$frame$dose)
    } else {
      at_used <- at
    }
    cvec  <- dr_contrast(x$shape, at_used, xref)
    label <- paste0(x$sm, " at ", at_used, " ", unit)
  }

  b  <- as.numeric(coef(x$fit))
  V  <- as.matrix(vcov(x$fit))
  est <- as.numeric(cvec %*% b)
  se  <- sqrt(as.numeric(cvec %*% V %*% t(cvec)))
  z   <- stats::qnorm(1 - (1 - x$provenance$level) / 2)
  list(est = est, se = se, lb = est - z * se, ub = est + z * se,
       p = 2 * stats::pnorm(-abs(est / se)), label = label, at = at_used)
}

#' Leave a margin round a range, multiplicatively on a ratio scale
#'
#' A log axis needs a multiplicative margin: adding a constant to a ratio makes
#' the gap at the bottom much larger than the gap at the top.
#'
#' @keywords internal
#' @noRd
dr_pad_range <- function(r, ratio, frac = 0.06) {
  if (!all(is.finite(r))) return(r)
  if (isTRUE(ratio)) {
    lr <- log(r)
    d  <- max(diff(lr), 1e-8) * frac
    return(exp(c(lr[1L] - d, lr[2L] + d)))
  }
  d <- max(diff(r), 1e-8) * frac
  c(r[1L] - d, r[2L] + d)
}

#' Spread the points of one column evenly, the same way every time
#'
#' Random jitter would move the picture between two runs of the same script,
#' which is the last thing a figure in a review should do.
#'
#' @keywords internal
#' @noRd
dr_spread <- function(group, width = 0.26) {
  out <- numeric(length(group))
  for (g in unique(group)) {
    i <- which(group == g)
    out[i] <- if (length(i) == 1L) 0 else seq(-width, width,
                                              length.out = length(i))
  }
  out
}

#' Egger's regression test for small-study effects
#'
#' Regresses each study's standard normal deviate on its precision. A non-zero
#' intercept is the asymmetry: small studies sitting away from the pooled value.
#'
#' The test is for the trend estimates, one per study, which the two-stage
#' first stage produces. It is underpowered below about ten studies, and that
#' is reported alongside rather than left to be remembered.
#'
#' @keywords internal
#' @noRd
dr_egger <- function(est, se, k_warn = 10L) {
  ok  <- is.finite(est) & is.finite(se) & se > 0
  est <- est[ok]
  se  <- se[ok]
  k   <- length(est)
  if (k < 3L) {
    return(list(estimable = FALSE, k = k,
                reason = paste0("Egger's test needs at least 3 studies with a ",
                                "standard error, and ", k, " ",
                                if (k == 1L) "is" else "are", " available")))
  }
  # Every study the same size means there is no small-study contrast, and the
  # precision column is constant: the regression is rank deficient and has no
  # intercept to read. That is a real answer about the data, not a failure.
  if (length(unique(se)) < 2L) {
    return(list(estimable = FALSE, k = k,
                reason = paste0("every study has the same standard error, so ",
                                "there is no contrast between large and small ",
                                "studies for the regression to pick up")))
  }
  m  <- stats::lm(I(est / se) ~ I(1 / se))
  cf <- summary(m)$coefficients
  if (nrow(cf) < 2L || anyNA(cf[, 1L])) {
    return(list(estimable = FALSE, k = k,
                reason = "the precisions do not vary enough to fit the line"))
  }
  list(estimable = TRUE, k = k,
       intercept = cf[1L, 1L], se = cf[1L, 2L],
       t = cf[1L, 3L], df = k - 2L, p = cf[1L, 4L],
       slope = cf[2L, 1L],
       underpowered = k < k_warn, k_warn = k_warn)
}

#' Format a p value the way a results table does
#' @keywords internal
#' @noRd
dr_p <- function(p, digits = 3) {
  if (length(p) != 1L || !is.finite(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  format(round(p, digits), nsmall = digits)
}

#' A p value with its own relational operator
#'
#' `dr_p()` returns "< 0.001" for a small p, so a caller that writes "p = "
#' in front of it produces "p = < 0.001". This writes the whole label.
#'
#' @keywords internal
#' @noRd
dr_p_lab <- function(p, digits = 3) {
  v <- dr_p(p, digits)
  if (startsWith(v, "<")) paste("p", v) else paste("p =", v)
}

#' The theme every panel shares
#'
#' Panels are composed side by side, so they must reserve identical bottom
#' space for the axis ink or the study rows will not line up. A panel with no
#' axis to show draws it in `NA` colour rather than dropping it.
#'
#' @keywords internal
#' @noRd
dr_gg_theme <- function(axis = FALSE, axis_text_size = 8) {
  base <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.y    = ggplot2::element_blank(),
      axis.text.y     = ggplot2::element_blank(),
      axis.ticks.y    = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      panel.border    = ggplot2::element_blank(),
      plot.margin     = ggplot2::margin(t = 5.5, r = 0, b = 2, l = 0),
      legend.position = "none"
    )
  # The axis title reserves its height on every panel, visible only on the one
  # that has an axis. Dropping it on the others makes their panel areas taller
  # and the study rows stop lining up across the figure.
  ink <- if (isTRUE(axis)) "black" else NA
  base + ggplot2::theme(
    axis.text.x  = ggplot2::element_text(colour = ink, size = axis_text_size),
    axis.ticks.x = ggplot2::element_line(colour = ink, linewidth = 0.4),
    axis.line.x  = ggplot2::element_line(colour = ink, linewidth = 0.4),
    axis.title.x = ggplot2::element_text(colour = ink, size = 9,
                                         margin = ggplot2::margin(t = 4)),
    axis.ticks.length = grid::unit(3, "pt"))
}

#' Stripe every other study row, and never a summary row
#' @keywords internal
#' @noRd
dr_zebra <- function(ypos) {
  df <- data.frame(ypos = ypos[seq_along(ypos) %% 2L == 1L])
  ggplot2::geom_rect(
    data = df, inherit.aes = FALSE,
    ggplot2::aes(xmin = -Inf, xmax = Inf,
                 ymin = .data$ypos - 0.5, ymax = .data$ypos + 0.5),
    fill = "grey92")
}

#' A panel of text columns, headed in bold
#'
#' Column positions come from character widths, so a long study label cannot
#' push a numeric column off the panel. Numeric columns are right-aligned
#' against their own edge and text columns left-aligned against their start.
#'
#' @keywords internal
#' @noRd
dr_text_panel <- function(cols, right, ypos, face, header_y, ylim, zebra,
                          notes = NULL, cell = 3.0, hdr = 3.2, xlab = " ") {
  w <- vapply(seq_along(cols), function(i)
    max(nchar(c(names(cols)[i], as.character(cols[[i]]))), na.rm = TRUE) +
      if (right[i]) 1.5 else 2.5, numeric(1))
  edge  <- cumsum(w) / sum(w)
  start <- c(0, utils::head(edge, -1))
  xs    <- ifelse(right, edge, start)

  cells <- do.call(rbind, lapply(seq_along(cols), function(i) data.frame(
    x = xs[i], ypos = ypos, label = as.character(cols[[i]]), face = face,
    hjust = if (right[i]) 1 else 0, stringsAsFactors = FALSE)))
  cells <- cells[!is.na(cells$label) & cells$label != "", , drop = FALSE]
  heads <- data.frame(x = xs, ypos = header_y, label = names(cols),
                      hjust = ifelse(right, 1, 0), stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(cells, ggplot2::aes(y = .data$ypos)) +
    zebra +
    ggplot2::geom_text(ggplot2::aes(x = .data$x, label = .data$label,
                                    fontface = .data$face,
                                    hjust = .data$hjust), size = cell) +
    ggplot2::geom_text(data = heads,
                       ggplot2::aes(x = .data$x, label = .data$label,
                                    hjust = .data$hjust),
                       fontface = "bold", size = hdr)
  for (n in notes) {
    p <- p + ggplot2::annotate("text", x = 0, y = n$y, label = n$label,
                               hjust = 0, size = cell)
  }
  p +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim) +
    ggplot2::labs(x = xlab) +
    dr_gg_theme(axis = FALSE)
}

#' Width of a text panel in inches, from the characters it has to hold
#' @keywords internal
#' @noRd
dr_panel_width <- function(cols, per_char = 0.070, min_in = 1.2,
                           notes = NULL) {
  chars <- sum(vapply(seq_along(cols), function(i)
    max(nchar(c(names(cols)[i], as.character(cols[[i]]))), na.rm = TRUE) + 2,
    numeric(1)))
  # A note is written across the whole panel, so a long one sets the width on
  # its own rather than running off the edge.
  if (length(notes) > 0L) {
    chars <- max(chars, max(nchar(vapply(notes, function(n) n$label,
                                         character(1)))) + 2)
  }
  max(min_in, per_char * chars)
}

#' A bordered statistics table, drawn under a plot
#'
#' easyDTA puts the numbers behind a picture in a bordered table beneath it
#' rather than in a caption, so the figure can be read on its own. This is that
#' table.
#'
#' @keywords internal
#' @noRd
dr_stats_grob <- function(df, base_size = 9) {
  tt <- gridExtra::ttheme_minimal(
    base_size = base_size,
    core    = list(fg_params = list(hjust = 0.5, x = 0.5),
                   bg_params = list(fill = c("white", "grey97"),
                                    col = "grey60", lwd = 0.6)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0.5, x = 0.5),
                   bg_params = list(fill = "grey90", col = "grey60",
                                    lwd = 0.6)))
  gridExtra::tableGrob(df, rows = NULL, theme = tt)
}

#' How tall a grob actually is, in points
#'
#' Guessing a table's height from its row count gets it wrong as soon as the
#' padding or the font changes, and the band then overlaps the plot above it.
#' The grob knows, so it is asked.
#'
#' @keywords internal
#' @noRd
dr_grob_height <- function(g, pad = 10) {
  h <- try(grid::convertHeight(sum(g$heights), "points", valueOnly = TRUE),
           silent = TRUE)
  if (inherits(h, "try-error") || !is.finite(h)) h <- 110
  grid::unit(h + pad, "points")
}

#' Stack a title, a body and an optional footer, and draw them
#'
#' A forest has a fixed height per row, so its composite is centred in the
#' device at exactly the height it needs. Everything else fills the device,
#' which is `body_h = NULL`: the title and footer take their points and the
#' body takes the rest. Mixing the two would push a null unit into the
#' viewport height, and the body then collapses under the footer.
#'
#' @keywords internal
#' @noRd
dr_draw_stack <- function(body, body_h = NULL, title = NULL, footer = NULL,
                          footer_h = NULL, title_size = 12, footer_size = 9,
                          body_w = NULL) {
  fill    <- is.null(body_h)
  grobs   <- list(body)
  heights <- list(if (fill) grid::unit(1, "null") else body_h)
  if (!is.null(title)) {
    grobs   <- c(list(grid::textGrob(
      title, x = grid::unit(4, "pt"), hjust = 0,
      gp = grid::gpar(fontface = "bold", fontsize = title_size))), grobs)
    heights <- c(list(grid::unit(1.9 * title_size, "points")), heights)
  }
  if (!is.null(footer)) {
    if (is.character(footer)) {
      n <- length(strsplit(paste(footer, collapse = "\n"), "\n",
                           fixed = TRUE)[[1]])
      footer <- grid::textGrob(
        paste(footer, collapse = "\n"), x = grid::unit(4, "pt"),
        y = grid::unit(1, "npc"), hjust = 0, vjust = 1,
        gp = grid::gpar(fontsize = footer_size, lineheight = 1.25))
      if (is.null(footer_h)) {
        footer_h <- grid::unit(1.35 * footer_size * n + 8, "points")
      }
    }
    if (is.null(footer_h)) footer_h <- grid::unit(1, "null")
    grobs   <- c(grobs, list(footer))
    heights <- c(heights, list(footer_h))
  }
  hv <- do.call(grid::unit.c, heights)
  g  <- gridExtra::arrangeGrob(grobs = grobs, ncol = 1, heights = hv)
  grid::grid.newpage()
  # A composite of fixed-width panels is centred at that width rather than
  # stretched across the device, or the value column drifts away from the
  # interval it belongs to on a wide device.
  w <- if (is.null(body_w)) grid::unit(1, "npc") else body_w
  vp <- if (fill) {
    grid::viewport(width = w, height = grid::unit(1, "npc"))
  } else {
    grid::viewport(x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
                   just = "centre", width = w, height = Reduce(`+`, heights))
  }
  grid::pushViewport(vp)
  grid::grid.draw(g)
  grid::popViewport()
  invisible(g)
}

#' Axis ticks on the reported scale, for an axis drawn on the model scale
#' @keywords internal
#' @noRd
dr_plot_ticks <- function(lo, hi, ratio) {
  if (!ratio) {
    v <- pretty(c(lo, hi))
    return(list(at = v, labels = format(v)))
  }
  v <- pretty(c(exp(lo), exp(hi)))
  v <- v[v > 0]
  list(at = log(v), labels = format(v))
}
