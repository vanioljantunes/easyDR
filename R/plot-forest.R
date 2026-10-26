#' Forest plot of the study trends and the pooled one
#'
#' @description
#' One row per study, laid out the way a published forest is. The study label
#' and what it brought sit on the left, a square sized by its weight carries
#' its confidence interval in the middle, and the value and the weight sit on
#' the right. The pooled estimate is the diamond underneath, with the
#' heterogeneity below it.
#'
#' @section What each row is:
#' Each row is the study's own slope per unit of dose, from the two-stage
#' first stage. A trend across its exposure categories is what a dose-response
#' study contributes. A shape with more than one term has no single slope, so
#' the row becomes the predicted effect at one dose against the fit's
#' reference dose, and the heading names that dose. Pass `at` to choose it,
#' for any shape.
#'
#' @section Why some fits cannot be drawn this way:
#' A one-stage fit estimates every study inside a single model and never
#' produces a per-study coefficient vector. Nothing can go on the rows, and
#' drawing the pooled line once per study would look like data, so this
#' stops.
#'
#' @param x a fitted result from [easydr()].
#' @param at a dose at which to compare studies, instead of their trends.
#'   Needed in effect for curves, where no single trend exists; the default is
#'   the first dose the fit already reports.
#' @param main plot title, drawn left-aligned above the figure. `NULL` writes
#'   one from the fit; `FALSE` draws none.
#' @param xlab label under the effect axis. `NULL` writes one from `sm`.
#' @param xlim limits for the effect axis, on the reported scale, such as
#'   `c(0.9, 1.5)` for a ratio.
#' @param order `"data"` keeps the order of the table, `"effect"` sorts by the
#'   estimate, `"weight"` sorts by weight, `"study"` sorts by label.
#' @param digits digits for the value column.
#' @param legend optional free text drawn under the figure, one line per
#'   element, in the manner of an added row.
#' @param columns `TRUE`, the default, shows the dose levels and dose range
#'   each study brought. `FALSE` shows the study label alone.
#' @param row.height height of one study row, in inches. Fixed rather than
#'   stretched, so a forest of six studies and one of thirty print with the
#'   same spacing.
#' @param effect.width width of the interval panel, in inches. The two text
#'   panels take the width their characters need and this takes the rest, so
#'   the figure has one width whatever the device is.
#'
#' @return The plotted table, invisibly: one row per study with `est`, `se`,
#'   `lb`, `ub` and `weight`, plus the `pooled` estimate and the `grob` that
#'   was drawn. Returned so the numbers behind the picture can be checked.
#'
#' @seealso [easydr_curve()] for the fitted dose-response itself, and
#'   [easydr_loo()] for what happens when a row is removed.
#'
#' @examples
#' fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
#'               dose.unit = "g/day")
#' easydr_forest(fit)
#'
#' # A curve has no single trend, so the rows become one dose.
#' rcs <- easydr(alcohol_cancer, sm = "IRR", shape = "rcs",
#'               dose.unit = "g/day")
#' easydr_forest(rcs, at = 25)
#'
#' @export
easydr_forest <- function(x, at = NULL, main = NULL, xlab = NULL, xlim = NULL,
                          order = c("data", "effect", "weight", "study"),
                          digits = 3, legend = NULL, columns = TRUE,
                          row.height = 0.26, effect.width = 3.1) {
  order <- match.arg(order)
  e   <- dr_study_effects(x, at, what = "a forest plot")
  tab <- e$table
  tab <- switch(order,
                data   = tab,
                effect = tab[order(tab$est), , drop = FALSE],
                weight = tab[order(-tab$weight), , drop = FALSE],
                study  = tab[order(tab$study), , drop = FALSE])
  k    <- nrow(tab)
  back <- if (e$ratio) exp else identity
  fmt  <- function(v) format(round(back(v), digits), nsmall = digits)

  # Rows run top to bottom: a header, the studies, then the pooled diamond
  # with the heterogeneity line under it.
  header_y <- k + 1
  study_y  <- rev(seq_len(k))
  pooled_y <- 0
  het_y    <- -1
  ylim     <- c(het_y - 0.5, header_y + 0.5)
  zebra    <- dr_zebra(study_y)

  lvl  <- format(100 * e$level)
  vals <- ifelse(is.na(tab$est), "not estimable",
                 paste0(fmt(tab$est), " (", fmt(tab$lb), " to ",
                        fmt(tab$ub), ")"))
  pooled_val <- paste0(fmt(e$pooled$est), " (", fmt(e$pooled$lb), " to ",
                       fmt(e$pooled$ub), ")")

  rows_face <- c(rep("plain", k), "bold")
  lab_cols  <- list(Study = c(tab$study, "Pooled"))
  lab_right <- FALSE
  if (isTRUE(columns)) {
    # "Levels" on its own told a reader nothing. The column counts the
    # exposure categories a study reported, and the range needs its unit or
    # the two numbers are bare.
    unit <- if (is.null(x$dose_unit)) "" else paste0(" (", x$dose_unit, ")")
    lab_cols <- c(lab_cols, list(
      `Dose levels` = c(as.character(tab$levels), ""),
      x = c(paste0(signif(tab$dose_min, 3), " to ",
                   signif(tab$dose_max, 3)), "")))
    names(lab_cols)[3L] <- paste0("Dose range", unit)
    lab_right <- c(FALSE, TRUE, TRUE)
  }

  het <- e$heterogeneity
  het_text <- if (isTRUE(het$estimable)) {
    sprintf("I2 = %.1f%%, Q = %.2f (df %d), %s",
            het$i2, het$Q, het$df, dr_p_lab(het$p))
  } else {
    paste0("I2: not estimable (", het$reason, ")")
  }

  het_note <- list(list(y = het_y, label = het_text))
  p_labels <- dr_text_panel(
    lab_cols, lab_right, ypos = c(study_y, pooled_y), face = rows_face,
    header_y = header_y, ylim = ylim, zebra = zebra, notes = het_note)

  val_cols <- list(`x` = c(vals, pooled_val),
                   Weight = c(paste0(format(round(tab$weight, 1), nsmall = 1),
                                     "%"), "100.0%"))
  names(val_cols)[1L] <- paste0(x$sm, " (", lvl, "% CI)")
  p_values <- dr_text_panel(
    val_cols, c(FALSE, TRUE), ypos = c(study_y, pooled_y), face = rows_face,
    header_y = header_y, ylim = ylim, zebra = zebra)

  # The effect panel. Squares have area proportional to weight, which is what
  # a reader reads off a square, and the diamond spans the pooled interval.
  lo <- min(c(tab$lb, e$pooled$lb, 0), na.rm = TRUE)
  hi <- max(c(tab$ub, e$pooled$ub, 0), na.rm = TRUE)
  if (!is.null(xlim)) {
    v  <- if (e$ratio) log(xlim) else xlim
    lo <- min(v); hi <- max(v)
  }
  pad <- 0.05 * (hi - lo)
  lo  <- lo - pad; hi <- hi + pad

  pts <- tab[!is.na(tab$est), , drop = FALSE]
  pts$ypos <- study_y[!is.na(tab$est)]
  pts$size <- 1.4 + 3.4 * sqrt(pts$weight / max(pts$weight, na.rm = TRUE))

  dh <- 0.34
  diamond <- data.frame(
    xx = c(e$pooled$lb, e$pooled$est, e$pooled$ub, e$pooled$est),
    yy = pooled_y + c(0, dh, 0, -dh))

  ticks <- dr_plot_ticks(lo, hi, e$ratio)
  keep  <- ticks$at >= lo & ticks$at <= hi
  if (is.null(xlab)) xlab <- e$label

  p_effect <- ggplot2::ggplot(pts, ggplot2::aes(x = .data$est,
                                                y = .data$ypos)) +
    zebra +
    ggplot2::geom_vline(xintercept = e$pooled$est, linetype = "dashed",
                        colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.4) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lb, xmax = .data$ub),
                           width = 0.26, orientation = "y", linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(size = .data$size), shape = 15L) +
    ggplot2::geom_polygon(data = diamond, inherit.aes = FALSE,
                          ggplot2::aes(x = .data$xx, y = .data$yy),
                          fill = "black") +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(breaks = ticks$at[keep],
                                labels = ticks$labels[keep]) +
    ggplot2::coord_cartesian(xlim = c(lo, hi), ylim = ylim) +
    ggplot2::labs(x = xlab) +
    dr_gg_theme(axis = TRUE)

  if (is.null(main)) {
    main <- paste0("Dose-response trend per study (", x$sm_label, ")")
  }
  if (identical(main, FALSE)) main <- NULL

  # The text panels take exactly the width their characters need, in inches,
  # and the effect panel takes whatever is left. Relative widths would spread
  # the columns of a narrow table across half the figure.
  lab_w <- dr_panel_width(lab_cols, notes = het_note)
  val_w <- dr_panel_width(val_cols)
  body <- gridExtra::arrangeGrob(
    p_labels, p_effect, p_values, ncol = 3,
    widths = grid::unit.c(grid::unit(lab_w, "inches"),
                          grid::unit(effect.width, "inches"),
                          grid::unit(val_w, "inches")),
    padding = grid::unit(0, "line"))

  body_h <- grid::unit((ylim[2] - ylim[1]) * row.height, "inches") +
    grid::unit(0.40, "inches")
  if (is.null(legend) && isTRUE(columns)) {
    legend <- paste0("Dose levels: exposure categories the study reported. ",
                     "Dose range: its lowest to its highest dose.")
  }
  g <- dr_draw_stack(body, body_h, title = main, footer = legend,
                     body_w = grid::unit(lab_w + effect.width + val_w,
                                         "inches"))

  # The composite has a fixed width, so the caller is told what it is: a
  # device narrower than this clips the first and last columns.
  invisible(list(table = tab, pooled = e$pooled, label = e$label, at = e$at,
                 heterogeneity = het, grob = g,
                 width = lab_w + effect.width + val_w,
                 height = (ylim[2] - ylim[1]) * row.height + 0.40))
}
