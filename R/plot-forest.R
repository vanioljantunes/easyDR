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
#' @param subgroup `TRUE` groups the rows by the subgroup the fit carries,
#'   `FALSE` lists every study in one block. The default, `NULL`, groups them
#'   when the fit carries a subgroup. A block is headed by the level on its own
#'   and closed by a subtotal: the model refitted within that level, drawn as a
#'   grey diamond, with its own heterogeneity beneath. The squares and the
#'   weights stay each study's share of the overall fit. See the `subgroup`
#'   argument of [easydr()].
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
#' fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
#'               dose.unit = "g/day")
#' easydr_forest(fit)
#'
#' # A curve has no single trend, so the rows become one dose.
#' rcs <- easydr(alcohol_cancer, sm = "RR", shape = "rcs",
#'               dose.unit = "g/day")
#' easydr_forest(rcs, at = 25)
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
easydr_forest <- function(x, at = NULL, main = NULL, xlab = NULL, xlim = NULL,
                          order = c("data", "effect", "weight", "study"),
                          subgroup = NULL, digits = 3, legend = NULL,
                          columns = TRUE, row.height = 0.26,
                          effect.width = 3.1,
                          filetype = c("png", "pdf", "tiff"),
                          path = NULL, file.name = NULL,
                          width = NULL, height = NULL) {
  order <- match.arg(order)
  # dr_study_effects() is what refuses a caller who passed something that is
  # not a fit, so it runs first: reading x$subgroup off a character vector
  # fails with a message about the `$` operator instead.
  e   <- dr_study_effects(x, at, what = "a forest plot")
  use_sub <- dr_use_subgroups(x, subgroup, "easydr_forest()")
  tab <- e$table
  tab <- switch(order,
                data   = tab,
                effect = tab[order(tab$est), , drop = FALSE],
                weight = tab[order(-tab$weight), , drop = FALSE],
                study  = tab[order(tab$study), , drop = FALSE])
  k    <- nrow(tab)
  back <- if (e$ratio) exp else identity
  fmt  <- function(v) format(round(back(v), digits), nsmall = digits)
  lvl  <- format(100 * e$level)

  # One row plan drives every panel. A forest with subgroups and a forest
  # without differ only in which rows the plan holds, so the three panels, the
  # widths and the zebra are built once instead of twice.
  mk <- function(kind, label, face = "plain", est = NA_real_, lb = NA_real_,
                 ub = NA_real_, weight = NA_real_, levels = "", dose = "",
                 note = NA_character_) {
    data.frame(kind = kind, label = label, face = face, est = est, lb = lb,
               ub = ub, weight = weight, levels = levels, dose = dose,
               note = note, stringsAsFactors = FALSE)
  }
  study_rows <- function(d, indent = "") {
    if (nrow(d) == 0L) return(NULL)
    mk("study", paste0(indent, d$study), "plain", d$est, d$lb, d$ub, d$weight,
       as.character(d$levels),
       paste0(signif(d$dose_min, 3), " to ", signif(d$dose_max, 3)))
  }
  het_line <- function(h) {
    if (isTRUE(h$estimable)) {
      sprintf("I2 = %.1f%%, Q = %.2f (df %d), %s",
              h$i2, h$Q, h$df, dr_p_lab(h$p))
    } else {
      paste0("I2: not estimable (", h$reason, ")")
    }
  }

  if (!use_sub) {
    plan <- rbind(study_rows(tab),
                  mk("pooled", "Pooled", "bold", e$pooled$est, e$pooled$lb,
                     e$pooled$ub, 100))
  } else {
    # subgroup.name = FALSE: the block header is the level on its own, never
    # "design = case-control". The column is named once, in the legend.
    sg   <- x$subgroup
    plan <- NULL
    for (lv in sg$levels) {
      members <- names(sg$by_study)[sg$by_study == lv]
      d <- tab[tab$study %in% members, , drop = FALSE]
      plan <- rbind(plan, mk("head", lv, "bold"), study_rows(d, "  "))
      f <- sg$fits[[lv]]
      if (is.null(f)) {
        reason <- sg$table$reason[sg$table$level == lv]
        plan <- rbind(plan, mk("note", "",
                               note = paste0("  not fitted: ", reason)))
      } else {
        pl <- dr_pooled_at(f, at)
        plan <- rbind(
          plan,
          mk("subtotal", "  Subtotal", "bold", pl$est, pl$lb, pl$ub,
             sum(d$weight, na.rm = TRUE)),
          mk("note", "", note = paste0("  ", het_line(f$heterogeneity))))
      }
    }
    plan <- rbind(plan, mk("pooled", "Overall", "bold", e$pooled$est,
                           e$pooled$lb, e$pooled$ub, 100))
  }

  # Rows run top to bottom: a header, the plan, then the overall heterogeneity
  # under the last row.
  n        <- nrow(plan)
  plan$y   <- rev(seq_len(n))
  header_y <- n + 1
  het_y    <- 0
  ylim     <- c(het_y - 0.5, header_y + 0.5)
  is_study <- plan$kind == "study"
  zebra    <- dr_zebra(plan$y[is_study])

  has_est <- !is.na(plan$est)
  vals <- ifelse(plan$kind %in% c("head", "note"), "",
                 ifelse(!has_est, "not estimable",
                        paste0(fmt(plan$est), " (", fmt(plan$lb), " to ",
                               fmt(plan$ub), ")")))
  wts <- ifelse(is.na(plan$weight), "",
                paste0(format(round(plan$weight, 1), nsmall = 1), "%"))

  lab_cols  <- list(Study = plan$label)
  lab_right <- FALSE
  if (isTRUE(columns)) {
    # "Levels" on its own told a reader nothing. The column counts the
    # exposure categories a study reported, and the range needs its unit or
    # the two numbers are bare.
    unit <- if (is.null(x$dose_unit)) "" else paste0(" (", x$dose_unit, ")")
    lab_cols <- c(lab_cols, list(`Dose levels` = plan$levels, x = plan$dose))
    names(lab_cols)[3L] <- paste0("Dose range", unit)
    lab_right <- c(FALSE, TRUE, TRUE)
  }

  het      <- e$heterogeneity
  het_text <- het_line(het)
  # A note is written across the label panel rather than into a cell, so a
  # long one cannot push the dose columns off the edge.
  notes <- c(
    lapply(which(plan$kind == "note"),
           function(i) list(y = plan$y[i], label = plan$note[i])),
    list(list(y = het_y, label = het_text)))

  p_labels <- dr_text_panel(
    lab_cols, lab_right, ypos = plan$y, face = plan$face,
    header_y = header_y, ylim = ylim, zebra = zebra, notes = notes)

  val_cols <- list(`x` = vals, Weight = wts)
  names(val_cols)[1L] <- paste0(x$sm, " (", lvl, "% CI)")
  p_values <- dr_text_panel(
    val_cols, c(FALSE, TRUE), ypos = plan$y, face = plan$face,
    header_y = header_y, ylim = ylim, zebra = zebra)

  # The effect panel. Squares have area proportional to weight, which is what
  # a reader reads off a square, and a diamond spans a pooled interval.
  lo <- min(c(plan$lb, 0), na.rm = TRUE)
  hi <- max(c(plan$ub, 0), na.rm = TRUE)
  if (!is.null(xlim)) {
    v  <- if (e$ratio) log(xlim) else xlim
    lo <- min(v); hi <- max(v)
  }
  pad <- 0.05 * (hi - lo)
  lo  <- lo - pad; hi <- hi + pad

  pts <- plan[is_study & has_est, , drop = FALSE]
  pts$ypos <- pts$y
  pts$size <- 1.4 + 3.4 * sqrt(pts$weight / max(pts$weight, na.rm = TRUE))

  dh   <- 0.34
  dia  <- plan[plan$kind %in% c("subtotal", "pooled") & has_est, ,
               drop = FALSE]
  diamond <- do.call(rbind, lapply(seq_len(nrow(dia)), function(i) {
    data.frame(xx = c(dia$lb[i], dia$est[i], dia$ub[i], dia$est[i]),
               yy = dia$y[i] + c(0, dh, 0, -dh),
               grp = i,
               # A subtotal is the level refitted, not the overall answer, so
               # it is drawn grey and the overall diamond stays black.
               fill = if (identical(dia$kind[i], "pooled")) "black" else
                 "grey45",
               stringsAsFactors = FALSE)
  }))

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
                          ggplot2::aes(x = .data$xx, y = .data$yy,
                                       group = .data$grp,
                                       fill = I(.data$fill))) +
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
  lab_w <- dr_panel_width(lab_cols, notes = notes)
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
  if (use_sub) {
    legend <- c(legend, paste0(
      "Grouped by ", x$subgroup$var, ". A subtotal is the model refitted ",
      "within that level, drawn as a grey diamond. Squares and weights are ",
      "each study's share of the overall fit."))
  }
  # The footer is wrapped to the width the figure actually has, and counted
  # into the height that is handed back. A caller opens a device at that
  # height, so a footer measured as one line and drawn as three is a footer
  # with its last two lines off the bottom of the page.
  total_w <- lab_w + effect.width + val_w
  if (!is.null(legend)) {
    legend <- unlist(lapply(legend, strwrap,
                            width = max(40L, floor(total_w / 0.068))))
  }
  foot_in  <- if (is.null(legend)) 0 else (1.35 * 9 * length(legend) + 8) / 72
  title_in <- if (is.null(main)) 0 else 1.9 * 12 / 72

  # The forest is the one figure whose size is exact rather than estimated:
  # the panels are in inches already and the rows have a fixed height, so the
  # device is opened at the size the composite will occupy.
  spec <- dr_file_spec(file.name, filetype, path, width, height,
                       what = "easydr_forest()")
  size <- dr_fig_size(total_w + 0.24,
                      dr_height_in(body_h) + title_in + foot_in + 0.20)
  con  <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)

  g <- dr_draw_stack(body, body_h, title = main, footer = legend,
                     body_w = grid::unit(total_w, "inches"))

  # The composite has a fixed width, so the caller is told what it is: a
  # device narrower than this clips the first and last columns.
  # Classed so printing it reads like the figure, with the values untouched.
  class(tab) <- c("easydr_contributions", "data.frame")
  invisible(list(table = tab, pooled = e$pooled, label = e$label, at = e$at,
                 heterogeneity = het, grob = g,
                 width = total_w,
                 height = (ylim[2] - ylim[1]) * row.height + 0.40 +
                   title_in + foot_in))
}
