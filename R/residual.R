#' Residual plots of a dose-response fit
#'
#' Plots the decorrelated residuals of the fit against dose, one panel per
#' shape. The reported effects of one study share its reference group, so
#' their residuals are correlated and a plot of raw residuals misleads; the
#' residuals here are decorrelated through the Cholesky factor of their
#' covariance, as dosresmeta's `gof()` returns them, following Discacciati,
#' Crippa and Orsini (2017) and their code at
#' <https://github.com/anddis/goodness-of-fit-meta-analysis>. A curve that
#' fits leaves residuals scattered around zero at every dose; a run of
#' residuals on one side of zero, which the smooth line traces, marks a dose
#' range the shape misses.
#'
#' @param x a fitted result from [easydr()].
#' @param shapes the shapes to plot, such as `shapes = c("rcs", "linear")`.
#'   The default is the shape `x` was fitted with. One shape draws one plot;
#'   two or more draw one panel per shape, each titled with the shape's name.
#' @param main a title. With one shape the default is no title; with two or
#'   more it is "Residual plots for different models in" followed by the
#'   name of the data `x` was fitted to. `FALSE` draws none.
#' @param smooth `TRUE`, the default, adds a lowess line through the
#'   residuals of each panel.
#' @param filetype,path,file.name,width,height where to save the figure, as
#'   in [easydr_curve()].
#'
#' @return A data frame, invisibly, with one row per non-reference dose and
#'   shape: `shape`, `studlab`, `dose` and `residual`. Called for the plot.
#'
#' @references Discacciati A, Crippa A, Orsini N. Goodness of fit tools for
#'   dose-response meta-analysis of binary outcomes. Research Synthesis
#'   Methods. 2017;8(2):149-160. doi:10.1002/jrsm.1194
#'
#' @examples
#' \donttest{
#' f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs")
#' easydr_residual(f)
#' easydr_residual(f, shapes = c("rcs", "linear", "quadratic"))
#' }
#' @export
easydr_residual <- function(x, shapes = x$shape$name, main = NULL,
                            smooth = TRUE,
                            filetype = c("png", "pdf", "tiff"),
                            path = NULL, file.name = NULL,
                            width = NULL, height = NULL) {
  dr_require_single(x, "easydr_residual()")
  if (!is.character(shapes) || !length(shapes)) {
    stop("easydr: `shapes` names one or more shapes, such as ",
         "shapes = c(\"rcs\", \"linear\").", call. = FALSE)
  }
  shapes <- unique(shapes)
  fits <- lapply(shapes, function(s) {
    if (identical(s, x$shape$name)) x else dr_refit_shape(x, s)
  })
  names(fits) <- shapes
  rows <- do.call(rbind, lapply(shapes, function(s) {
    f <- fits[[s]]
    if (is.null(f)) {
      warning("easydr: the ", s, " shape could not be fitted to this table ",
              "and is left out.", call. = FALSE)
      return(NULL)
    }
    dr_residual_rows(f, s)
  }))
  if (is.null(rows) || !nrow(rows)) {
    stop("easydr: none of those shapes gave residuals for this table.",
         call. = FALSE)
  }
  rows$panel <- factor(dr_shape_title(rows$shape),
                       levels = unique(dr_shape_title(rows$shape)))

  many <- nlevels(rows$panel) > 1L
  if (is.null(main)) {
    main <- if (many) {
      paste("Residual plots for different models in", dr_data_name(x))
    } else NULL
  }
  if (isFALSE(main)) main <- NULL

  unit <- if (is.null(x$dose_unit)) "" else paste0(" (", x$dose_unit, ")")
  p <- ggplot2::ggplot(rows, ggplot2::aes(x = .data$dose, y = .data$residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_point(shape = 21, size = 1.8, fill = "grey85",
                        colour = "black", stroke = 0.4) +
    ggplot2::labs(x = paste0("Dose", unit), y = "Decorrelated residuals") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 10, face = "bold"),
                   plot.margin = ggplot2::margin(6, 10, 2, 6))
  if (isTRUE(smooth)) {
    sm <- do.call(rbind, lapply(split(rows, rows$panel, drop = TRUE), function(d) {
      if (length(unique(d$dose)) < 3L) return(NULL)
      l <- stats::lowess(d$dose, d$residual)
      data.frame(panel = d$panel[1L], dose = l$x, residual = l$y)
    }))
    if (!is.null(sm)) {
      p <- p + ggplot2::geom_line(data = sm, linewidth = 0.9, colour = "black")
    }
  }
  # Up to three across; four sit as two by two.
  k <- nlevels(rows$panel)
  ncol <- if (!many) 1L else if (k <= 3L) k else min(3L, ceiling(k / 2))
  if (many) p <- p + ggplot2::facet_wrap(~ panel, ncol = ncol)

  spec <- dr_file_spec(file.name, filetype, path, width, height,
                       what = "easydr_residual()")
  nrow <- ceiling(nlevels(rows$panel) / ncol)
  title_in <- if (is.null(main)) 0 else 1.9 * 12 / 72
  size <- dr_fig_size(if (many) 3.4 * ncol + 0.6 else 6.5,
                      3.1 * nrow + 0.6 + title_in)
  con <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)
  dr_draw_stack(p, NULL, title = main, title_hjust = 0.5)

  invisible(rows[, c("shape", "studlab", "dose", "residual")])
}

#' Decorrelated residuals of one fit, one row per non-reference dose
#'
#' dosresmeta's `gof()` drops the reference rows, so its residuals line up
#' with the non-reference rows of the fit's frame, in the frame's order.
#'
#' @keywords internal
#' @noRd
dr_residual_rows <- function(f, shape) {
  g <- try(dosresmeta::gof(f$fit), silent = TRUE)
  if (inherits(g, "try-error") || is.null(g$tdata$tresiduals)) return(NULL)
  fr <- f$frame[!f$frame$is_ref, , drop = FALSE]
  res <- as.numeric(g$tdata$tresiduals)
  if (length(res) != nrow(fr)) {
    warning("easydr: the residuals of the ", shape, " shape do not line up ",
            "with its doses and are left out.", call. = FALSE)
    return(NULL)
  }
  data.frame(shape = shape, studlab = fr$studlab, dose = fr$dose,
             residual = res, stringsAsFactors = FALSE)
}

#' The name of the data a fit was called on, for a title
#'
#' `ma$Larsson_2011` reads as "Larsson_2011": the sheet name is what the
#' reader knows the data by.
#'
#' @keywords internal
#' @noRd
dr_data_name <- function(x) {
  d <- if (is.call(x$call)) x$call$data else NULL
  if (is.null(d)) return("the data")
  nm <- paste(deparse(d), collapse = "")
  sub("^.*\\$", "", gsub("[`\"]", "", nm))
}
