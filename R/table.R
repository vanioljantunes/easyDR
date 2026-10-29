# The table a dose-response paper puts next to its curve.
#
# A curve answers "what is the shape". A reader planning a dose asks "what is
# the effect at 10, and at 25, and at 50". That is a different object, it is
# the one that goes in Table 2 of a manuscript, and it was previously built by
# hand out of `fit$predicted`. This builds it from the fit, at whatever doses
# the reviewer names, and writes it out.

#' The effect at the doses you name
#'
#' Takes a fitted dose-response model and reports the summary measure, its
#' confidence interval and its prediction interval at any doses asked for. The
#' reference dose is always included, because an effect is a comparison and the
#' table has to say what it is compared against.
#'
#' This is the table a manuscript reports beside its curve. The doses belong in
#' the protocol rather than in the output: choosing them after seeing the shape
#' is choosing the result. They are usually the intakes existing guidance
#' already speaks about.
#'
#' A dose outside the range the studies reported is still estimated, because
#' the model is a function and a function has a value everywhere. It is marked
#' as extrapolated, and an extrapolated row should not be reported as a
#' finding.
#'
#' @param x a fitted result from [easydr()].
#' @param doses doses to report, such as `doses = c(10, 25, 50)`. `NULL` uses
#'   the doses the fit already reports.
#' @param main a title. It heads a PDF table and names the sheet of an
#'   `.xlsx` file. `NULL` writes a title from the summary measure and the
#'   reference dose.
#' @param digits digits for the reported effects.
#' @param filetype `"xlsx"` or `"pdf"`. Only used when `file.name` is given.
#' @param path the directory to write into. `NULL` is the working directory.
#' @param file.name the name of the file to write. `NULL`, the default,
#'   writes nothing and returns the table.
#' @param width,height size of a PDF table in inches. `NULL` lets easyDR work
#'   it out from the table itself.
#'
#' @return A data frame of one row per dose, with the summary measure, its
#'   confidence limits, its prediction limits and whether the row is
#'   extrapolated. Returned invisibly when a file was written.
#'
#' @examples
#' data(alcohol_cvd)
#' easydr.object <- easydr(alcohol_cvd, sm = "OR", shape = "rcs",
#'                         dose.unit = "g/day")
#'
#' # The effect at the intakes guidance speaks about.
#' easydr_table(easydr.object, doses = c(10, 25, 50))
#'
#' \donttest{
#' # The same table, written to a spreadsheet in the working directory.
#' easydr_table(easydr.object, doses = c(10, 25, 50),
#'              main = "Alcohol and cardiovascular disease",
#'              filetype = "xlsx", file.name = "table2")
#' }
#'
#' @export
easydr_table <- function(x, doses = NULL, main = NULL, digits = 2,
                         filetype = c("xlsx", "pdf"), path = NULL,
                         file.name = NULL, width = NULL, height = NULL) {
  if (!inherits(x, "easydr_result")) {
    stop("easydr: easydr_table() takes a fitted result from easydr().",
         call. = FALSE)
  }
  xref <- x$reference$curve
  if (is.null(doses)) {
    doses <- x$predicted$dose
  }
  if (!is.numeric(doses) || !length(doses) || anyNA(doses) ||
      !all(is.finite(doses))) {
    stop("easydr: `doses` must be a numeric vector, such as ",
         "doses = c(10, 25, 50).", call. = FALSE)
  }
  # The reference joins the grid whether it was asked for or not. Every number
  # in the table is measured against it, so a table without it is a column of
  # effects with nothing to compare them to.
  grid <- sort(unique(c(xref, as.numeric(doses))))

  pred <- dr_predicted(x$fit, grid, xref, x$provenance$level, x$sm,
                       x$provenance$single_study, x$shape,
                       range(x$frame$dose))$table

  lvl  <- format(100 * x$provenance$level)
  unit <- if (is.null(x$dose_unit)) NULL else x$dose_unit
  dose_col <- if (is.null(unit)) "Dose" else paste0("Dose (", unit, ")")

  out <- data.frame(pred$dose, round(pred$est, digits),
                    round(pred$ci.lb, digits), round(pred$ci.ub, digits),
                    round(pred$pi.lb, digits), round(pred$pi.ub, digits),
                    pred$extrapolated,
                    check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- c(dose_col, x$sm,
                  paste0(lvl, "% CI low"),  paste0(lvl, "% CI high"),
                  paste0(lvl, "% PI low"),  paste0(lvl, "% PI high"),
                  "Extrapolated")
  rownames(out) <- NULL
  # A column of blanks tells the reader nothing. It appears when there is
  # something in it to read.
  if (!any(pred$extrapolated)) out[["Extrapolated"]] <- NULL

  if (is.null(main)) {
    main <- paste0(x$sm_label, " against ", xref,
                   if (!is.null(unit)) paste0(" ", unit))
  }
  attr(out, "main") <- main
  attr(out, "reference") <- xref

  if (is.null(file.name)) return(out)

  filetype <- match.arg(tolower(filetype), c("xlsx", "pdf"))
  if (identical(filetype, "xlsx")) {
    dr_table_xlsx(out, main, path, file.name)
  } else {
    dr_table_pdf(out, main, path, file.name, width, height)
  }
  invisible(out)
}

#' Write a table to a spreadsheet
#'
#' writexl is a one-function package with no Java and no system library behind
#' it, which is why it is the one asked for. It is suggested rather than
#' imported: a reviewer who never exports a table should not have to install
#' it.
#'
#' @keywords internal
#' @noRd
dr_table_xlsx <- function(df, main, path, file.name) {
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("easydr: writing an .xlsx file needs the writexl package. ",
         "Install it with install.packages(\"writexl\"), or use ",
         "filetype = \"pdf\".", call. = FALSE)
  }
  if (is.null(path)) path <- getwd()
  if (!dir.exists(path)) {
    stop("easydr: the directory `path` names does not exist: ", path,
         call. = FALSE)
  }
  if (!grepl("\\.xlsx$", file.name, ignore.case = TRUE)) {
    file.name <- paste0(file.name, ".xlsx")
  }
  file <- file.path(path, file.name)
  # A sheet name is capped at 31 characters and cannot hold \ / ? * [ ] :
  sheet <- gsub("[\\\\/?*\\[\\]:]", " ", main)
  sheet <- substr(sheet, 1L, 31L)
  if (!nzchar(trimws(sheet))) sheet <- "easyDR"
  plain <- as.data.frame(df)
  attr(plain, "main") <- NULL
  attr(plain, "reference") <- NULL
  writexl::write_xlsx(stats::setNames(list(plain), sheet), path = file)
  message("easydr: wrote ", file)
  invisible(file)
}

#' Draw a table into a PDF, at the size the table needs
#' @keywords internal
#' @noRd
dr_table_pdf <- function(df, main, path, file.name, width, height) {
  spec <- dr_file_spec(file.name, "pdf", path, width, height,
                       what = "easydr_table()")
  show <- df
  if (!is.null(df[["Extrapolated"]])) {
    show[["Extrapolated"]] <- ifelse(df[["Extrapolated"]], "yes", "")
  }
  tg <- dr_stats_grob(show)
  title_in <- if (is.null(main)) 0 else 1.9 * 12 / 72
  size <- dr_fig_size(dr_grob_width_in(tg) + 0.6,
                      dr_height_in(dr_grob_height(tg)) + title_in + 0.3,
                      min.w = 3, min.h = 1.2)
  con <- dr_begin_file(spec, size[["width"]], size[["height"]])
  on.exit(dr_end_file(con), add = TRUE)
  dr_draw_stack(tg, dr_grob_height(tg), title = main)
  invisible(spec$file)
}
