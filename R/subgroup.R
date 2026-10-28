# Subgroups: the same model, fitted again inside each level of one column.
#
# A subgroup is not a new model. It is the model already fitted, run on fewer
# studies, so every setting travels through dr_refit(): the shape as a resolved
# object with its fitted knots, the pooling stage, the covariance
# approximation, the columns by the name the first fit resolved. A level that
# re-derived any of those would print cleanly and answer a different question.
#
# Each level keeps its own tau, the way meta's subgroup analysis does. There is
# no test for subgroup differences: the table reports the levels side by side
# and leaves the comparison to the reader.

#' Fit the inherited model inside each level of a subgroup column
#'
#' @param object a fitted result.
#' @param var the name of the column holding the subgroup value.
#'
#' @return A list with `var`, the `levels` in order, `by_study` naming each
#'   study's level, the per-level `fits` (`NULL` where the refit failed), and a
#'   `table` of one row per level carrying the study count, the status and the
#'   reason for any failure.
#'
#' @keywords internal
#' @noRd
dr_subgroup_fits <- function(object, var) {
  if (!is.character(var) || length(var) != 1L || is.na(var)) {
    stop("easydr: `subgroup` is the name of one column, given as a string, ",
         "for example subgroup = \"design\".", call. = FALSE)
  }
  inh <- dr_inherited(object)

  # The same rule moderators follow: a subgroup is a characteristic of a
  # study, so it takes one value per study, and a column with one value across
  # every study splits nothing. dr_mod_study_table() refuses both, and refuses
  # a column that is not in the data.
  tab <- dr_mod_study_table(inh$data, inh$studlab, var, what = "subgroup")
  value <- as.character(tab[[var]])
  levels_in_order <- unique(value[order(value)])

  fits <- vector("list", length(levels_in_order))
  names(fits) <- levels_in_order
  rows <- vector("list", length(levels_in_order))

  # Every level reports the same doses and is evaluated on the same grid as
  # the overall fit. Left to its own defaults a level picks doses from its own
  # range, and the figure then marks each curve at a different place, which
  # reads as a difference between the levels rather than between the ranges.
  same_doses <- list(at = object$predicted$dose[!object$predicted$is_ref],
                     grid = object$curve$dose)

  for (i in seq_along(levels_in_order)) {
    lv   <- levels_in_order[i]
    keep <- tab$.study[value == lv]
    res  <- tryCatch(dr_refit(inh, keep, extra = same_doses),
                     error = function(e) e)

    if (inherits(res, "error")) {
      # A level too small or too sparse for the inherited model is a fact
      # about the data, recorded against its row. It does not stop the levels
      # that did fit from being reported.
      rows[[i]] <- data.frame(
        level = lv, studies = length(keep), status = "failed",
        reason = dr_one_line(conditionMessage(res)),
        stringsAsFactors = FALSE)
    } else {
      fits[[i]] <- res
      rows[[i]] <- data.frame(
        level = lv, studies = length(keep), status = "fitted",
        reason = NA_character_, stringsAsFactors = FALSE)
    }
  }

  out <- list(var = var, levels = levels_in_order,
              by_study = stats::setNames(value, tab$.study),
              fits = fits, table = do.call(rbind, rows))

  bad <- out$table$status == "failed"
  if (any(bad)) {
    warning("easydr: the model could not be fitted in ", sum(bad), " of ",
            nrow(out$table), " subgroup levels (",
            paste(out$table$level[bad], collapse = ", "),
            "). The subgroup table gives the reason against each row.",
            call. = FALSE)
  }
  out
}

#' An error message as one line, short enough for a table cell
#' @keywords internal
#' @noRd
dr_one_line <- function(msg) {
  msg <- gsub("[\r\n]+", " ", msg)
  msg <- gsub("^easydr: ", "", msg)
  msg <- gsub("[[:space:]]+", " ", msg)
  if (nchar(msg) > 90L) msg <- paste0(substr(msg, 1L, 87L), "...")
  msg
}

#' Does this result carry subgroups that were fitted
#' @keywords internal
#' @noRd
dr_has_subgroups <- function(x) {
  !is.null(x$subgroup) && any(x$subgroup$table$status == "fitted")
}

#' Resolve the `subgroup` argument of a plotting function
#'
#' `TRUE` with nothing to draw is a mistake worth naming: silently drawing the
#' overall fit would answer a question the caller did not ask.
#'
#' @keywords internal
#' @noRd
dr_use_subgroups <- function(x, subgroup, what) {
  if (isFALSE(subgroup)) return(FALSE)
  if (is.null(x$subgroup)) {
    if (isTRUE(subgroup)) {
      stop("easydr: this fit carries no subgroup, so ", what, " has nothing ",
           "to group by. Refit with easydr(..., subgroup = \"<column>\") or ",
           "call this with subgroup = FALSE.", call. = FALSE)
    }
    return(FALSE)
  }
  if (!dr_has_subgroups(x)) {
    stop("easydr: no subgroup level could be fitted, so ", what, " has ",
         "nothing to draw. The reason for each level is in ",
         "`fit$subgroup$table`.", call. = FALSE)
  }
  TRUE
}

#' The colour, line type and point shape each subgroup level is drawn with
#'
#' Colour alone fails in print and for a colour-blind reader, so a level is
#' identified twice over: by colour and by line type for the curves, by colour
#' and by shape for the marked points.
#'
#' @keywords internal
#' @noRd
dr_subgroup_style <- function(levels) {
  n <- length(levels)
  cols <- c("#1F2A44", "#C0392B", "#1B7F5A", "#8E44AD", "#B8860B",
            "#2E86C1", "#D35400")
  ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  # 22 square, 24 triangle, 23 diamond, 8 star, 25 down triangle, 21 circle.
  # The circle comes last because the overall fit's marked doses are circles,
  # and a level sharing that shape reads as the overall curve.
  shps <- c(22L, 24L, 23L, 8L, 25L, 21L)
  list(
    colour   = stats::setNames(rep_len(cols, n), levels),
    linetype = stats::setNames(rep_len(ltys, n), levels),
    shape    = stats::setNames(rep_len(shps, n), levels)
  )
}
