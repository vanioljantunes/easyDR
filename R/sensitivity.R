# Leave-one-out and cumulative analyses. Both ask the same thing of the same
# machinery: refit the inherited model on a sequence of subsets, and report how
# far the answer moved. They differ only in which subsets.

#' Which quantity the rows are compared on
#'
#' A line has a slope, so the slope is the thing to compare. A curve has no
#' single slope, so there is nothing to put in that column, and inventing one
#' would be worse than naming a dose. For more than one term the comparison is
#' the predicted effect at a stated dose, and the dose is reported rather than
#' assumed.
#'
#' @keywords internal
#' @noRd
dr_sens_quantity <- function(fit, at) {
  back <- if (fit$ratio) exp else identity
  if (!is.null(fit$slope)) {
    return(list(kind = "slope", dose = NA_real_,
                raw = fit$slope$estimate,
                est = back(fit$slope$estimate),
                lb  = back(fit$slope$lb),
                ub  = back(fit$slope$ub)))
  }
  row <- fit$predicted[!fit$predicted$is_ref & fit$predicted$dose == at, ,
                       drop = FALSE]
  if (nrow(row) == 0L) {
    stop("easydr: the comparison dose ", at, " is not among the predicted ",
         "doses of this refit.", call. = FALSE)
  }
  list(kind = "predicted", dose = at,
       raw = if (fit$ratio) log(row$est[1]) else row$est[1],
       est = row$est[1], lb = row$ci.lb[1], ub = row$ci.ub[1])
}

#' Pick the comparison dose for a curve
#'
#' The fitted result's own reported doses come first, because those are the
#' doses the reviewer already chose to look at. Failing that, the middle of the
#' curve grid, which is the least arbitrary single dose available.
#'
#' @keywords internal
#' @noRd
dr_sens_dose <- function(fit, at) {
  if (!is.null(fit$slope)) return(NA_real_)
  if (!is.null(at)) return(as.numeric(at)[1])
  own <- fit$predicted$dose[!fit$predicted$is_ref]
  if (length(own) > 0L) return(own[1])
  stats::median(fit$curve$dose[!fit$curve$is_ref])
}

#' Refit the inherited model once per subset and tabulate the movement
#'
#' A refit that fails records its reason on its own row and the loop carries on
#' (rule V24). The alternative, stopping, would throw away every row that did
#' work because of one study the shape cannot survive without.
#'
#' @keywords internal
#' @noRd
dr_sens_loop <- function(inh, subsets, labels, full, at, quiet, kind) {
  dose <- dr_sens_dose(full, at)
  base <- dr_sens_quantity(full, dose)
  n    <- length(subsets)
  rows <- vector("list", n)
  fits <- vector("list", n)
  t0   <- Sys.time()

  for (i in seq_len(n)) {
    keep <- subsets[[i]]
    if (!quiet) {
      cat(sprintf("  [%d/%d] %s\n", i, n, labels[i]))
      utils::flush.console()
    }
    extra <- if (is.na(dose)) NULL else list(at = dose)
    res <- tryCatch(dr_refit(inh, keep, extra), error = function(e) e)
    if (inherits(res, "error")) {
      rows[[i]] <- dr_sens_row(labels[i], length(keep), NULL, base,
                               conditionMessage(res))
      next
    }
    q <- tryCatch(dr_sens_quantity(res, dose), error = function(e) e)
    if (inherits(q, "error")) {
      rows[[i]] <- dr_sens_row(labels[i], length(keep), NULL, base,
                               conditionMessage(q))
      next
    }
    fits[[i]] <- res
    rows[[i]] <- dr_sens_row(labels[i], length(keep), list(q = q, fit = res),
                             base, NA_character_)
  }

  table <- do.call(rbind, rows)
  structure(
    list(
      kind       = kind,
      table      = table,
      quantity   = base$kind,
      dose       = dose,
      full       = base,
      ordered_on = if (identical(kind, "loo")) {
        "the absolute shift in the comparison quantity"
      } else {
        "the order of accumulation"
      },
      refits     = stats::setNames(fits, labels),
      inherited  = inh,
      sm         = full$sm,
      sm_label   = full$sm_label,
      ratio      = full$ratio,
      dose_unit  = full$dose_unit,
      level      = full$provenance$level,
      tried      = n,
      fitted     = sum(table$status == "fitted"),
      failed     = sum(table$status == "failed"),
      elapsed    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ),
    class = "easydr_sensitivity"
  )
}

#' One row of a sensitivity table
#' @keywords internal
#' @noRd
dr_sens_row <- function(label, studies, got, base, reason) {
  if (is.null(got)) {
    return(data.frame(
      label = label, studies = studies, estimate = NA_real_,
      lb = NA_real_, ub = NA_real_, estimate_raw = NA_real_,
      shift = NA_real_, abs_shift = NA_real_,
      Q = NA_real_, df = NA_integer_, p = NA_real_, i2 = NA_real_,
      status = "failed", reason = reason, stringsAsFactors = FALSE
    ))
  }
  q   <- got$q
  het <- got$fit$heterogeneity
  ok  <- isTRUE(het$estimable)
  data.frame(
    label = label, studies = studies, estimate = q$est,
    lb = q$lb, ub = q$ub, estimate_raw = q$raw,
    shift = q$raw - base$raw, abs_shift = abs(q$raw - base$raw),
    Q  = if (ok) het$Q else NA_real_,
    df = if (ok) as.integer(het$df) else NA_integer_,
    p  = if (ok) het$p else NA_real_,
    i2 = if (ok) het$i2 else NA_real_,
    status = "fitted", reason = NA_character_, stringsAsFactors = FALSE
  )
}

#' Leave one study out at a time
#'
#' Refits the model once without each study and reports how far the pooled
#' result moved, largest movement first. It describes which study is carrying
#' the result; it does not test anything, and no p value is attached to a
#' shift. Everything about the model is inherited from the fit you pass.
#'
#' @section What this answers:
#' A pooled dose-response is a weighted summary of studies that disagree. The
#' question a referee asks is whether one of them is carrying the result. This
#' refits the model once without each study and reports how far the pooled
#' answer moved each time, largest movement first.
#'
#' @section What it does not license:
#' This is a description, not a hypothesis test. No p value is attached to a
#' shift, because none belongs there: the refits are not independent of each
#' other or of the full fit, and whichever study moves the result most will
#' always be first, even when nothing is wrong. A large shift is a reason to go
#' and read that study, not a finding about it. Dropping a study because its
#' removal changes the answer turns a protocol into a conclusion.
#'
#' @section What is inherited:
#' Everything about the model comes from the fit you pass: the shape with its
#' knots or powers, the summary measure, the pooling stage, the covariance
#' approximation, the reference dose and the confidence level. The printed
#' output repeats them, so you can see they were inherited rather than trust
#' it. A refit the shape no longer suits is reported as a failed row with its
#' reason, never quietly fitted with a simpler shape.
#'
#' @param object a fitted result from [easydr()].
#' @param at for a shape with more than one term, the dose at which predicted
#'   effects are compared, since a curve has no single slope. Defaults to the
#'   first dose the fit already reports. Ignored for a one-term shape, which
#'   compares slopes.
#' @param quiet `FALSE`, the default, prints one line per refit.
#' @param plot `TRUE`, the default, draws the dose-response once per omitted
#'   study over the pooled one, letters each curve, and puts the table
#'   underneath. It then returns the result invisibly, so `print()` still gives
#'   the table when it is asked for. `FALSE` draws nothing.
#' @param main,xlab,ylab,ylim passed to that plot. `main = FALSE` drops the
#'   title; the rest default to the fit's own labels.
#'
#' @return An object of class `easydr_sensitivity`: the ordered table in
#'   `$table`, the refits themselves in `$refits` so any row can be inspected
#'   without repeating the work, and the inherited settings in `$inherited`.
#'
#' @seealso [easydr_forest()] for what each study contributed before any of
#'   them was left out.
#'
#' @examples
#' easydr.object <- easydr(alcohol_cvd, sm = "RR", shape = "linear",
#'                         dose.unit = "g/day")
#'
#' # One cohort is carrying this result, and it is the first row. The curve
#' # lettered A is the fit without it.
#' easydr_loo(easydr.object, quiet = TRUE)
#'
#' # The table on its own, with nothing drawn.
#' loo <- easydr_loo(easydr.object, quiet = TRUE, plot = FALSE)
#' loo$table[, c("label", "estimate", "shift", "i2")]
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
easydr_loo <- function(object, at = NULL, quiet = FALSE, plot = TRUE,
                       main = NULL, xlab = NULL, ylab = NULL, ylim = NULL,
                       filetype = c("png", "pdf", "tiff"),
                       path = NULL, file.name = NULL,
                       width = NULL, height = NULL) {
  dr_require_single(object, "easydr_loo()")
  inh <- dr_inherited(object)
  if (length(inh$studies) < 2L) {
    stop("easydr: leave-one-out needs at least two studies; this fit has ",
         length(inh$studies), ".", call. = FALSE)
  }
  subsets <- lapply(inh$studies, function(s) setdiff(inh$studies, s))
  out <- dr_sens_loop(inh, subsets, inh$studies, object, at, quiet,
                      kind = "loo")
  ord <- order(out$table$abs_shift, decreasing = TRUE, na.last = TRUE)
  out$table  <- out$table[ord, , drop = FALSE]
  out$refits <- out$refits[ord]
  rownames(out$table) <- NULL
  out$call <- match.call()
  # The picture is the point of a leave-one-out: a table of six near-identical
  # numbers hides the one curve that moved. Drawing it returns invisibly, so
  # `print()` still gives the table when it is asked for.
  if (isTRUE(plot)) {
    # A figure that cannot be drawn must not cost the caller the analysis:
    # when every refit failed the table is the whole answer, and it says why.
    if (any(out$table$status == "fitted")) {
      out$plot <- dr_loo_plot(out, object, main = main, xlab = xlab,
                              ylab = ylab, ylim = ylim,
                              filetype = filetype, path = path,
                              file.name = file.name,
                              width = width, height = height)
      return(invisible(out))
    }
    warning("easydr: no refit succeeded, so nothing was drawn. The table ",
            "gives the reason against each row.", call. = FALSE)
  }
  out
}

#' Print a leave-one-out or cumulative analysis
#'
#' The block order of features 001 and 002: what the quantity is, what it is
#' measured against, what was inherited, then the table, then every default
#' that was applied rather than chosen.
#'
#' @param x an object from [easydr_loo()].
#' @param digits digits for the reported quantity.
#' @param ... ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.easydr_sensitivity <- function(x, digits = 4, ...) {
  unit <- if (is.null(x$dose_unit)) "dose unit" else x$dose_unit
  fmt  <- function(v) dr_fmt(v, digits, na = "not estimable")
  loo  <- identical(x$kind, "loo")

  cat(if (loo) "Leave-one-out sensitivity (" else "Cumulative analysis (",
      x$sm_label, ")\n", sep = "")
  cat(x$tried, if (loo) " refits, " else " steps, ", x$fitted, " fitted, ",
      x$failed, " failed, ", format(round(x$elapsed, 1), nsmall = 1),
      " seconds\n", sep = "")
  cat("Inherited: ", dr_inherited_line(x$inherited), "\n", sep = "")

  quantity <- if (identical(x$quantity, "slope")) {
    paste0("the slope, ", x$sm, " per 1 ", unit)
  } else {
    paste0("the predicted ", x$sm, " at ", x$dose, " ", unit)
  }
  if (loo) {
    cat("Ordered by the absolute shift in ", quantity, ". Full fit: ",
        fmt(x$full$est), "\n\n", sep = "")
  } else {
    stated <- if (isTRUE(x$order_stated)) {
      "as stated"
    } else {
      "as the studies appear in the data (no order was stated)"
    }
    cat("Comparing ", quantity, ". Full fit: ", fmt(x$full$est), "\n", sep = "")
    cat("Order: ", stated, "\n\n", sep = "")
  }

  tab <- x$table
  body <- data.frame(
    a = if (loo) tab$label else paste0(seq_len(nrow(tab)), "  ", tab$label),
    studies = tab$studies,
    est = fmt(tab$estimate),
    ci = ifelse(is.na(tab$lb), "not estimable",
                paste0(fmt(tab$lb), " to ", fmt(tab$ub))),
    shift = ifelse(is.na(tab$shift), "",
                   sprintf("%+.*f", digits, tab$shift)),
    i2 = ifelse(is.na(tab$i2), "not estimable",
                paste0(format(round(tab$i2, 1), nsmall = 1), "%")),
    stringsAsFactors = FALSE
  )
  names(body) <- c(if (loo) "omitted" else "step and study", "n", x$sm,
                   paste0(format(100 * x$level), "% CI"), "shift", "I-squared")
  print(body, row.names = FALSE, right = FALSE)

  failed <- tab[tab$status == "failed", , drop = FALSE]
  if (nrow(failed) > 0L) {
    cat("\nRefits that could not be fitted:\n")
    for (i in seq_len(nrow(failed))) {
      cat("  - ", failed$label[i], ": ", failed$reason[i], "\n", sep = "")
    }
  }

  if (!loo) {
    cat("\nThe final step reproduces the full fit, which is the check that ",
        "the sequence is coherent.\n", sep = "")
  }

  cat("\nDefaults applied:\n")
  for (d in dr_sens_defaults(x)) cat("  - ", d, "\n", sep = "")
  invisible(x)
}

#' Every default a sensitivity analysis applied on the reviewer's behalf
#' @keywords internal
#' @noRd
dr_sens_defaults <- function(x) {
  unit <- if (is.null(x$dose_unit)) "dose unit" else x$dose_unit
  out <- if (identical(x$quantity, "slope")) {
    "comparison quantity: the slope, because this shape has one term"
  } else {
    paste0("comparison quantity: the predicted effect at ", x$dose, " ", unit,
           ", because a curve has no single slope")
  }
  scale <- if (x$ratio) "the log scale" else "the reported scale"
  if (identical(x$kind, "loo")) {
    out <- c(out, paste0("the shift is on ", scale,
                         ", which the model is fitted on, and the rows are ",
                         "ordered on that same quantity"))
    out <- c(out, paste0(
      "a shift carries no p value: the refits are not independent, so this ",
      "describes the data rather than testing it"))
  } else {
    out <- c(out, paste0("the shift is on ", scale,
                         ", which the model is fitted on, and is measured ",
                         "against the full fit rather than the previous step"))
    out <- c(out, paste0(
      "the steps share almost all their data, so consecutive estimates are ",
      "dependent and the sequence is not a trend"))
  }
  out
}

# Cumulative pooling. Kept as an internal worker: the sequence is easy to
# misread as a trend over time, and it is not part of the public family yet.
# Exposing it is a decision about the published surface, not about the code.

#' Add the studies up in order
#'
#' Refits the model on the first study, then the first two, and so on. The
#' steps share almost all their data, so the sequence is not a trend and not a
#' test of accumulating evidence. The order is the caller's to state; the
#' output says which one was used.
#'
#' @keywords internal
#' @noRd
dr_cumul <- function(object, order = NULL, at = NULL, quiet = FALSE) {
  inh     <- dr_inherited(object)
  stated  <- !is.null(order)
  ord     <- if (stated) as.character(order) else inh$studies
  dr_check_order(ord, inh$studies)

  subsets <- lapply(seq_along(ord), function(i) ord[seq_len(i)])
  out <- dr_sens_loop(inh, subsets, ord, object, at, quiet,
                      kind = "cumulative")
  out$order        <- ord
  out$order_stated <- stated
  out$call         <- match.call()
  out
}

#' Rule V22: a stated cumulative order must match the studies in the data
#'
#' Stops rather than accumulating a sequence the reviewer did not describe. A
#' quietly reordered or silently shortened sequence would still print one row
#' per step and still look right.
#'
#' @keywords internal
#' @noRd
dr_check_order <- function(ord, studies) {
  if (anyDuplicated(ord)) {
    stop("easydr: the stated order repeats ",
         paste(unique(ord[duplicated(ord)]), collapse = ", "),
         ". Each study appears once in a cumulative sequence.", call. = FALSE)
  }
  unknown <- setdiff(ord, studies)
  missing <- setdiff(studies, ord)
  if (length(unknown) == 0L && length(missing) == 0L) return(invisible(TRUE))
  parts <- character()
  if (length(unknown) > 0L) {
    parts <- c(parts, paste0("not in the data: ",
                             paste(unknown, collapse = ", ")))
  }
  if (length(missing) > 0L) {
    parts <- c(parts, paste0("in the data but not in the order: ",
                             paste(missing, collapse = ", ")))
  }
  stop("easydr: the stated order does not match the studies in the fit. ",
       paste(parts, collapse = "; "), ".", call. = FALSE)
}
