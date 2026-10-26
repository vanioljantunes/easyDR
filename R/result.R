# The result object and how it prints. One fitted object answers the review
# question, so nothing here asks the reviewer to call anything else.

#' Build the result object
#' @keywords internal
#' @noRd
dr_result <- function(fit, prep, sm, dose_unit, at, level, proc, method, call,
                      shape = NULL, grid = NULL, data = NULL, entry = NULL,
                      args = NULL) {
  if (is.null(shape)) shape <- dr_shape_linear()
  curved <- shape$terms > 1L
  frame <- prep$frame
  prov  <- prep$provenance
  ratio <- dr_is_ratio(sm)

  # The pooled curve needs one reference dose. Per-study references are honoured
  # inside the engine; this is the dose the printed curve is measured against.
  # When studies disagree, the lowest reference dose actually used is chosen, not
  # the lowest dose in the table: a dose no study treated as its reference would
  # make the printed contrasts answer a question nobody asked. The disagreement
  # is recorded and printed.
  ref_doses <- frame$dose[frame$is_ref]
  xref <- min(ref_doses)
  if (length(unique(ref_doses)) > 1L) {
    prov$warnings <- c(prov$warnings, paste0(
      "per-study reference doses differ, from ", min(ref_doses), " to ",
      max(ref_doses), "; predicted effects are reported against ", xref))
  }

  # A curve has no single slope, so the slot is NULL rather than removed: a
  # feature 001 reader of a linear fit keeps working (FR-011, data model).
  slope <- if (curved) {
    NULL
  } else {
    b     <- as.numeric(coef(fit))[1]
    se    <- sqrt(as.numeric(vcov(fit))[1])
    zcrit <- stats::qnorm(1 - (1 - level) / 2)
    list(label = dr_label(sm, dose_unit), estimate = b, se = se,
         lb = b - zcrit * se, ub = b + zcrit * se,
         p = 2 * stats::pnorm(-abs(b / se)),
         scale = if (ratio) "log" else "raw")
  }

  doses <- if (is.null(at)) {
    dr_default_doses(frame$dose, xref)
  } else {
    sort(unique(c(xref, at)))
  }
  pred <- dr_predicted(fit, doses, xref, level, sm, prov$single_study, shape,
                       range(frame$dose))
  het  <- dr_heterogeneity(fit, prov$single_study)

  # The curve: a denser grid than the reported doses, for reading off and for a
  # later plotting feature to consume unchanged.
  curve_doses <- if (is.null(grid)) {
    sort(unique(c(xref, seq(min(frame$dose), max(frame$dose), length.out = 20))))
  } else {
    sort(unique(c(xref, as.numeric(grid))))
  }
  curve <- dr_predicted(fit, curve_doses, xref, level, sm, prov$single_study,
                        shape, range(frame$dose))$table

  structure(
    list(
      fit      = fit,
      slope    = slope,
      sm       = sm,
      sm_label = dr_sm_label(sm),
      ratio    = ratio,
      dose_unit = dose_unit,
      reference = list(
        curve = xref,
        per_study = stats::setNames(ref_doses, frame$studlab[frame$is_ref])
      ),
      heterogeneity = het,
      Psi = fit$Psi,
      shape = shape,
      nonlinearity = dr_nonlinearity_test(fit, shape, prov),
      gof = dr_gof(fit),
      curve = curve,
      prediction = pred$pi,
      predicted  = pred$table,
      counts = list(
        studies    = length(prep$studies),
        categories = nrow(frame),
        dose_min   = min(frame$dose),
        dose_max   = max(frame$dose)
      ),
      provenance = c(prov, list(
        proc = proc,
        method = method,
        covariance = attr(fit, "easydr_covariance"),
        covariance_label = dr_covariance_label(attr(fit, "easydr_covariance")),
        level = level
      )),
      frame = frame,
      # The table as supplied, the entry point that was called and the column
      # and model arguments it was called with. A sensitivity analysis refits by
      # subsetting this table and calling that entry point again with these
      # arguments, so nothing about the model is re-derived from the subset.
      data  = data,
      entry = entry,
      args  = args,
      call  = call
    ),
    class = "easydr_result"
  )
}

#' Five doses spanning the observed range, with the reference dose included
#' @keywords internal
#' @noRd
dr_default_doses <- function(dose, xref) {
  rng <- range(dose)
  sort(unique(c(xref, signif(seq(rng[1], rng[2], length.out = 5), 4))))
}

#' Round and pad a number for printing
#'
#' Three print methods formatted numbers with the same two lines. `na` gives
#' the text to use where a value is not estimable, which only the sensitivity
#' tables need.
#'
#' @keywords internal
#' @noRd
dr_fmt <- function(v, digits, na = NULL) {
  out <- format(round(v, digits), nsmall = digits)
  if (is.null(na)) return(out)
  ifelse(is.na(v), na, out)
}

#' Predicted effects and the prediction interval, back-transformed
#'
#' The engine supplies the point estimate and the confidence limits. It requires
#' the reference dose to be present in `newdata` and takes the first row as the
#' reference unless told otherwise, so the reference dose is always included.
#'
#' @keywords internal
#' @noRd
dr_predicted <- function(fit, doses, xref, level, sm, single_study,
                         shape = NULL, observed = NULL) {
  if (is.null(shape)) shape <- dr_shape_linear()

  # The basis is rebuilt from the shape stored parameters, never re-derived
  # from this grid, or the model would change between fitting and prediction
  # (research decision D4).
  newdata <- if (identical(shape$name, "linear")) {
    data.frame(dose = doses)
  } else {
    as.data.frame(shape$basis(doses))
  }
  # The reference dose is always one of the predicted doses, so its basis value
  # is read from the grid rather than from a separate one-value basis call.
  xref_row <- if (identical(shape$name, "linear")) {
    xref
  } else {
    as.numeric(newdata[[1L]][match(xref, doses)])
  }

  ci <- stats::predict(fit, newdata = newdata, xref = xref_row,
                       expo = FALSE, ci.level = level)
  pi <- dr_prediction_interval(fit, doses, xref, level, single_study, shape)

  back <- if (dr_is_ratio(sm)) exp else identity
  tab <- data.frame(
    dose   = doses,
    est    = back(ci$pred),
    ci.lb  = back(ci$ci.lb),
    ci.ub  = back(ci$ci.ub),
    pi.lb  = if (isTRUE(pi$estimable)) back(pi$lb) else NA_real_,
    pi.ub  = if (isTRUE(pi$estimable)) back(pi$ub) else NA_real_,
    is_ref = doses == xref,
    extrapolated = if (is.null(observed)) {
      rep(FALSE, length(doses))
    } else {
      tol <- 1e-6 * max(1, abs(observed[2] - observed[1]))
      doses < observed[1] - tol | doses > observed[2] + tol
    },
    stringsAsFactors = FALSE
  )
  list(table = tab, pi = pi)
}

#' Print a dose-response meta-analysis
#'
#' Prints the whole answer: the pooled slope on the reported scale, the dose unit
#' and reference dose it is measured against, the departure-from-linearity
#' status, the engine's heterogeneity test, the prediction interval, the
#' predicted effects, and every default easydr applied on the reviewer's behalf.
#'
#' @param x an object from [easydr()].
#' @param digits significant digits for reported effects.
#' @param ... ignored.
#'
#' @return `x`, invisibly. Called for the printed output.
#'
#' @examples
#' fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear", dose.unit = "g/day",
#'                  at = c(10, 25, 50))
#' print(fit)
#'
#' # Every default easydr applied is printed under "Defaults applied", so the
#' # methods paragraph of a manuscript can be written from this output alone.
#'
#' @export
print.easydr_result <- function(x, digits = 3, ...) {
  unit <- if (is.null(x$dose_unit)) "dose unit" else x$dose_unit
  p    <- x$provenance
  back <- if (x$ratio) exp else identity
  fmt  <- function(v) dr_fmt(v, digits)

  cat("Dose-response meta-analysis (", x$sm_label, ")\n", sep = "")
  cat("Pooling: ", if (identical(p$proc, "2stage")) "two-stage" else "one-stage",
      # A common-effect fit estimates no heterogeneity parameter, so naming a
      # tau estimator beside it would describe something that was not done.
      if (identical(p$method, "fixed")) {
        " common effect"
      } else {
        paste0(" random effects, ", p$method)
      },
      ". Covariance: ", p$covariance_label,
      ".\n\n", sep = "")

  if (is.null(x$slope)) {
    # A curve has no single slope, and saying so is the point (FR-011).
    cat("Shape: ", x$shape$label, "
", sep = "")
    cat("No single slope describes this shape. Read the predicted effects ",
        "below.
", sep = "")
  } else {
    # The slope, labelled by easydr and never by the engine.
    cat(x$slope$label, ":  ", x$sm, " ", fmt(back(x$slope$estimate)),
        "  (", format(100 * p$level), "% CI ", fmt(back(x$slope$lb)), " to ",
        fmt(back(x$slope$ub)), ")
", sep = "")
  }
  cat("Reference dose: ", x$reference$curve, " ", unit, "\n", sep = "")
  # Only easydr() names the curve; the internal workers leave the slot empty,
  # so a refit inside a sensitivity analysis prints exactly what it always did.
  if (!is.null(x$curve_label) && !is.na(x$curve_label)) {
    cat("Curve shape: ", x$curve_label, "\n", sep = "")
  }
  cat("\n")

  if (isTRUE(x$nonlinearity$estimable)) {
    cat("Departure from linearity: X2 = ", round(x$nonlinearity$statistic, 2),
        " (df = ", x$nonlinearity$df, "), p = ",
        format.pval(x$nonlinearity$p, digits = 3), "
", sep = "")
    cat("  ", x$nonlinearity$comparison, "
", sep = "")
  } else {
    cat("Departure from linearity: not estimable
")
    cat("  ", x$nonlinearity$reason, "
", sep = "")
    if (!is.null(x$nonlinearity$studies)) {
      cat("  ", x$nonlinearity$studies_with_three_levels, " of ",
          x$nonlinearity$studies, " studies have three or more dose levels
",
          sep = "")
    }
  }
  if (isTRUE(x$gof$available) && isTRUE(x$nonlinearity$estimable)) {
    cat("  Goodness of fit: deviance D = ", round(x$gof$deviance, 2),
        " (df = ", x$gof$df, "), p = ", format.pval(x$gof$p, digits = 3),
        "; R-squared ", round(x$gof$r2, 3), "
", sep = "")
  }
  cat("
")

  if (isTRUE(x$heterogeneity$estimable)) {
    cat("Heterogeneity: Q = ", round(x$heterogeneity$Q, 2), " (df = ",
        x$heterogeneity$df, "), p = ",
        format.pval(x$heterogeneity$p, digits = 3), "; I-squared = ",
        round(x$heterogeneity$i2, 1), "%\n", sep = "")
  } else {
    cat("Heterogeneity: not estimable (", x$heterogeneity$reason, ")\n", sep = "")
  }
  if (!isTRUE(x$prediction$estimable)) {
    cat("Prediction interval: not estimable (", x$prediction$reason, ")\n",
        sep = "")
  }
  cat("\n")

  cat("Predicted ", x$sm, "\n", sep = "")
  tab <- x$predicted
  lines <- data.frame(
    dose = format(tab$dose),
    est  = ifelse(tab$is_ref, "reference", fmt(tab$est)),
    ci   = ifelse(tab$is_ref, "reference",
                  paste0(fmt(tab$ci.lb), " to ", fmt(tab$ci.ub))),
    pi   = if (isTRUE(x$prediction$estimable)) {
      ifelse(tab$is_ref, "reference",
             paste0(fmt(tab$pi.lb), " to ", fmt(tab$pi.ub)))
    } else {
      rep("not estimable", nrow(tab))
    },
    stringsAsFactors = FALSE
  )
  names(lines) <- c(paste0("dose (", unit, ")"), x$sm,
                    paste0(format(100 * p$level), "% CI"),
                    paste0(format(100 * p$level), "% PI"))
  if (any(tab$extrapolated)) {
    lines$note <- ifelse(tab$extrapolated, "beyond observed doses", "")
  }
  print(lines, row.names = FALSE, right = FALSE)
  cat("\n")

  cat(x$counts$studies, " studies, ", x$counts$categories,
      " dose categories, dose range ", x$counts$dose_min, " to ",
      x$counts$dose_max, " ", unit, "\n", sep = "")

  if (!is.null(x$viability)) cat(dr_viability_line(x$viability), "\n", sep = "")

  defaults <- dr_defaults_applied(x)
  if (length(defaults) > 0L) {
    cat("Defaults applied:\n")
    for (d in defaults) cat("  - ", d, "\n", sep = "")
  }
  invisible(x)
}

#' List every default and warning the reviewer should know about
#' @keywords internal
#' @noRd
dr_defaults_applied <- function(x) {
  p <- x$provenance
  out <- paste0("intake: ", switch(p$intake,
    counts        = "event counts, with the effect computed from them",
    effect        = "reported effect with limits, no denominators supplied",
    effect_counts = "reported effect with limits, counts used for the covariance",
    continuous    = "means with standard deviations",
    p$intake))
  if (length(p$inferred_ref) > 0L) {
    out <- c(out, paste0("reference category inferred in ",
                         length(p$inferred_ref), " stud",
                         if (length(p$inferred_ref) == 1L) "y" else "ies",
                         ": ", paste(p$inferred_ref, collapse = ", ")))
  }
  if (length(p$continuity) > 0L) {
    out <- c(out, paste0("zero-cell correction applied in: ",
                         paste(p$continuity, collapse = ", ")))
  }
  out <- c(out, dr_mixed_designs_note(p))
  out <- c(out, setdiff(p$warnings, ""))
  unique(out)
}

#' Say so when the studies do not share one design
#'
#' Each study keeps its own design, which is what selects its covariance
#' approximation. When they differ, the pooled ratio is a mixture of the
#' measures the authors reported, so the print says which designs are present
#' rather than letting one `sm` imply they all share it.
#'
#' @keywords internal
#' @noRd
dr_mixed_designs_note <- function(p) {
  d <- p$designs
  if (length(d) == 0L) return(character())
  kinds <- unique(unname(d))
  if (length(kinds) < 2L) return(character())
  words <- c(cc = "case-control", ir = "incidence rate",
             ci = "cumulative incidence")
  counts <- vapply(kinds, function(k) sum(d == k), integer(1))
  parts  <- paste0(counts, " ", ifelse(is.na(words[kinds]), kinds, words[kinds]))
  paste0("designs differ across studies (", paste(parts, collapse = ", "),
         "); each study keeps its own covariance approximation and the pooled ",
         "ratio is the mixture of the measures reported")
}

#' Summarise a dose-response meta-analysis
#'
#' Adds the per-study contribution table to everything [print.easydr_result()]
#' shows: how many dose levels each study brought, its dose range, its design,
#' and its own slope with a standard error.
#'
#' @param object an object from [easydr()].
#' @param ... ignored.
#'
#' @return An object of class `summary.easydr_result`.
#'
#' @examples
#' # The per-study table shows how many dose levels each study brought and what
#' # slope it contributed, which is what a sensitivity discussion needs.
#' summary(easydr(alcohol_cancer, sm = "IRR", shape = "linear", dose.unit = "g/day"))
#'
#' @export
summary.easydr_result <- function(object, ...) {
  frame   <- object$frame
  studies <- unique(frame$studlab)
  rows <- lapply(studies, function(s) {
    sub    <- frame[frame$studlab == s, , drop = FALSE]
    fitted <- dr_study_slope(sub, object)
    data.frame(
      study      = s,
      levels     = length(unique(sub$dose)),
      dose_range = paste0(signif(min(sub$dose), 4), " to ",
                          signif(max(sub$dose), 4)),
      design     = if (all(is.na(sub$type))) "-" else sub$type[1],
      slope      = fitted$estimate,
      se         = fitted$se,
      stringsAsFactors = FALSE
    )
  })
  out <- list(result = object, per_study = do.call(rbind, rows))
  class(out) <- "summary.easydr_result"
  out
}

#' Study-specific slope, read from the engine's first-stage estimates
#'
#' A two-stage fit already estimated one slope per study, so those are read off
#' the fitted object rather than refitted here.
#'
#' @keywords internal
#' @noRd
dr_study_slope <- function(sub, object) {
  id <- sub$id[1]
  bi <- object$fit$bi
  Si <- object$fit$Si
  if (is.null(bi) || is.null(nrow(bi)) || id > nrow(bi)) {
    return(list(estimate = NA_real_, se = NA_real_))
  }
  est <- as.numeric(bi[id, 1])
  se  <- if (!is.null(Si) && length(Si) >= id) {
    sqrt(as.numeric(Si[[id]])[1])
  } else {
    NA_real_
  }
  list(estimate = est, se = se)
}

#' @param x a `summary.easydr_result` object.
#' @param digits significant digits for reported effects.
#' @rdname summary.easydr_result
#' @export
print.summary.easydr_result <- function(x, digits = 3, ...) {
  print(x$result, digits = digits)
  if (isTRUE(x$result$gof$available)) {
    cat("
Goodness of fit: deviance D = ", round(x$result$gof$deviance, 2),
        " (df = ", x$result$gof$df, "), p = ",
        format.pval(x$result$gof$p, digits = 3), "
", sep = "")
    cat("  R-squared ", round(x$result$gof$r2, 3), ", adjusted ",
        round(x$result$gof$r2adj, 3), "
", sep = "")
    cat("  ", x$result$gof$comparison, "
", sep = "")
  } else if (!is.null(x$result$gof$reason)) {
    cat("
Goodness of fit: not available (", x$result$gof$reason, ")
",
        sep = "")
  }
  cat("\nPer-study contribution\n")
  tab <- x$per_study
  tab$slope <- signif(tab$slope, digits + 1)
  tab$se    <- signif(tab$se, digits + 1)
  names(tab) <- c("study", "levels", "dose range", "design",
                  "slope (log scale)", "SE")
  print(tab, row.names = FALSE, right = FALSE)
  invisible(x)
}
