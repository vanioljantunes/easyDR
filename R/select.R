#' Rank candidate dose-response shapes
#'
#' @description
#' Fits a bounded set of candidate shapes to the same data and returns them
#' ranked, with the statistics the ranking used. It hands back the best-ranked
#' fit so you need not refit, and it changes nothing else: a shape you stated in
#' your own call is never rewritten, replaced or refitted by this function.
#'
#' @section What ranking shapes does and does not license:
#' An information criterion says which candidate describes this data most
#' economically. It does not say the winning shape is true, that the runner-up is
#' wrong, or that a curve is real. With dozens of candidates the best score is
#' partly luck, and the gap between the top few is usually smaller than the gap
#' between any of them and the truth. That is why the whole table comes back
#' rather than one answer, why ties are shown, and why the printed output calls
#' the result a proposal.
#'
#' Use it to narrow the field and to see whether a curve is warranted at all,
#' which is what including the straight line among the candidates is for. Then
#' choose a shape for reasons you can write down: the biology, the dose range
#' your readers care about, what comparable reviews did. Report that you
#' searched, and over what.
#'
#' @section Why the candidates are refitted one-stage by maximum likelihood:
#' Two things make a naive ranking wrong, and both are corrected here.
#'
#' Restricted maximum likelihood estimates heterogeneity after conditioning on
#' the fixed effects, so its likelihood depends on what those fixed effects are.
#' Comparing a line against a curve by REML information criteria compares two
#' numbers that are not on the same scale.
#'
#' Two-stage pooling is worse for this purpose. Its second stage treats each
#' study's fitted coefficients as the outcome, and those coefficients shrink or
#' grow with the scale of the design. Multiplying the dose by 1000, which changes
#' no model at all, moves the two-stage AIC on the bundled alcohol data from
#' -64.5 to -175.1. A criterion that rewards a change of units cannot rank
#' shapes. One-stage pooling keeps the same outcome for every candidate, and the
#' same comparison gives -6.1 either way.
#'
#' So every candidate is refitted one-stage by maximum likelihood unless you pass
#' `proc` or `method.tau` yourself, and the printed output says so. The fit handed
#' back in `best` is that one-stage fit, matching the ranking. Refit it two-stage
#' yourself if that is how you intend to report it.
#'
#' @section The candidate set:
#' Bounded and documented, so the work is reproducible and the count is
#' reportable: the straight line, a quadratic, a logarithmic dose, the 36
#' fractional-polynomial power pairs of the engine's own grid, and restricted
#' cubic splines with 3, 4 and 5 knots. Candidates the data cannot support fail
#' and are listed with their reason rather than quietly dropped.
#'
#' @param data A data frame in long form, one row per dose category per study,
#'   exactly as [easydr()] takes it.
#' @param ... Further arguments passed unchanged to [easydr()], such as `sm`,
#'   `dose.unit` or `studlab`. `sm` decides whether the outcome is binary or
#'   continuous, as it does everywhere else. Do not pass `shape`, `knots` or
#'   `powers`: those are what this function varies.
#' @param candidates A list of candidate specifications, each a list with a
#'   `shape` and optionally `knots` or `powers`. `NULL` uses the documented set
#'   above.
#' @param rank_on `"AIC"` or `"BIC"`. Named in the printed output, never
#'   implicit.
#' @param tie_tolerance Candidates within this much of the best on `rank_on` are
#'   reported as tied. Two is the conventional threshold for an information
#'   criterion.
#' @param quiet `FALSE` prints progress as the loop runs, because this is the one
#'   slow operation in the package.
#'
#' @return An object of class `easydr_selection` with the ranked table
#'   (`ranked`), the statistic ranked on (`ranked_on`), the tie group (`ties`),
#'   the best-ranked fit (`best`), the count attempted (`tried`) and the elapsed
#'   seconds (`elapsed`).
#'
#' @seealso [easydr()], which fits the one shape you have decided on.
#'
#' @references
#' Royston P, Altman DG (1994). Regression using fractional polynomials of
#' continuous covariates. Applied Statistics 43, 429-467.
#'
#' @examples
#' \donttest{
#' # Around forty fits, so a few seconds.
#' sel <- easydr_select(alcohol_cancer, sm = "IRR", dose.unit = "g/day")
#' sel
#'
#' # The whole field, successes and failures alike.
#' head(sel$ranked, 10)
#'
#' # The best-ranked fit is ready to use, but adopting it is your decision.
#' sel$best
#' }
#'
#' @export
easydr_select <- function(data, ..., candidates = NULL, rank_on = "AIC",
                          tie_tolerance = 2, quiet = FALSE) {
  call    <- match.call()
  rank_on <- match.arg(rank_on, c("AIC", "BIC"))
  dots    <- list(...)

  # `sm` picks the outcome kind here exactly as it does in easydr(), so the
  # caller never names a fitting function.
  sm      <- if (is.null(dots$sm)) "RR" else match.arg(
               dots$sm, c("RR", "OR", "IRR", "MD", "SMD"))
  fn      <- if (sm %in% c("MD", "SMD")) easydrcont else easydrbin
  fn_name <- "easydr"

  if (any(c("shape", "knots", "powers") %in% names(dots))) {
    stop("easydr: easydr_select varies the shape itself, so `shape`, `knots` and ",
         "`powers` cannot be passed to it. Drop them, or fit the shape you ",
         "already have in mind with easydr().", call. = FALSE)
  }

  # Comparing information criteria across different fixed effects needs a
  # likelihood that does not condition on them, so candidates are refitted by
  # maximum likelihood unless the caller insists otherwise.
  # A common-effect request beats the forced estimator, because the pair
  # chooses the model. Saying "refitted by maximum likelihood" afterwards would
  # describe a fit that was not made, so the forcing is recorded as not having
  # happened and the caller is told why.
  common_asked <- isTRUE(dots$common) ||
    ("random" %in% names(dots) && !isTRUE(dots$random))
  forced_ml  <- !"method.tau" %in% names(dots) && !common_asked
  forced_1st <- !"proc" %in% names(dots)
  if (common_asked) {
    warning("easydr: you asked for a common-effect fit, so the candidates are ",
            "not refitted by maximum likelihood and the ranking compares ",
            "fits that assume no heterogeneity at all. Drop `common` and ",
            "`random` to rank random-effects fits.", call. = FALSE)
  }
  if (forced_ml) dots$method.tau <- "ML"
  if (forced_1st) dots$proc <- "1stage"
  if (!forced_1st && identical(dots$proc, "2stage")) {
    warning("easydr: ranking two-stage fits by an information criterion is not ",
            "reliable. The second stage treats each study coefficient vector as ",
            "the outcome, so the criterion moves with the scale of the dose ",
            "terms and rewards a change of units. Drop `proc` to rank ",
            "one-stage fits.", call. = FALSE)
  }

  cands <- if (is.null(candidates)) dr_candidates() else candidates
  n     <- length(cands)
  start <- proc.time()[["elapsed"]]
  rows  <- vector("list", n)
  fits  <- vector("list", n)

  for (i in seq_len(n)) {
    spec <- cands[[i]]
    if (!quiet) {
      cat("\rFitting candidate ", i, " of ", n, ": ",
          format(dr_candidate_label(spec), width = 30), sep = "")
      utils::flush.console()
    }
    args <- c(list(data = data), dots, spec)
    fit  <- try(suppressWarnings(do.call(fn, args)), silent = TRUE)

    if (inherits(fit, "try-error")) {
      # Rule V18: a candidate that fails is recorded and the loop continues.
      rows[[i]] <- dr_candidate_row(
        spec, i, status = "failed",
        reason = conditionMessage(attr(fit, "condition")))
      next
    }
    rows[[i]] <- dr_candidate_row(spec, i, status = "fitted", fit = fit)
    fits[[i]] <- fit
  }
  if (!quiet) cat("\r", strrep(" ", 60), "\r", sep = "")

  ranked      <- do.call(rbind, rows)
  ok          <- ranked$status == "fitted"
  ranked$rank <- NA_integer_
  ranked$rank[ok] <- rank(ranked[[rank_on]][ok], ties.method = "min")
  ranked <- ranked[order(ranked$rank, na.last = TRUE), , drop = FALSE]
  rownames(ranked) <- NULL

  best <- if (any(ok)) fits[[ranked$candidate[1L]]] else NULL
  ties <- dr_ties(ranked, rank_on, tie_tolerance)

  structure(list(ranked = ranked, ranked_on = rank_on, ties = ties,
                 best = best, tried = n, fitted = sum(ok), failed = sum(!ok),
                 elapsed = round(proc.time()[["elapsed"]] - start, 1),
                 forced_ml = forced_ml, forced_1stage = forced_1st,
                 tie_tolerance = tie_tolerance,
                 fn_name = fn_name, call = call),
            class = "easydr_selection")
}

#' The documented candidate set
#' @keywords internal
#' @noRd
dr_candidates <- function() {
  out <- list(list(shape = "linear"),
              list(shape = "quadratic"),
              list(shape = "log"))
  grid <- dosresmeta::fpgrid()
  for (i in seq_len(nrow(grid))) {
    out[[length(out) + 1L]] <- list(shape = "fp",
                                    powers = as.numeric(grid[i, ]))
  }
  for (k in 3:5) {
    out[[length(out) + 1L]] <- list(shape = "rcs", knots = k)
  }
  out
}

#' @keywords internal
#' @noRd
dr_candidate_label <- function(spec) {
  # Every branch must collapse to one string: a knot vector left uncollapsed
  # would produce one table row per knot.
  switch(spec$shape,
         linear    = "linear",
         quadratic = "quadratic",
         log       = "log dose",
         fp        = paste0("fp powers ", paste(spec$powers, collapse = ", ")),
         rcs       = if (length(spec$knots) == 1L) {
           paste0("rcs ", spec$knots, " knots")
         } else {
           paste0("rcs knots ", paste(signif(spec$knots, 4), collapse = ", "))
         },
         spec$shape)
}

#' One row of the ranked table
#' @keywords internal
#' @noRd
dr_candidate_row <- function(spec, index, status, reason = NA_character_,
                             fit = NULL) {
  data.frame(
    candidate = index,
    shape     = dr_candidate_label(spec),
    terms     = if (is.null(fit)) NA_integer_ else fit$shape$terms,
    AIC       = if (is.null(fit)) NA_real_ else as.numeric(stats::AIC(fit$fit)),
    BIC       = if (is.null(fit)) NA_real_ else as.numeric(stats::BIC(fit$fit)),
    deviance_p = if (is.null(fit) || !isTRUE(fit$gof$available)) {
      NA_real_
    } else {
      fit$gof$p
    },
    nonlinear_p = if (is.null(fit) || !isTRUE(fit$nonlinearity$estimable)) {
      NA_real_
    } else {
      fit$nonlinearity$p
    },
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

#' Candidates indistinguishable from the best on the ranking statistic
#' @keywords internal
#' @noRd
dr_ties <- function(ranked, rank_on, tolerance) {
  ok <- ranked$status == "fitted"
  if (!any(ok)) return(character())
  best <- min(ranked[[rank_on]][ok], na.rm = TRUE)
  tied <- ranked$shape[ok & (ranked[[rank_on]] - best) <= tolerance]
  if (length(tied) < 2L) character() else tied
}

#' Print a shape selection
#'
#' @param x an object from [easydr_select()].
#' @param n how many ranked candidates to show.
#' @param ... ignored.
#'
#' @return `x`, invisibly. Called for the printed output.
#'
#' @examples
#' \donttest{
#' print(easydr_select(alcohol_cancer, sm = "IRR"), n = 5)
#' }
#'
#' @export
print.easydr_selection <- function(x, n = 10, ...) {
  cat("Shape selection for a dose-response meta-analysis\n")
  cat(x$tried, " candidates tried, ", x$fitted, " fitted, ", x$failed,
      " failed, ", x$elapsed, " seconds\n", sep = "")
  cat("Ranked on ", x$ranked_on, ". Lower is better.\n", sep = "")
  # Say what was actually forced. Claiming maximum likelihood when a
  # common-effect request overrode it would describe a fit that was not made.
  if (isTRUE(x$forced_ml) && isTRUE(x$forced_1stage)) {
    cat("Refitted one-stage by maximum likelihood. Restricted maximum
")
    cat("  likelihood cannot compare models with different dose terms, and a
")
    cat("  two-stage criterion moves with the scale of those terms.
")
  } else if (isTRUE(x$forced_1stage)) {
    cat("Refitted one-stage, but not by maximum likelihood, because a
")
    cat("  common-effect fit was asked for. A two-stage criterion moves with
")
    cat("  the scale of the dose terms.
")
  } else if (isTRUE(x$forced_ml)) {
    cat("Refitted by maximum likelihood. Restricted maximum likelihood cannot
")
    cat("  compare models with different dose terms.
")
  }
  cat("\n")

  fitted_rows <- x$ranked[x$ranked$status == "fitted", , drop = FALSE]
  if (nrow(fitted_rows) == 0L) {
    cat("No candidate fitted. See x$ranked for the reason against each.\n")
    return(invisible(x))
  }
  show <- utils::head(fitted_rows, n)
  tab <- data.frame(
    rank  = show$rank,
    shape = show$shape,
    terms = show$terms,
    AIC   = round(show$AIC, 1),
    BIC   = round(show$BIC, 1),
    `deviance p` = ifelse(is.na(show$deviance_p), "-",
                          format(round(show$deviance_p, 3))),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  print(tab, row.names = FALSE, right = FALSE)
  if (nrow(fitted_rows) > n) {
    cat("... ", nrow(show), " of ", nrow(fitted_rows),
        " fitted candidates shown\n", sep = "")
  }
  cat("\n")

  if (length(x$ties) > 1L) {
    cat("Within ", x$tie_tolerance, " ", x$ranked_on, " of the best: ",
        paste(x$ties, collapse = "; "), "\n", sep = "")
  }

  cat("This ranking is a proposal, not a decision. To adopt the best-ranked\n")
  cat("shape, state it yourself:\n")
  cat("  ", dr_adopt_call(x), "\n", sep = "")
  if (x$failed > 0L) {
    cat(x$failed, " candidate", if (x$failed == 1L) "" else "s",
        " failed. See x$ranked for the reason against each.\n", sep = "")
  }
  invisible(x)
}

#' The call that would adopt the best-ranked shape
#' @keywords internal
#' @noRd
dr_adopt_call <- function(x) {
  lbl <- x$ranked$shape[1L]
  shp <- if (grepl("^fp", lbl)) {
    paste0('shape = "fp", powers = c(', sub("^fp powers ", "", lbl), ")")
  } else if (grepl("^rcs knots ", lbl)) {
    paste0('shape = "rcs", knots = c(', sub("^rcs knots ", "", lbl), ")")
  } else if (grepl("^rcs", lbl)) {
    paste0('shape = "rcs", knots = ',
           sub(" knots$", "", sub("^rcs ", "", lbl)))
  } else if (identical(lbl, "log dose")) {
    'shape = "log"'
  } else {
    paste0('shape = "', lbl, '"')
  }
  paste0(x$fn_name, "(data, ..., ", shp, ")")
}
