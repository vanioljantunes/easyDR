# Target doses.
#
# A curve says what happens at every dose. A guideline has to name one. The
# quantity that bridges them is the target dose: the smallest dose reaching a
# stated fraction of the largest effect the curve achieves.
#
# Three things have to be said about it carefully.
#
# It is defined inside the range the studies reported, because outside that
# range the curve is the model talking to itself.
#
# The largest effect may sit at the end of that range rather than at a turning
# point. A straight line has no maximum; it has an edge. Reporting the edge as
# ED max would dress a range limit up as a finding, so it is flagged.
#
# It has no standard error. ED gamma is a non-linear functional of the
# coefficients, and for a spline read off a grid there is no derivative for
# the delta method to use. The interval comes from the bootstrap.

#' Read the target doses off one curve
#'
#' The direction is taken from the largest absolute effect, so a treatment
#' that lowers a score and an exposure that raises a risk are handled the same
#' way without the caller declaring which they have.
#'
#' @keywords internal
#' @noRd
dr_ed_from_curve <- function(dose, est, gamma) {
  ok   <- is.finite(dose) & is.finite(est)
  dose <- dose[ok]
  est  <- est[ok]
  if (length(dose) < 3L) return(NULL)
  o    <- order(dose)
  dose <- dose[o]
  est  <- est[o]

  i     <- which.max(abs(est))
  emax  <- est[i]
  edmax <- dose[i]
  if (!is.finite(emax) || isTRUE(all.equal(emax, 0))) return(NULL)

  ed <- vapply(gamma, function(g) {
    target  <- g * emax
    # The first dose that reaches the target, in whichever direction the
    # curve runs, and never past the dose of the maximum itself.
    reached <- if (emax > 0) est >= target else est <= target
    reached[dose > edmax] <- FALSE
    j <- which(reached)
    if (!length(j)) return(NA_real_)
    j <- j[1L]
    if (j == 1L) return(dose[1L])
    # Interpolated between the two grid points either side, so the answer
    # does not move with the density of the grid.
    x0 <- dose[j - 1L]
    x1 <- dose[j]
    y0 <- est[j - 1L]
    y1 <- est[j]
    if (isTRUE(all.equal(y1, y0))) return(x1)
    x0 + (target - y0) * (x1 - x0) / (y1 - y0)
  }, numeric(1))

  list(ed = ed, emax = emax, ed_max = edmax,
       at_edge = isTRUE(all.equal(edmax, max(dose))))
}

#' The dose grid a target dose is read on
#' @keywords internal
#' @noRd
dr_ed_grid <- function(x, n = 501L) {
  seq(x$counts$dose_min, x$counts$dose_max, length.out = n)
}

#' Check a vector of fractions
#' @keywords internal
#' @noRd
dr_check_gamma <- function(gamma) {
  if (!is.numeric(gamma) || !length(gamma) || anyNA(gamma) ||
      any(gamma <= 0) || any(gamma > 1)) {
    stop("easydr: `gamma` is the fraction of the maximum effect to reach, ",
         "above 0 and at most 1, such as gamma = c(0.5, 0.8).",
         call. = FALSE)
  }
  invisible(NULL)
}

#' The fitted effect at any doses, on the model scale
#' @keywords internal
#' @noRd
dr_curve_at <- function(x, doses) {
  b    <- as.numeric(stats::coef(x$fit))
  xref <- x$reference$curve
  vapply(doses, function(d) {
    as.numeric(dr_contrast(x$shape, d, xref) %*% b)
  }, numeric(1))
}

#' The dose reaching a fraction of the maximum effect
#'
#' `ed_gamma()` reads one number off a fitted curve: the smallest dose at
#' which the effect reaches `gamma` times the largest effect the curve
#' achieves inside the range the studies reported. `gamma = 0.5` gives the
#' ED50, `gamma = 0.8` the ED80, and the dose of the largest effect itself is
#' reported alongside as ED max.
#'
#' The direction is read from the curve, so an exposure that raises a risk and
#' a treatment that lowers a score are both handled without being told which
#' is which.
#'
#' This gives the point estimate only. For an interval use
#' [easydr_target_dose()], which bootstraps it.
#'
#' @param x a fitted result from [easydr()].
#' @param gamma one or more fractions of the maximum effect, above 0 and at
#'   most 1.
#'
#' @return A list with `ed`, one dose per `gamma`; `ed_max`, the dose of the
#'   largest effect; `emax`, that effect; and `at_edge`, `TRUE` when the
#'   largest effect sits at the top of the observed range rather than at a
#'   turning point.
#'
#' @examples
#' data(aripiprazole)
#' easydr.object <- easydr(aripiprazole, sm = "MD",
#'                         shape = "hyperbolic_emax", dose.unit = "mg/day")
#' ed_gamma(easydr.object, gamma = c(0.5, 0.8))
#'
#' @export
ed_gamma <- function(x, gamma = c(0.5, 0.8)) {
  dr_require_single(x, "ed_gamma()")
  dr_check_gamma(gamma)
  grid <- dr_ed_grid(x)
  out  <- dr_ed_from_curve(grid, dr_curve_at(x, grid), gamma)
  if (is.null(out)) {
    stop("easydr: this curve has no effect to take a fraction of. A flat ",
         "fit has no target dose.", call. = FALSE)
  }
  names(out$ed) <- paste0("ED", format(100 * gamma, trim = TRUE))
  out$gamma <- gamma
  out$unit  <- x$dose_unit
  out
}

#' Target doses, with bootstrap confidence intervals
#'
#' Reports ED max, the dose of the largest effect the curve reaches inside the
#' observed range, and one target dose per `gamma`: the smallest dose reaching
#' that fraction of it. Each comes with a percentile interval from a
#' parametric bootstrap.
#'
#' A target dose has no standard error, being a non-linear functional of the
#' coefficients with no derivative for a spline read off a grid, so the
#' bootstrap is not a refinement here but the only route to an interval. Every
#' draw refits the whole model, and for a shape that profiles a parameter it
#' profiles again, so the uncertainty in where the curve bends is carried
#' rather than conditioned away.
#'
#' @param x a fitted result from [easydr()].
#' @param gamma fractions of the maximum effect, above 0 and at most 1.
#' @param B bootstrap draws. Several hundred is usual. A profiled shape refits
#'   its whole grid on every draw, so this is the slow argument.
#' @param level the confidence level for the percentile interval. `NULL` uses
#'   the level the model was fitted at.
#' @param reprofile re-run the profile on each draw, for a shape that has one.
#' @param seed a seed, so a reported interval can be reproduced.
#' @param quiet `FALSE` reports progress.
#'
#' @return An object of class `easydr_target`, holding `table`, one row per
#'   reported dose.
#'
#' @examples
#' \donttest{
#' data(aripiprazole)
#' easydr.object <- easydr(aripiprazole, sm = "MD",
#'                         shape = "hyperbolic_emax", dose.unit = "mg/day")
#' easydr_target_dose(easydr.object, B = 200, seed = 1)
#' }
#'
#' @export
easydr_target_dose <- function(x, gamma = c(0.5, 0.8), B = 500L,
                               level = NULL, reprofile = TRUE, seed = NULL,
                               quiet = TRUE) {
  dr_require_single(x, "easydr_target_dose()")
  dr_check_gamma(gamma)
  dr_check_boot(B)
  if (is.null(level)) level <- x$provenance$level

  point <- ed_gamma(x, gamma)
  grid  <- dr_ed_grid(x)
  xref  <- x$reference$curve

  if (!is.null(seed)) set.seed(seed)
  draws <- matrix(NA_real_, nrow = B, ncol = length(gamma) + 1L)
  for (b in seq_len(B)) {
    d <- dr_boot_once(x, reprofile)
    if (is.null(d)) next
    est <- try(dr_boot_predict(d, grid, xref), silent = TRUE)
    if (inherits(est, "try-error")) next
    got <- dr_ed_from_curve(grid, est, gamma)
    if (is.null(got)) next
    draws[b, ] <- c(got$ed, got$ed_max)
    if (!isTRUE(quiet) && b %% 50L == 0L) {
      message("easydr: bootstrap draw ", b, " of ", B)
    }
  }
  draws <- draws[stats::complete.cases(draws), , drop = FALSE]
  if (nrow(draws) < 20L) {
    stop("easydr: only ", nrow(draws), " of ", B, " draws gave a target ",
         "dose, which is too few to read an interval from. A curve that is ",
         "nearly flat has no stable target dose.", call. = FALSE)
  }
  lim <- dr_boot_limits(draws, level)

  tab <- data.frame(
    Quantity = c(names(point$ed), "ED max"),
    Dose     = c(point$ed, point$ed_max),
    Lower    = lim$lb,
    Upper    = lim$ub,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  rownames(tab) <- NULL

  structure(list(table = tab, gamma = gamma, emax = point$emax,
                 at_edge = point$at_edge, unit = x$dose_unit,
                 sm_label = x$sm_label, level = level,
                 draws = nrow(draws), B = B,
                 # Whether the profile actually ran again, not merely
                 # whether it was asked for: a shape whose parameter
                 # the caller supplied has nothing to profile.
                 reprofiled = isTRUE(reprofile) &&
                   !is.null(x$shape$profiled_over),
                 degenerate = all(abs(lim$ub - lim$lb) < 1e-8),
                 # A percentile interval need not contain the point
                 # estimate. When it does not, the bootstrap
                 # distribution is shifted away from it, which is
                 # worth saying rather than leaving to be noticed.
                 excludes_point = any(tab$Dose < lim$lb - 1e-8 |
                                      tab$Dose > lim$ub + 1e-8),
                 range = c(x$counts$dose_min, x$counts$dose_max)),
            class = "easydr_target")
}

#' Print target doses
#'
#' @param x an object from [easydr_target_dose()].
#' @param digits digits for the reported doses.
#' @param ... ignored.
#'
#' @return `x`, invisibly. Called for the printed output.
#'
#' @export
print.easydr_target <- function(x, digits = 2, ...) {
  unit <- if (is.null(x$unit)) "dose unit" else x$unit
  cat("Target doses (", unit, ")\n", sep = "")
  cat("Largest fitted effect: ", format(round(x$emax, 3)), " ",
      x$sm_label, " at ", format(round(x$table$Dose[nrow(x$table)], digits)),
      " ", unit, "\n", sep = "")
  cat("Read inside the observed range, ", x$range[1L], " to ", x$range[2L],
      " ", unit, "\n\n", sep = "")

  show <- x$table
  for (nm in c("Dose", "Lower", "Upper")) {
    show[[nm]] <- format(round(show[[nm]], digits), nsmall = digits)
  }
  names(show)[3:4] <- paste0(format(100 * x$level), "% ", c("lower", "upper"))
  print(show, row.names = FALSE, right = FALSE)

  cat("\nPercentile interval from ", x$draws, " of ", x$B,
      " bootstrap draws", if (x$reprofiled) ", re-profiled on each draw",
      ".\n", sep = "")
  if (isTRUE(x$degenerate)) {
    cat("\nEvery draw gave the same target dose, so the interval has",
        "\nno width. This shape has one coefficient, and a target",
        "\ndose depends on the shape of the curve rather than on its",
        "\nheight, so nothing that varies between draws can move it.",
        "\nProfile the parameter rather than supplying it if you want",
        "\nan interval that carries where the curve bends.\n")
  }
  if (isTRUE(x$excludes_point)) {
    cat("\nAt least one interval does not contain its own estimate. A",
        "\npercentile interval is read off where the draws fell, and",
        "\nwhen they fall mostly to one side of the estimate it can",
        "\nexclude it. That is a sign the estimate sits near an edge",
        "\nof what this data can support, not an arithmetic fault.\n")
  }
  if (isTRUE(x$at_edge)) {
    cat("\nThe largest effect sits at the top of the observed range rather",
        "\nthan at a turning point, so ED max is a range limit and not an",
        "\nestimate of where the curve stops rising. Report it as such.\n")
  }
  invisible(x)
}

#' The goodness of fit, as one row a manuscript can quote
#'
#' The deviance against a saturated model asks whether the curve could have
#' produced the doses that were reported. The R-squared says how much of the
#' variation between them it accounts for. The Wald test for departure from
#' linearity answers a different question, which is why both belong in the
#' same report rather than one standing for the other.
#'
#' @keywords internal
#' @noRd
dr_gof_line <- function(x) {
  g <- x$gof
  if (is.null(g) || is.null(g$deviance)) {
    return(data.frame(Statistic = "Goodness of fit",
                      Value = "not estimable",
                      check.names = FALSE, stringsAsFactors = FALSE))
  }
  nl <- x$nonlinearity
  data.frame(
    Statistic = c("Deviance (df), p", "R-squared",
                  "Departure from linearity"),
    Value = c(sprintf("%.2f (%d), %s", g$deviance, g$df, dr_p(g$p)),
              sprintf("%.3f", g$r2),
              if (isTRUE(nl$estimable)) dr_p_lab(nl$p) else "not estimable"),
    check.names = FALSE, stringsAsFactors = FALSE)
}
