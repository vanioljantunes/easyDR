# The only file that talks to dosresmeta. Engine argument names live here and
# nowhere else, so the exported grammar stays the one reviewers already know.

#' Choose the covariance approximation
#'
#' Counts get Greenland-Longnecker, a reported effect gets Hamling
#' reconstruction, continuous outcomes get the mean-difference or standardised
#' mean-difference approximation. The reviewer can override the binary choice.
#'
#' @keywords internal
#' @noRd
dr_covariance <- function(intake, sm, covariance) {
  if (!is.null(covariance)) {
    allowed <- c("gl", "h", "md", "smd", "indep", "user")
    if (!covariance %in% allowed) {
      stop("easydr: covariance must be one of ",
           paste(dQuote(allowed), collapse = ", "), ".", call. = FALSE)
    }
    return(covariance)
  }
  switch(intake,
         counts        = "gl",
         effect_counts = "gl",
         effect        = "indep",
         continuous    = if (identical(sm, "SMD")) "smd" else "md",
         "gl")
}

#' Human name for a covariance approximation
#' @keywords internal
#' @noRd
dr_covariance_label <- function(x) {
  switch(x,
         gl          = "Greenland and Longnecker",
         h           = "Hamling count reconstruction",
         md          = "mean differences",
         smd         = "standardised mean differences",
         indep       = "independent contrasts, because no denominators were supplied",
         user        = "user supplied",
         x)
}

#' Build the coefficient label
#'
#' The engine prints the dose slope under the name "(Intercept)", because the
#' linear model is fitted without an intercept and the first column inherits
#' that name. A reviewer copying that would report a slope as an intercept, so
#' easydr always builds its own label and never reads the engine's coefficient
#' names (research decision D6).
#'
#' @keywords internal
#' @noRd
dr_label <- function(sm, dose_unit) {
  unit <- if (is.null(dose_unit)) "dose unit" else dose_unit
  paste0("Per 1 ", unit, " increase")
}

#' Scale name for a summary measure
#' @keywords internal
#' @noRd
dr_sm_label <- function(sm) {
  switch(sm,
         RR  = "risk ratio",
         OR  = "odds ratio",
         IRR = "incidence rate ratio",
         MD  = "mean difference",
         SMD = "standardised mean difference",
         sm)
}

#' Is this summary measure reported on a ratio scale
#' @keywords internal
#' @noRd
dr_is_ratio <- function(sm) sm %in% c("RR", "OR", "IRR")

#' Fit the dose-response model
#'
#' Translates the easydr grammar into one `dosresmeta` call. Every estimate in
#' the package comes from this call.
#'
#' @keywords internal
#' @noRd
dr_fit <- function(frame, sm, intake, proc, method, covariance, method.smd,
                   shape = NULL, mod = NULL, mod_data = NULL) {
  if (!proc %in% c("2stage", "1stage")) {
    stop("easydr: proc must be ", dQuote("2stage"), " or ", dQuote("1stage"),
         ".", call. = FALSE)
  }
  cov_code <- dr_covariance(intake, sm, covariance)

  if (is.null(shape) || identical(shape$name, "linear")) {
    # The feature 001 path, untouched: one term, named dose.
    dat <- data.frame(y = frame$y, dose = frame$dose, id = frame$id,
                      v = frame$v, stringsAsFactors = FALSE)
    form <- y ~ dose
  } else {
    design <- as.data.frame(shape$basis(frame$dose))
    dat <- data.frame(y = frame$y, id = frame$id, v = frame$v, design,
                      stringsAsFactors = FALSE)
    form <- stats::reformulate(colnames(design), response = "y")
  }

  # A moderator is a study-level covariate, so its columns ride along on the
  # analysis frame and the engine regresses the dose-response coefficients on
  # them. The formula stays the dose-response formula; `mod` is separate.
  if (!is.null(mod)) {
    dat <- cbind(dat, mod_data)
  }

  args <- list(formula = form, id = dat$id, data = dat,
               covariance = cov_code, method = method, proc = proc)
  if (!is.null(mod)) args$mod <- mod

  if (identical(intake, "continuous")) {
    # The engine derives the covariance of the mean differences from the
    # standard deviations and group sizes, so no variance is supplied.
    args$sd <- frame$cases
    args$n  <- frame$n
    if (identical(cov_code, "smd")) args$method.smd <- method.smd
  } else if (cov_code %in% c("gl", "h")) {
    args$v     <- dat$v
    if (all(is.na(frame$cases)) || all(is.na(frame$n))) {
      stop("easydr: covariance = ", dQuote(cov_code), " needs the event counts ",
           "and their denominators. Supply event.e and n.e, or let easydr ",
           "choose the approximation.", call. = FALSE)
    }
    args$type  <- frame$type
    args$cases <- frame$cases
    args$n     <- frame$n
  } else {
    args$v <- dat$v
  }

  fit <- try(do.call(dosresmeta::dosresmeta, args), silent = TRUE)
  if (inherits(fit, "try-error")) {
    stop("easydr: the dose-response model did not fit. The engine reported: ",
         conditionMessage(attr(fit, "condition")), call. = FALSE)
  }
  attr(fit, "easydr_covariance") <- cov_code
  fit
}

#' Wald test on the moderator columns only
#'
#' The right tool here, and not an information criterion. Adding a moderator
#' changes the fixed effects, so a restricted likelihood is not comparable
#' across the two models, and a two-stage criterion additionally moves with the
#' scale of the dose terms, which feature 002 documented at length. A Wald test
#' stays inside one fitted design and is free of both objections.
#'
#' @keywords internal
#' @noRd
dr_moderator_test <- function(fit, terms) {
  if (length(terms) == 0L) {
    return(list(estimable = FALSE,
                reason = "no moderator columns to test"))
  }
  w <- try(dosresmeta::waldtest(b = coef(fit), Sigma = vcov(fit),
                                Terms = terms), silent = TRUE)
  if (inherits(w, "try-error")) {
    return(list(estimable = FALSE,
                reason = "the engine could not test the moderator columns"))
  }
  list(estimable = TRUE,
       chi2  = as.numeric(w$chitest["chi2"]),
       df    = as.integer(w$chitest["df"]),
       p     = as.numeric(w$chitest["P"]),
       terms = terms)
}

#' Departure from linearity, from the engine's own Wald test
#'
#' With more than one term the nonlinear terms can be tested jointly against
#' zero, which is the test a reviewer means by departure from linearity. It
#' depends on the first basis column being the dose itself, which is why the
#' spline basis is Harrell's and not a natural spline (research decision D1).
#'
#' A one-term fit admits no such test, so it reports the reason instead, exactly
#' as feature 001 does.
#'
#' @keywords internal
#' @noRd
dr_nonlinearity_test <- function(fit, shape, prov = NULL) {
  if (is.null(shape) || shape$terms < 2L) {
    # The feature 001 report, unchanged: the reason, plus how many studies
    # could have informed a test had one been admissible.
    return(list(
      estimable = FALSE,
      reason = if (is.null(shape) || identical(shape$name, "linear")) {
        "a linear fit has one term, so there is nothing to test it against"
      } else {
        paste0("the ", shape$label,
               " fit has one term, so there is nothing to test it against")
      },
      studies_with_three_levels = if (is.null(prov$n_with_three)) 0L else prov$n_with_three,
      studies = if (is.null(prov$dose_levels)) 0L else length(prov$dose_levels)))
  }
  terms_tested <- 2:shape$terms
  w <- try(dosresmeta::waldtest(b = coef(fit), Sigma = vcov(fit),
                                Terms = terms_tested), silent = TRUE)
  if (inherits(w, "try-error")) {
    return(list(estimable = FALSE,
                reason = paste0("the engine could not compute the Wald test: ",
                                conditionMessage(attr(w, "condition")))))
  }
  list(estimable = TRUE,
       statistic = as.numeric(w$chitest["chi2"]),
       df        = as.numeric(w$chitest["df"]),
       p         = as.numeric(w$chitest["P"]),
       terms     = terms_tested,
       comparison = paste0("Wald test on the ", length(terms_tested),
                           " nonlinear term",
                           if (length(terms_tested) == 1L) "" else "s",
                           " against the linear fit"))
}

#' Goodness of fit, from the engine
#' @keywords internal
#' @noRd
dr_gof <- function(fit) {
  g <- try(dosresmeta::gof(fit), silent = TRUE)
  if (inherits(g, "try-error") || is.null(g$deviance)) {
    return(list(available = FALSE,
                reason = "the engine returned no goodness-of-fit statistics for this fit"))
  }
  # gof() returns deviance as a list of D, df and p, not a named vector.
  d <- g$deviance
  list(available = TRUE,
       deviance = as.numeric(d$D),
       df       = as.numeric(d$df),
       p        = as.numeric(d$p),
       r2       = as.numeric(g$R2),
       r2adj    = as.numeric(g$R2adj),
       comparison = "deviance of the fitted model against a saturated fit")
}

#' Pull the heterogeneity summary the engine already computes
#'
#' The engine's own `summary()` reports the Cochran Q test and I-squared.
#' easydr reads them rather than computing a competing statistic
#' (research decision D2, constitution Principle V).
#'
#' @keywords internal
#' @noRd
dr_heterogeneity <- function(fit, single_study) {
  if (single_study) {
    return(list(
      estimable = FALSE,
      reason = "one study only, so between-study heterogeneity is not estimable"))
  }
  # qtest is the mixmeta generic; dosresmeta supplies the method for its own fits.
  q <- try(qtest(fit), silent = TRUE)
  if (inherits(q, "try-error") || is.null(q$Q)) {
    return(list(
      estimable = FALSE,
      reason = "the engine returned no heterogeneity test for this fit"))
  }
  Q  <- as.numeric(q$Q)[1]
  df <- as.numeric(q$df)[1]
  p  <- as.numeric(q$pvalue)[1]
  list(estimable = TRUE, Q = Q, df = df, p = p,
       i2 = 100 * max(0, (Q - df) / Q))
}
