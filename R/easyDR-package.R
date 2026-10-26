#' easyDR: dose-response meta-analysis, argued like meta
#'
#' @description
#' One call turns a table of exposure categories into a pooled dose-response.
#' easyDR is a front end for \code{\link[dosresmeta]{dosresmeta}}: that package
#' does the estimation, easyDR prepares the data, calls it once, and reports the
#' answer in the words a systematic reviewer uses.
#'
#' The arguments are the ones a `meta` user already knows. There is no
#' `formula`. You write `sm = "RR"`, not `logrr ~ dose`.
#'
#' @section The family:
#' \describe{
#'   \item{[easydr()]}{fit it. One entry point for every outcome; `sm` decides
#'     which columns are read.}
#'   \item{[easydr_forest()]}{the study trends and the pooled one, laid out as
#'     a forest plot.}
#'   \item{[easydr_curve()]}{the fitted dose-response with its band.}
#'   \item{[easydr_loo()]}{leave one study out at a time.}
#'   \item{[easydr_funnel()]}{small-study effects, with Egger's test.}
#' }
#'
#' @section What your table must contain:
#' One row per exposure category per study, in long form. Nothing has to be
#' reshaped, sorted or aggregated.
#'
#' \describe{
#'   \item{\code{studlab}}{study label. Rows sharing a label are one study.}
#'   \item{\code{dose}}{the numeric dose you assigned to that category. A dose
#'     is never inferred from a category range; that judgement is yours.}
#'   \item{\code{ref}}{optional logical, marking each study's reference
#'     category. When absent one is inferred and the guesses are printed.}
#'   \item{outcome columns}{\code{event.e} and \code{n.e} for counts,
#'     \code{TE} with \code{lower} and \code{upper} for a reported effect, or
#'     \code{mean.e}, \code{sd.e} and \code{n.e} for continuous outcomes.}
#' }
#'
#' @section Two tables to copy the shape from:
#' [alcohol_cancer] is a binary-outcome extraction, [aripiprazole] a continuous
#' one. Both are real published reviews, and both run in one call:
#'
#' \preformatted{
#' easydr(alcohol_cancer, sm = "RR", shape = "linear", dose.unit = "g/day")
#' easydr(aripiprazole, sm = "SMD", shape = "linear", dose.unit = "mg/day")
#' }
#'
#' @section Fitting a curve instead of a line:
#' Change one argument. `shape = "rcs"` fits a restricted cubic spline,
#' `shape = "fp"` a fractional polynomial, and `"quadratic"` and `"log"` do what
#' they say.
#'
#' A curve prints a test for departure from linearity, which a straight line
#' cannot: there is nothing for one term to be tested against. A curve also
#' prints no single slope, because no single number describes it. Read the
#' predicted effects instead.
#'
#' @section What easyDR does not do:
#' Individual participant data, network dose-response and Bayesian fitting are
#' out of scope. Dose units are assumed already harmonised across studies:
#' easyDR warns when the per-study dose ranges look inconsistent, but it never
#' converts units.
#'
#' @references
#' Orsini N, Bellocco R, Greenland S (2006). Generalized least squares for trend
#' estimation of summarized dose-response data. Stata Journal 6, 40-57.
#'
#' Crippa A, Orsini N (2016). Multivariate dose-response meta-analysis: the
#' dosresmeta R package. Journal of Statistical Software 72, 1-15.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats qt pchisq setNames coef vcov quantile
#' @importFrom stats AIC BIC reformulate lm
#' @importFrom Hmisc rcspline.eval
#' @importFrom utils head flush.console
#' @importFrom mixmeta qtest
#' @importFrom rlang .data
NULL
