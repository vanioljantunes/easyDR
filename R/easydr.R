#' Dose-response meta-analysis in one call
#'
#' @description
#' Pools a dose-response relationship across studies. Give it one long table,
#' one row per exposure category per study, say which summary measure the
#' studies report, and name the shape you want fitted. It returns the pooled
#' dose-response, its heterogeneity, the prediction interval, the predicted
#' effects, a viability check on the studies, and a name for the curve.
#'
#' The arguments are the ones a `meta` user already knows: `sm`, `studlab`,
#' `event.e`, `n.e`, `TE`, `lower`, `upper`, `method.tau`, `level`,
#' `common` and `random`. You write `sm = "RR"` where another package would
#' ask for `logrr ~ dose`, and easyDR derives the contrasts, the standard
#' errors and the within-study covariances from your table.
#'
#' @section Which measure, and so which columns:
#' `sm` decides what kind of outcome this is and therefore which columns are
#' read.
#'
#' | `sm` | Outcome | Columns read |
#' |:---|:---|:---|
#' | `"RR"` | risk ratio | `event.e`, `n.e`, or `TE` with `lower`/`upper` |
#' | `"OR"` | odds ratio | the same, from a case-control design |
#' | `"IRR"` | incidence rate ratio | `event.e` with person-time in `n.e` |
#' | `"MD"` | mean difference | `mean.e`, `sd.e`, `n.e` |
#' | `"SMD"` | standardised mean difference | `mean.e`, `sd.e`, `n.e` |
#'
#' Column names are guessed from the table when you leave them `NULL`, and the
#' guesses are printed. Naming them is always allowed and never penalised.
#'
#' @section Why the categories inside one study are correlated:
#' A study does not report independent effects. Every category is compared
#' with the same reference category, so the reference group's random error is
#' shared by all of that study's estimates. Treat them as independent and the
#' slope looks more certain than it is. Published reviews get this wrong
#' often.
#'
#' Published tables almost never include those covariances, so they are
#' reconstructed. With counts, the Greenland and Longnecker method recovers a
#' set of counts consistent with the reported effects and derives the
#' covariance from them; `covariance = "h"` uses Hamling's method instead. Both
#' need denominators, so a table carrying only an effect with its limits
#' supports neither: the contrasts are then treated as independent and the
#' output says so.
#'
#' @section Why you must name the shape:
#' `shape` has no default and `easydr()` will not choose one. The shape is a
#' claim about the biology, and it belongs in your protocol with a reason
#' beside it: what the dose range is, what comparable reviews fitted, whether a
#' threshold or a flattening is plausible for this exposure. Letting software
#' pick it turns a reviewable decision into an invisible one.
#'
#' The six you can name:
#'
#' | `shape` | Fits | Terms |
#' |:---|:---|:---|
#' | `"linear"` | a straight line through the dose | 1 |
#' | `"log"` | a straight line through log dose | 1 |
#' | `"quadratic"` | dose and dose squared | 2 |
#' | `"cubic"` | dose, dose squared and dose cubed | 3 |
#' | `"rcs"` | a restricted cubic spline, `knots` to place them | 2 or more |
#' | `"fp"` | a fractional polynomial, `powers` to choose them | 1 or 2 |
#'
#' A shape with more than one term prints a departure-from-linearity test and
#' no single slope, because no single number describes a curve. A one-term
#' shape that is not the straight line, `"log"` or a one-power `"fp"`, reports
#' its coefficient per unit of its own term rather than per unit of dose, and
#' says which.
#'
#' @section Reading the slope and its reference dose:
#' A linear fit reports one slope per unit of dose, measured against a
#' reference dose. A rate ratio of 1.07 per 10 g/day does not mean the exposed
#' run 1.07 times the rate. Each 10 g/day step multiplies the rate by about
#' 1.07, relative to the reference, and moving that reference moves every
#' predicted effect. easyDR prints the reference for that reason. A curve has
#' no single slope, so it prints none, and the predicted effects carry the
#' interpretation.
#'
#' @param data A data frame in long form, one row per dose category per study.
#'   It is never required to be reshaped, sorted or aggregated.
#' @param sm Summary measure, and so the scale every reported number is on.
#'   One of `"RR"`, `"OR"`, `"IRR"`, `"MD"`, `"SMD"`. It must suit the design:
#'   an odds ratio belongs to a case-control study, a rate ratio needs a
#'   person-time denominator.
#' @param shape The dose-response form, required: `"linear"`, `"log"`,
#'   `"quadratic"`, `"cubic"`, `"rcs"` for a restricted cubic spline, or
#'   `"fp"` for a fractional polynomial. There is no default. A shape the data
#'   cannot support stops the fit rather than falling back to a simpler one.
#'   Note that `"cubic"` is the third-degree polynomial in dose, which is a
#'   different model from the restricted cubic spline under `"rcs"`.
#' @param dose Name of the column holding the numeric dose you assigned to each
#'   category. Defaults to the column named `dose`. A dose is never inferred
#'   from a category range: that judgement is yours.
#' @param studlab Name of the study label column. Rows sharing a label are one
#'   study. Defaults to `studlab`.
#' @param event.e,n.e Names of the count columns: events, and the denominator.
#'   For incidence-rate data the denominator is person-time, and for `"MD"` or
#'   `"SMD"` it is the group size.
#' @param TE,lower,upper Names of the reported-effect columns: the effect on the
#'   reported scale with its confidence limits. When present this is the
#'   estimate that is pooled, because it is the adjusted estimate the authors
#'   intended. Supply counts alongside it whenever the table has them.
#' @param mean.e,sd.e Names of the mean and standard deviation columns, for
#'   `sm = "MD"` or `sm = "SMD"`.
#' @param ref Name of a logical column marking each study's reference category.
#'   When absent it is inferred, preferring a row whose effect is exactly 1,
#'   and the studies it was guessed for are printed.
#' @param design Name of a column naming each study's design: case-control,
#'   incidence rate, or cumulative incidence. When absent, a table of counts is
#'   treated as cumulative incidence and the assumption is printed.
#' @param method.tau Heterogeneity estimator. `"REML"` by default; `"ML"`,
#'   `"mm"` and `"fixed"` are also available.
#' @param proc `"2stage"` (default) fits a dose-response within each study and
#'   pools the study-specific slopes, which is what most published reviews
#'   report. `"1stage"` fits every study in a single model, which behaves
#'   better with few categories or sparse cells.
#' @param covariance Override the covariance approximation with `"gl"`, `"h"`
#'   or `"indep"`. `NULL`, the default, chooses from what the table supplies.
#' @param method.smd For `sm = "SMD"`: `"cohens"`, `"hedges"` or `"glass"`.
#' @param knots For `shape = "rcs"`: knot positions, or a single number giving
#'   how many to place. The default is 3 knots at the 10th, 50th and 90th
#'   percentiles of the non-zero doses, and whatever is used is printed.
#' @param powers For `shape = "fp"`: one power for a first-order fractional
#'   polynomial, such as `powers = -0.5`, or a pair for a second-order one.
#'   The default is the pair `c(-1, 0)`. A power of 0 is read as the logarithm,
#'   following Royston and Altman.
#' @param at Doses to report predicted effects at. `NULL` uses a five-point
#'   grid across the observed dose range.
#' @param grid Doses for the curve table on the result, which [easydr_curve()]
#'   draws. `NULL` uses 20 points across the observed range.
#' @param dose.unit A label such as `"g/day"`, printed beside the slope and the
#'   doses. It is a label only and changes no number.
#' @param level Confidence and prediction level. `0.95` by default.
#' @param random,common Which model to fit. The default fits random effects.
#'   Either `common = TRUE` or `random = FALSE` fits the common-effect model
#'   instead, as does `method.tau = "fixed"`. One model is fitted rather than
#'   two displayed, and the printed output names the one it fitted. A
#'   common-effect fit reports no prediction interval, because it assumes there
#'   is no between-study variation for one to describe.
#' @param subgroup The name of one column holding the subgroup each study
#'   belongs to, as a string, for example `subgroup = "design"`. The whole
#'   model is then fitted again inside each level, inheriting the shape, the
#'   pooling stage, the covariance approximation and the columns, so a level
#'   differs from the overall fit only in which studies it contains. Each level
#'   keeps its own between-study variance. A subgroup is a characteristic of a
#'   study, so the column must take one value per study. Nothing is tested
#'   between levels: [easydr_forest()] and [easydr_curve()] report them side by
#'   side and leave the comparison to you.
#' @param min.levels Exposure levels per study the viability check expects.
#'   Three is the conventional floor for a two-stage fit and the default.
#'
#' @return An object of class `easydr_result`. Print it for the whole answer;
#'   [summary()][summary.easydr_result()] adds the per-study contributions.
#'   Pass it to [easydr_curve()], [easydr_forest()], [easydr_loo()] or
#'   [easydr_funnel()].
#'
#' @references
#' Greenland S, Longnecker MP (1992). Methods for trend estimation from
#' summarized dose-response data, with applications to meta-analysis. American
#' Journal of Epidemiology 135, 1301-1309.
#'
#' Orsini N, Li R, Wolk A, Khudyakov P, Spiegelman D (2012). Meta-analysis for
#' linear and nonlinear dose-response relations. American Journal of
#' Epidemiology 175, 66-73.
#'
#' Crippa A, Orsini N (2016). Multivariate dose-response meta-analysis: the
#' dosresmeta R package. Journal of Statistical Software 72, 1-15.
#'
#' @examples
#' # Eight cohort studies of alcohol intake and colorectal cancer, one row per
#' # intake category: the published rate ratio with its limits, the cases, and
#' # the person-years.
#' head(alcohol_cancer)
#'
#' easydr.object <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
#'                         dose.unit = "g/day", at = c(10, 25, 50))
#' easydr.object
#'
#' # A spline on the same table. The shape is named, never chosen for you.
#' easydr(alcohol_cancer, sm = "RR", shape = "rcs", dose.unit = "g/day")
#'
#' # Continuous outcomes: means with standard deviations.
#' easydr(aripiprazole, sm = "SMD", shape = "linear", dose.unit = "mg/day")
#'
#' @export
easydr <- function(data,
                   sm = "RR",
                   shape,
                   dose = NULL,
                   studlab = NULL,
                   event.e = NULL,
                   n.e = NULL,
                   TE = NULL,
                   lower = NULL,
                   upper = NULL,
                   mean.e = NULL,
                   sd.e = NULL,
                   ref = NULL,
                   design = NULL,
                   method.tau = "REML",
                   proc = "2stage",
                   covariance = NULL,
                   method.smd = "cohens",
                   knots = NULL,
                   powers = NULL,
                   at = NULL,
                   grid = NULL,
                   dose.unit = NULL,
                   level = 0.95,
                   random = TRUE,
                   common = FALSE,
                   subgroup = NULL,
                   min.levels = 3L) {
  call <- match.call()
  sm <- match.arg(sm, c("RR", "OR", "IRR", "MD", "SMD"))

  # No default, and no guess. See the "Why you must name the shape" section.
  if (missing(shape)) {
    stop("easydr: name the `shape` you want fitted. There is no default, ",
         "because the shape is a claim about the biology and belongs in your ",
         "protocol with a reason beside it. Name one of:
",
         "  shape = \"linear\"      a straight line through the dose
",
         "  shape = \"log\"         a straight line through log dose
",
         "  shape = \"quadratic\"   dose and dose squared
",
         "  shape = \"cubic\"       dose, dose squared and dose cubed
",
         "  shape = \"rcs\"         a restricted cubic spline, see `knots`
",
         "  shape = \"fp\"          a fractional polynomial, see `powers`",
         call. = FALSE)
  }

  # An argument that belongs to the other kind of outcome is an error, not a
  # silently ignored name. Passing sd.e with sm = "RR" means the caller has
  # misread which measure they are fitting, and the fit would otherwise look
  # like it had used the column.
  cont <- sm %in% c("MD", "SMD")
  wrong <- if (cont) {
    c(event.e = "event.e", TE = "TE", lower = "lower", upper = "upper",
      design = "design", covariance = "covariance")
  } else {
    c(mean.e = "mean.e", sd.e = "sd.e")
  }
  given <- vapply(names(wrong), function(nm) !is.null(get(nm)), logical(1))
  if (any(given)) {
    stop("easydr: sm = \"", sm, "\" is a ",
         if (cont) "continuous" else "binary",
         " outcome, so ", paste(wrong[given], collapse = ", "),
         " ", if (sum(given) == 1L) "is" else "are",
         " not read. Drop ", if (sum(given) == 1L) "it" else "them",
         " or change sm.", call. = FALSE)
  }

  out <- if (cont) {
    easydrcont(data = data, dose = dose, studlab = studlab, sm = sm,
               mean.e = mean.e, sd.e = sd.e, n.e = n.e, ref = ref,
               method.tau = method.tau, proc = proc, method.smd = method.smd,
               shape = shape, knots = knots, powers = powers, grid = grid,
               at = at, dose.unit = dose.unit, level = level,
               random = random, common = common)
  } else {
    easydrbin(data = data, dose = dose, studlab = studlab, sm = sm,
              event.e = event.e, n.e = n.e, TE = TE, lower = lower,
              upper = upper, ref = ref, design = design,
              method.tau = method.tau, proc = proc, covariance = covariance,
              shape = shape, knots = knots, powers = powers, grid = grid,
              at = at, dose.unit = dose.unit, level = level,
              random = random, common = common)
  }

  # The two checks that belong either side of the model, folded in so nobody
  # has to remember to run them. The viability check reads the prepared frame
  # rather than the raw table, so it reports the reference category the fit
  # actually used instead of guessing at a zero dose.
  # The count columns are only worth checking when the fit used them. A table
  # of reported effects with no denominators carries no counts by design, and
  # the printed output already says the contrasts were treated as independent.
  counts_used <- out$provenance$intake %in% c("counts", "effect_counts")
  out$viability <- dr_viability(
    out$frame, study = "studlab", dose = "dose", ref = "is_ref",
    cases = if (counts_used) "cases" else NULL,
    n = if (counts_used) "n" else NULL,
    min_levels = min.levels)
  out$curve_label <- dr_curve_label(out)

  # The subgroup levels are fitted here rather than at drawing time, so the
  # fit a reader prints and the fit a figure draws are the same objects. The
  # refits inherit through dr_refit(), and they are not themselves given a
  # subgroup, so this does not recurse.
  if (!is.null(subgroup)) out$subgroup <- dr_subgroup_fits(out, subgroup)

  out$call <- call
  out
}

#' Name the fitted curve, on the scale the model was fitted on
#'
#' `curve_shape()` describes how the effect changes from one dose to the next,
#' so a monotone curve gets the same name on the log scale and the ratio scale.
#' The log scale is passed because that is the scale the coefficients live on,
#' and because its null value is that function's documented default.
#'
#' Fewer than three grid doses cannot be described, which happens only when a
#' caller passes a shorter `grid` than the default.
#'
#' @keywords internal
#' @noRd
dr_curve_label <- function(object) {
  tab <- object$curve
  if (is.null(tab) || nrow(tab) < 3L) return(NA_character_)
  est <- if (isTRUE(object$ratio)) log(tab$est) else tab$est
  if (anyNA(est) || any(!is.finite(est))) return(NA_character_)
  curve_shape(tab$dose, est)
}
