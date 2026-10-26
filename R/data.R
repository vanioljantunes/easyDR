#' Alcohol intake and colorectal cancer
#'
#' Eight prospective cohort studies of alcohol intake and colorectal cancer,
#' already extracted into the shape easydr reads: one row per intake category
#' per study, with the published rate ratio, its confidence limits, the case
#' count, the person-time denominator, and the reference category flag.
#'
#' Use it to see what your own extraction table should look like before you
#' call [easydr()].
#'
#' @format A data frame with 48 rows and 9 columns.
#' \describe{
#'   \item{studlab}{study label.}
#'   \item{dose}{assigned alcohol intake for the category, grams per day.}
#'   \item{TE}{published incidence rate ratio against the study reference.}
#'   \item{lower, upper}{95 percent confidence limits of `TE`. Missing on the
#'     reference row, which has no uncertainty of its own.}
#'   \item{event.e}{cases in the category.}
#'   \item{n.e}{person-years in the category.}
#'   \item{design}{study design, incidence rate for all eight studies.}
#'   \item{ref}{`TRUE` on each study's reference category.}
#' }
#'
#' @source Derived from the `alcohol_crc` dataset of the `dosresmeta` package
#'   (GPL-2), which reproduces the analysis in Orsini N, Li R, Wolk A,
#'   Khudyakov P, Spiegelman D (2012). Meta-analysis for linear and nonlinear
#'   dose-response relations. American Journal of Epidemiology 175, 66-73. The
#'   rate ratios and their limits are rounded to four decimals; the counts are
#'   unchanged. The build script is `data-raw/make-datasets.R`.
#'
#' @examples
#' head(alcohol_cancer)
#' easydr(alcohol_cancer, sm = "RR", shape = "linear", dose.unit = "g/day", at = c(10, 25, 50))
"alcohol_cancer"

#' Alcohol intake and cardiovascular disease
#'
#' Six studies of alcohol intake and cardiovascular disease, in the shape
#' easydr reads: one row per intake category per study, with the published
#' ratio, its confidence limits, the case count, the denominator, the study
#' design and the reference category flag.
#'
#' Two things make this the table to reach for when a result needs probing
#' rather than simply fitting.
#'
#' The studies disagree. Between-study variance is well away from zero, so the
#' pooled slope moves when a study is removed and a prediction interval has
#' something to describe. [alcohol_cancer] is the opposite case, agreeing
#' almost exactly, which makes it a poor table for a sensitivity analysis and a
#' good one for a first fit.
#'
#' The designs differ. Four studies are case-control and two report cumulative
#' incidence. Each keeps its own covariance approximation, `sm` names the scale
#' the reported ratios are pooled on, and the printed output says which designs
#' are present. Design is therefore a real study-level covariate here, which is
#' a study-level characteristic to compare the cohorts on.
#'
#' @format A data frame with 25 rows and 9 columns.
#' \describe{
#'   \item{studlab}{study label, the first author as the source names them.}
#'   \item{dose}{assigned alcohol intake for the category, grams per day.}
#'   \item{TE}{published ratio against the study reference.}
#'   \item{lower, upper}{95 percent confidence limits of `TE`. Missing on the
#'     reference row, which has no uncertainty of its own.}
#'   \item{event.e}{cases in the category.}
#'   \item{n.e}{participants in the category.}
#'   \item{design}{study design: case-control for four studies, cumulative
#'     incidence for two.}
#'   \item{ref}{`TRUE` on each study's reference category.}
#' }
#'
#' @source Derived from the `alcohol_cvd` dataset of the `dosresmeta` package
#'   (GPL-2). The design comes from that dataset's `type` column rather than
#'   being assumed, because it is what selects each study's covariance
#'   approximation. The ratios are not rounded, unlike [alcohol_cancer], so the
#'   examples here reproduce the engine exactly. The build script is
#'   `data-raw/make-datasets.R`.
#'
#' @examples
#' head(alcohol_cvd)
#' easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day", at = c(10, 25, 50))
"alcohol_cvd"

#' Aripiprazole dose and symptom change
#'
#' Five trials of aripiprazole reporting the mean change in a symptom score at
#' several doses, extracted into the shape easydr reads: one row per dose group
#' per trial, with the group mean, its standard deviation, the group size, and
#' the reference category flag.
#'
#' Use it as the worked example for [easydr()].
#'
#' @format A data frame with 18 rows and 6 columns.
#' \describe{
#'   \item{studlab}{trial label.}
#'   \item{dose}{assigned dose for the group, milligrams per day.}
#'   \item{mean.e}{mean change in the symptom score.}
#'   \item{sd.e}{standard deviation of the change.}
#'   \item{n.e}{number of participants in the group.}
#'   \item{ref}{`TRUE` on each trial's placebo group, which is its reference.}
#' }
#'
#' @source Derived from the `ari` dataset of the `dosresmeta` package (GPL-2),
#'   which reproduces the analysis in Crippa A, Orsini N (2016). Dose-response
#'   meta-analysis of differences in means. BMC Medical Research Methodology
#'   16, 91. Values are unchanged; only the column names and types differ. The
#'   build script is `data-raw/make-datasets.R`.
#'
#' @examples
#' head(aripiprazole)
#' easydr(aripiprazole, sm = "SMD", shape = "linear", dose.unit = "mg/day", at = c(5, 10, 20))
"aripiprazole"
