# easyDR

**Dose-response meta-analysis in one call, argued the way the `meta` package
argues. The book that teaches it is free.**

Sponsored by [Meta-Analysis Academy](https://metaanalysisacademy.com).

## Contents

| Section | Jump to |
|:---|:---|
| **[1 · The book](#1-the-book)**<br><sub>Free, six chapters, the whole method in the order of the work</sub> | [Read it](https://vanioljantunes.github.io/easyDR/) · [Chapters](#what-is-in-it) |
| **[2 · What this package is](#2-what-this-package-is)**<br><sub>One entry point, four things to do with what it returns</sub> | |
| **[3 · Install, data, and what <code>sm</code> decides](#3-install-data-and-what-sm-decides)**<br><sub>From a spreadsheet to a table easyDR reads</sub> | |
| **[4 · Fit it](#4-fit-it)**<br><sub><code>easydr()</code>, and why <code>shape</code> has no default</sub> | |
| **[5 · What to do with the fit](#5-what-to-do-with-the-fit)**<br><sub>Forest, leave-one-out, funnel</sub> | |
| **[6 · Validation](#6-validation-reproduce-a-published-analysis)**<br><sub>Paste it in, and check the numbers against this page</sub> | [Run it](#run-it) · [What you should see](#what-you-should-see) |

<sub>[References](#references) · [License](#license)</sub>

---

## 1. The book

The book in this repository teaches dose-response meta-analysis end to end,
and it is free to read:
**[open the book](https://vanioljantunes.github.io/easyDR/)**.

### What is in it

Six chapters, in the order the work actually happens:

| Chapter | What it covers |
|:---|:---|
| [Welcome](book/index.qmd) | What the book is for, and who it is for |
| [1 · What it is, goals and examples](book/01-what-is-it.qmd) | The question a dose-response answers, the five curve shapes, worked examples |
| [2 · Viability, search strategy and data extraction](book/02-viability.qmd) | Whether your studies can support one, and the extraction sheet that follows |
| [3 · Registering on PROSPERO and risk of bias](book/03-prospero.qmd) | The protocol, and appraising the studies you found |
| [4 · Statistics theory and code](book/04-statistics.qmd) | The two-stage and one-stage models, splines, and the R to fit them |
| [5 · Checklist, writing and interpretation](book/05-checklist-writing.qmd) | Reporting it, and what the curve does and does not license |

The source is in [`book/`](book) as Quarto files; `quarto render book` builds
it locally, after `R CMD INSTALL .` so the chunks can call the package.

---

## 2. What this package is

![One call and four things to do with it: easyDR has no formula argument, and folds the viability check and the curve label into the fit](man/figures/overview.png)

<details open>
<summary><code>{r the-one-call}</code> &mdash; the fit, the check that preceded it, the name of the curve</summary>

```r
library(easyDR)

fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
              dose.unit = "g/day", at = c(10, 25, 50))
fit                 # the whole answer, printed at once
summary(fit)        # the same, plus the coefficient and the per-study rows
fit$viability       # the check that runs before the model
fit$curve_label     # the name the curve earned: "monotonic increasing"
```

</details>

---

## 3. Install, data, and what `sm` decides

![Install easyDR, lay the data out one row per exposure level per study, and let sm choose which columns are read](man/figures/data.png)

<details open>
<summary><code>{r install-and-look}</code> &mdash; install it, then look at the two bundled tables</summary>

```r
# install.packages("remotes")
remotes::install_github("vanioljantunes/easyDR")
library(easyDR)

head(alcohol_cancer)   # a binary outcome, one row per intake category
head(aripiprazole)     # a continuous one: mean.e, sd.e, n.e
```

</details>

---

## 4. Fit it

![easydr() in one call, the five shapes you can name, and easydr_curve() with its arguments](man/figures/api.png)

<details open>
<summary><code>{r fit-and-draw}</code> &mdash; the call, the refusal when <code>shape</code> is missing, and the curve</summary>

```r
fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
              dose.unit = "g/day", at = c(10, 25, 50))

# shape has no default. Leave it out and the error lists what you may write:
# "linear", "log", "quadratic", "rcs", "fp".
easydr(alcohol_cancer, sm = "RR")

# the curve, with every argument the picture used
easydr_curve(fit,
             main  = "Alcohol and colorectal cancer",
             xlab  = "Alcohol intake (g/day)",
             ylab  = "Risk ratio",
             ylim  = c(0.9, 1.8),
             round = 2,
             prediction = TRUE)
```

</details>

---

## 5. What to do with the fit

![easydr_forest, easydr_loo and easydr_funnel, all taking the object easydr returned](man/figures/plots.png)

<details open>
<summary><code>{r three-questions}</code> &mdash; what each study gave, whether one is carrying it, whether the small ones differ</summary>

```r
cvd <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")

# what each study contributed
easydr_forest(fit)$table                 # draws it, and the table prints itself
easydr_forest(rcs, at = 25)              # a curve has no single trend

# is one study carrying it
easydr_loo(cvd, quiet = TRUE)            # draws the curves and the table
easydr_loo(cvd, quiet = TRUE, plot = FALSE)$table   # the table alone

# are the small studies different
easydr_funnel(fit, labels = TRUE)$egger  # Egger's test, printed as a test
```

</details>

---

## 6. Validation: reproduce a published analysis

Everything above is a picture. This section is not. Paste it in and you should
get the numbers printed here, on your own machine.

The example is the one the book uses: eight prospective cohorts of alcohol
intake and colorectal cancer, analysed in Orsini et al. (2012) and bundled as
`alcohol_cancer`.

### Run it

<details open>
<summary><code>{r setup}</code> &mdash; install and attach</summary>

```r
# install.packages("remotes")
remotes::install_github("vanioljantunes/easyDR")
library(easyDR)
```

</details>

<details open>
<summary><code>{r step-1-fit}</code> &mdash; one call, and the whole answer prints</summary>

```r
fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
              dose.unit = "g/day", at = c(10, 25, 50))
fit
```

```text
Dose-response meta-analysis (risk ratio)
Pooling: two-stage random effects, reml. Covariance: Greenland and Longnecker.

Per 1 g/day increase:  RR 1.006  (95% CI 1.004 to 1.009)
Reference dose: 0 g/day
Curve shape: monotonic increasing

Departure from linearity: not estimable
  a linear fit has one term, so there is nothing to test it against
  8 of 8 studies have three or more dose levels

Heterogeneity: Q = 4.78 (df = 7), p = 0.687; I-squared = 0%

Predicted RR
 dose (g/day) RR        95% CI         95% PI
  0           reference reference      reference
 10           1.066     1.039 to 1.095 1.033 to 1.101
 25           1.175     1.099 to 1.255 1.084 to 1.272
 50           1.380     1.208 to 1.575 1.176 to 1.619

8 studies, 48 dose categories, dose range 0 to 60.6 g/day
Viability: 8 of 8 studies carry 3 or more dose levels.
Defaults applied:
  - intake: reported effect with limits, counts used for the covariance
```

</details>

No `formula` argument, and no `logrr` column. `sm = "RR"` tells easyDR to read
the counts and the reported risk ratios, then derive the contrasts and
reconstruct the within-study covariances. The risk ratios at 10, 25 and 50
g/day are in that same print, each one measured against the reference dose,
which is why the reference row stays in the table.

<details open>
<summary><code>{r step-2-summary}</code> &mdash; the coefficient, the fit, and what each cohort gave</summary>

```r
summary(fit)
```

```text
[everything the fit printed above, then:]

Goodness of fit: deviance D = 48.21 (df = 39), p = 0.148
  R-squared 0.32, adjusted 0.302
  deviance of the fitted model against a saturated fit

Pooled trend on the log scale, per 1 g/day
  b = 0.0064373,  SE = 0.0013521,  p < 0.001
  no intercept: every contrast runs against the study's own reference category

Per-study contribution
 study levels dose range design slope (log scale) SE
 atm   6      0 to 58.43 ir      0.007141         0.003912
 hpm   6      0 to 59.4  ir      0.009169         0.003044
 io    6      0 to 60.6  ir      0.002359         0.004106
 n1    6      0 to 58.4  ir     -0.002127         0.006536
 n2    6      0 to 57.6  ir      0.004905         0.004073
 nt    6      0 to 50.33 ir      0.012630         0.005995
 ntm   6      0 to 56.23 ir      0.006131         0.003314
 nym   6      0 to 57    ir      0.006880         0.002990
```

</details>

The coefficient on the log scale is the number to quote in a methods section,
and the one to check an install against. `summary()` reports it, so nothing has
to be pulled out of the fitted object and rounded by hand.

<details open>
<summary><code>{r step-3-forest}</code> &mdash; the rows behind the forest plot, with weights</summary>

```r
easydr_forest(fit)$table     # draws the forest, and the table prints itself
```

```text
 study levels dose range trend      SE        weight %
 atm   6      0 to 58.4   0.0071406 0.0039123 11.94
 hpm   6      0 to 59.4   0.0091691 0.0030437 19.73
 io    6      0 to 60.6   0.0023590 0.0041061 10.84
 n1    6      0 to 58.4  -0.0021265 0.0065361  4.28
 n2    6      0 to 57.6   0.0049053 0.0040725 11.02
 nt    6      0 to 50.3   0.0126300 0.0059954  5.09
 ntm   6      0 to 56.2   0.0061313 0.0033139 16.65
 nym   6      0 to 57     0.0068797 0.0029901 20.45
```

</details>

Each plot returns the table it drew. These rows are the engine's own
first-stage estimates. The printing rounds, the stored values do not, so
`as.data.frame()` gives them back untouched.

<details open>
<summary><code>{r step-4-funnel}</code> &mdash; small-study effects</summary>

```r
easydr_funnel(fit)$egger
```

```text
Egger's test for small-study effects
  intercept -1.0879 (SE 1.2204)
  t = -0.89 on 6 df, p = 0.407
  8 studies: fewer than 10, so the test has little power and a null result means very little
```

</details>

Egger's regression. Deeks' test belongs to diagnostic test accuracy, so it does
not apply here. The figure prints the same power warning above its panel.

<details open>
<summary><code>{r step-5-checks}</code> &mdash; the two checks easyDR ran for you</summary>

```r
fit$curve_label
fit$viability
```

```text
[1] "monotonic increasing"
Dose-response viability

 study levels has_reference two_stage_ok
 atm   6      TRUE          TRUE
 hpm   6      TRUE          TRUE
 io    6      TRUE          TRUE
 n1    6      TRUE          TRUE
 n2    6      TRUE          TRUE
 nt    6      TRUE          TRUE
 ntm   6      TRUE          TRUE
 nym   6      TRUE          TRUE

No blocking problems found.
```

</details>


### What you should see

| Quantity | Value |
|---|---|
| Studies, exposure levels each | 8, 6 |
| Dose slope, log scale, per g/day | 0.0064373 |
| Standard error | 0.0013521 |
| Risk ratio at 10 g/day | 1.0665 (1.0386 to 1.0951) |
| Risk ratio at 25 g/day | 1.1746 (1.0993 to 1.2551) |
| Risk ratio at 50 g/day | 1.3797 (1.2085 to 1.5752) |
| Cochran Q | 4.78 on 7 df, p = 0.687 |
| I-squared | 0% |
| Largest contribution | nym, 20.4% |
| Egger intercept | -1.0879 (SE 1.2204), p = 0.407 |
| Curve shape | monotonic increasing |

Captured on R 4.6.0 with `dosresmeta` 2.2.0. A different BLAS or engine
version can move the last digit or two, but the rate ratios should agree to
three decimals. If they do not, something in your install differs, and the rest
of the book will differ with it. Worth finding out before you trust a curve of
your own.

---

## References

Greenland S, Longnecker MP (1992). Methods for trend estimation from summarized
dose-response data, with applications to meta-analysis. *American Journal of
Epidemiology* 135(11), 1301-1309.

Crippa A, Discacciati A, Bottai M, Spiegelman D, Orsini N (2019). One-stage
dose-response meta-analysis for aggregated data. *Statistical Methods in
Medical Research* 28(5), 1579-1596.

Orsini N, Li R, Wolk A, Khudyakov P, Spiegelman D (2012). Meta-analysis for
linear and nonlinear dose-response relations: examples, an evaluation of
approximations, and software. *American Journal of Epidemiology* 175(1), 66-73.

Crippa A, Orsini N (2016). Multivariate dose-response meta-analysis: the
dosresmeta R package. *Journal of Statistical Software* 72(1), 1-15.

Egger M, Davey Smith G, Schneider M, Minder C (1997). Bias in meta-analysis
detected by a simple, graphical test. *BMJ* 315(7109), 629-634.

The book's full bibliography is in [`book/references.bib`](book/references.bib).

---

## License

GPL-2 for the code in `R/`, `man/` and `tests/`. See [`LICENSE`](LICENSE).
easyDR imports `dosresmeta`, which is GPL-2, so easyDR has to be GPL-2 too.

CC BY-NC-SA 4.0 for the book text, figures and examples in `book/`. See
[`LICENSE-text.md`](LICENSE-text.md).

---

## Related

Part of the MetaRank 4 catalog from the [Meta-Analysis
Academy](https://metaanalysisacademy.com): one design, taught end to end, per
repo. The book is free. The video course, templates and support are the paid
half.

The four pictures above are built from [`tutorial/`](tutorial): `Rscript
tutorial/figures.R` draws the figures and captures the console blocks, then
`node tutorial/render.mjs <page>` for each of `overview`, `data`, `api` and
`plots`. The layout follows easyDTA, and the funnel follows `metafor`.

Questions and corrections:
[issues](https://github.com/vanioljantunes/easyDR/issues).
