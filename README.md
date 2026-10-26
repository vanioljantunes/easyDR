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

<details>
<summary><b>The one call, and what it returns</b>: <code>{r what-this-is}</code></summary>

```r
# One entry point. Paste the whole block into R.
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

<details>
<summary><b>Install it, then look at the data</b>: <code>{r install-and-data}</code></summary>

```r
# Paste the whole block into R.
# install.packages("remotes")
remotes::install_github("vanioljantunes/easyDR")

library(easyDR)

head(alcohol_cancer)   # a binary outcome, one row per intake category
head(aripiprazole)     # a continuous one: mean.e, sd.e, n.e

# sm decides which columns are read, how the within-study covariance is
# reconstructed, and whether the effect is drawn on a log axis.
```

</details>

---

## 4. Fit it

![easydr() in one call, the five shapes you can name, and easydr_curve() with its arguments](man/figures/api.png)

<details>
<summary><b>The call, the refusal, and the curve</b>: <code>{r fit-and-draw}</code></summary>

```r
# Paste the whole block into R.
library(easyDR)

# at = names the doses to report effects at. The curve marks and ticks them.
fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
              dose.unit = "g/day", at = c(10, 25, 50))
fit

# shape has no default. Leave it out and the error lists what you may write:
# "linear", "log", "quadratic", "rcs", "fp".
try(easydr(alcohol_cancer, sm = "RR"))

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

<details>
<summary><b>Forest, leave-one-out, funnel</b>: <code>{r three-questions}</code></summary>

```r
# Paste the whole block into R. Three questions, three calls, each taking the
# object easydr() returned.
library(easyDR)

fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
              dose.unit = "g/day", at = c(10, 25, 50))
rcs <- easydr(alcohol_cancer, sm = "RR", shape = "rcs", dose.unit = "g/day")
cvd <- easydr(alcohol_cvd,    sm = "RR", shape = "linear", dose.unit = "g/day")

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

Everything above is a picture. This section is not. Run the script and you
should get the values printed below, on your own machine.

The example is the one the book uses: eight prospective cohorts of alcohol
intake and colorectal cancer, analysed in Orsini et al. (2012) and bundled as
`alcohol_cancer`.

### Run it yourself

Paste the script into a fresh R session, or run it from the package root with
`Rscript validation/reproduce.R`. It needs easyDR and nothing else: no data
files, and it runs in a couple of seconds.

<details>
<summary><b>The whole analysis</b>: <code>validation/reproduce.R</code></summary>

```r
# easyDR: reproduce the published alcohol and colorectal cancer analysis.
# Eight prospective cohorts, analysed in Orsini et al. (2012) and bundled as
# alcohol_cancer. Paste the whole file into R; no data files are needed.
# install.packages("remotes"); remotes::install_github("vanioljantunes/easyDR")

library(easyDR)

# One call. No formula argument and no logrr column: sm = "RR" tells easyDR to
# read the counts and the reported risk ratios, derive the contrasts, and
# reconstruct the within-study covariances.
fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
              dose.unit = "g/day", at = c(10, 25, 50))
fit

# The coefficient on the log scale, the goodness of fit, and what each cohort
# contributed. This is the number to quote in a methods section.
summary(fit)

# The rows behind the forest plot, with weights. The figure is drawn and the
# table prints itself; the stored values keep full precision.
easydr_forest(fit)$table

# Small-study effects. Egger's regression, printed as a test result with its
# power warning. Deeks' test belongs to diagnostic test accuracy.
easydr_funnel(fit)$egger

# The two checks easyDR ran before and after the model.
fit$curve_label
fit$viability

# The curve, with every argument the README figure used.
easydr_curve(fit,
             main  = "Alcohol and colorectal cancer",
             xlab  = "Alcohol intake (g/day)",
             ylab  = "Risk ratio",
             round = 2)
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
