# easyDR

**Dose-response meta-analysis, made teachable: a free book, and the two checks
that belong either side of the model.**

## Contents

| Section | Jump to |
|:---|:---|
| **[1 · The book](#1-the-book)**<br><sub>Free, six chapters, the whole method in the order of the work</sub> | [Read it](https://vanioljantunes.github.io/easyDR/) · [Chapters](#what-is-in-it) |
| **[2 · What this package is](#2-what-this-package-is)**<br><sub>Two functions, and why it is only two</sub> | |
| **[3 · How to install it](#3-how-to-install-it)**<br><sub>One line from GitHub</sub> | |
| **[4 · Data preparation](#4-data-preparation)**<br><sub>One row per exposure level per study</sub> | [Columns](#the-columns) · [Doses](#doses-must-already-be-comparable) |
| **[5 · Before the model](#5-before-the-model-dr_viability)**<br><sub>Can these studies support a dose-response at all?</sub> | [`dr_viability()`](#5-before-the-model-dr_viability) · [One-stage](#when-the-answer-is-no-it-may-be-an-argument-for-one-stage) |
| **[6 · The model itself](#6-the-model-itself-dosresmeta)**<br><sub>Delegated to dosresmeta, with one label to watch</sub> | [Fit](#6-the-model-itself-dosresmeta) · [Gotcha](#the-coefficient-is-not-an-intercept) |
| **[7 · After the model](#7-after-the-model-curve_shape)**<br><sub>Name the curve in the words the literature uses</sub> | [`curve_shape()`](#7-after-the-model-curve_shape) · [Scales](#ratio-scales) |
| **[8 · Build the book locally](#8-build-the-book-locally)**<br><sub>Quarto, and why the package has to be installed first</sub> | |

<sub>[References](#references) · [License](#license)</sub>

---

## 1. The book

The book is the main thing in this repository. It is free, it is online, and it
teaches dose-response meta-analysis end to end:

### **[Read it: vanioljantunes.github.io/easyDR](https://vanioljantunes.github.io/easyDR/)**

Most meta-analysis teaching stops at a single pooled estimate, which answers
"treated versus untreated" and nothing else. A dose-response meta-analysis asks
a different question: how does the effect change as the dose changes? The answer
is a curve, not a number, and almost everything that makes it hard sits around
the model rather than in it.

The source is in [`book/`](book) in this same folder, as Quarto `.qmd` files.
The rendered site is what the link above serves.

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

The book text is CC BY-NC-SA 4.0. The package code is MIT. See
[License](#license).

---

## 2. What this package is

`easyDR` is the package the book leans on, and it is deliberately small. It
carries the two steps students repeatedly get wrong, and leaves the modelling
to [dosresmeta](https://github.com/alecri/dosresmeta) by Crippa and Orsini,
which is validated, peer reviewed, and already does that job well.

| Function | Does |
|---|---|
| `dr_viability()` | Checks whether a dataset can support a dose-response model at all, before any fitting |
| `curve_shape()` | Names the shape of a fitted curve in the vocabulary the literature uses |

**What you get:**

- **A check that runs before the model.** The commonest reason a dose-response
  meta-analysis fails is not the model; it is that the included studies never
  could have supported one. `dr_viability()` reports that per study, with the
  reason, rather than letting a fit fail later for a reason that reads as
  numerical.
- **A name for the curve.** `curve_shape()` returns one of five labels, so the
  results paragraph and the figure agree with each other, and with the words
  the field already uses.
- **No estimator of its own.** Every estimate in the book comes from
  `dosresmeta`. This package prepares and interprets; it does not compute a
  pooled slope, and it does not wrap `dosresmeta` in a layer you then have to
  learn.

**What you do not get**, so you are not surprised: no fitting function, no
plotting, no bundled datasets. The fitting examples below use `dosresmeta` and
the data `dosresmeta` ships.

---

## 3. How to install it

```r
# install.packages("remotes")
remotes::install_github("vanioljantunes/easyDR")
```

The package itself needs only base R (>= 4.1). To run the fitting section of
this page, and the code in the book's statistics chapter, install the engine
too:

```r
install.packages("dosresmeta")
```

It is in `Suggests` rather than `Imports`, because the two functions here do not
call it.

---

## 4. Data preparation

One row per exposure level per study, in long form. This is the layout the
book's extraction sheet produces and the layout `dosresmeta` expects, so
nothing has to be reshaped between them.

### The columns

| Column | Holds |
|---|---|
| `study` | Study label. Rows sharing a label are one study |
| `dose` | The exposure assigned to that category, in a unit shared across studies |
| `cases` | Events in the category |
| `n` | The denominator: group size, or person-time for rate data |

Each study needs a reference category, conventionally the row with the lowest
dose, against which its other categories are compared.

### Doses must already be comparable

`dr_viability()` checks structure, not units. If one study reports grams per day
and another cups per week, the frame will pass every check and the pooled curve
will be meaningless. Converting to a common unit is extraction work, and the
book's second chapter covers it.

---

## 5. Before the model: `dr_viability()`

Two-stage methods need at least three exposure levels per study. A study with
two points carries no information about curvature, and its within-study
covariance cannot be reconstructed (Greenland and Longnecker, 1992). That is
the floor `dr_viability()` checks, and three is its default.

```r
library(easyDR)

good <- data.frame(
  study = rep(c("Allen 2011", "Brown 2014", "Chen 2019"), each = 4),
  dose  = rep(c(0, 5, 15, 30), 3),
  cases = c(12, 19, 28, 41, 8, 14, 22, 33, 15, 21, 30, 44),
  n     = rep(c(200, 190, 185, 180), 3)
)

dr_viability(good)
```

```text
Dose-response viability

       study levels has_reference two_stage_ok
1 Allen 2011      4          TRUE         TRUE
2 Brown 2014      4          TRUE         TRUE
3  Chen 2019      4          TRUE         TRUE

No blocking problems found.
```

A dataset with a two-level study and a missing count reports both, and says what
each one costs you:

```r
bad <- data.frame(
  study = c("Allen 2011", "Allen 2011", "Allen 2011",
            "Brown 2014", "Brown 2014",
            "Chen 2019", "Chen 2019", "Chen 2019"),
  dose  = c(0, 10, 20, 0, 15, 0, 8, 16),
  cases = c(10, 20, 31, 5, 9, 11, NA, 25),
  n     = c(100, 100, 100, 50, 50, 120, 120, 120)
)

dr_viability(bad)
```

```text
Dose-response viability

       study levels has_reference two_stage_ok
1 Allen 2011      3          TRUE         TRUE
2 Brown 2014      2          TRUE        FALSE
3  Chen 2019      3          TRUE         TRUE

Problems:
 - 1 study(ies) with fewer than 3 exposure levels (Brown 2014): a two-stage fit drops these. Consider one-stage.
 - missing values in 'cases': the within-study covariance cannot be reconstructed without it.
```

The object is a plain list, so the per-study table is there to use:

```r
v <- dr_viability(good)

v$ok        # TRUE when nothing blocking was found
v$summary   # one row per study: levels, has_reference, two_stage_ok
v$problems  # character vector, empty when ok
```

Column names other than the defaults are given as arguments, which is how you
point it at a frame written for the engine:

```r
library(dosresmeta)
data(alcohol_crc)

dr_viability(alcohol_crc, study = "id", dose = "dose",
             cases = "cases", n = "peryears")
```

```text
Dose-response viability

  study levels has_reference two_stage_ok
1   atm      6          TRUE         TRUE
2   hpm      6          TRUE         TRUE
3    io      6          TRUE         TRUE
4    n1      6          TRUE         TRUE
5    n2      6          TRUE         TRUE
6    nt      6          TRUE         TRUE
7   ntm      6          TRUE         TRUE
8   nym      6          TRUE         TRUE

No blocking problems found.
```

### When the answer is no, it may be an argument for one-stage

A dataset that fails the two-stage check is not useless. One-stage methods keep
studies with two exposure levels, and that is precisely their selling point
(Crippa et al., 2019). Lowering the floor shows what a one-stage model would
have to work with:

```r
dr_viability(bad, min_levels = 2)
```

```text
Dose-response viability

       study levels has_reference two_stage_ok
1 Allen 2011      3          TRUE         TRUE
2 Brown 2014      2          TRUE         TRUE
3  Chen 2019      3          TRUE         TRUE

Problems:
 - missing values in 'cases': the within-study covariance cannot be reconstructed without it.
```

The missing count is still reported, because no choice of model recovers it.

---

## 6. The model itself: `dosresmeta`

This package does not fit anything. The book's statistics chapter fits with
`dosresmeta`, and so does this page:

```r
library(dosresmeta)
data(alcohol_crc)

fit <- dosresmeta(formula = logrr ~ dose, type = type, id = id,
                  se = se, cases = cases, n = peryears,
                  data = alcohol_crc)

summary(fit)$coefficients
```

```text
               Estimate  Std. Error        z     Pr(>|z|)    95%ci.lb    95%ci.ub
(Intercept) 0.006437538 0.001352042 4.761344 1.923078e-06 0.003787584 0.009087492
```

### The coefficient is not an intercept

Look at that row name. The model has no intercept: it is suppressed, because
every contrast is against the study's own reference category. So the number
labelled `(Intercept)` is the **dose slope**, on the log scale, per unit of
dose. Here it is 0.006438 per gram per day, which exponentiates to a rate ratio
of 1.0065 per gram per day.

Reporting it as an intercept is a mistake the label invites, and it is worth
checking in any output you copy.

---

## 7. After the model: `curve_shape()`

A fitted curve has to be described in words at some point, and the words should
be the ones the literature uses. `curve_shape()` takes the predicted effects
and returns one of five labels:

```r
grid <- data.frame(dose = seq(0, 60, by = 10))
predicted <- predict(fit, newdata = grid, expo = FALSE)

curve_shape(grid$dose, predicted$pred)
#> [1] "positive association"
```

The five labels, and the behaviour that earns each one:

```r
curve_shape(0:4, c(0,  0.20,  0.40,  0.60,  0.80))  # rises throughout
curve_shape(0:4, c(0,  0.50,  0.80,  0.85,  0.86))  # rises, then flattens
curve_shape(0:4, c(0,  0.50,  0.90,  0.40, -0.10))  # rises, then turns down
curve_shape(0:4, c(0, -0.20, -0.40, -0.60, -0.80))  # falls throughout
curve_shape(0:4, c(0,  0.01,  0.00,  0.01,  0.00))  # flat
```

```text
rises throughout         -> positive association
rises then flattens      -> Emax
rises then turns down    -> non-monotonic
falls throughout         -> inverse association
flat                     -> no association
```

`Emax` is judged against the curve's own steepest step rather than its overall
range, which is what separates a genuine plateau from a curve that merely
decelerates.

### Ratio scales

The labels describe how the effect **changes** from one dose to the next, not
where the curve sits. So the answer is the same whether you pass log effects or
the ratios themselves:

```r
rr <- c(1.00, 1.22, 1.49, 1.82, 2.23)

curve_shape(0:4, log(rr))   # "positive association"
curve_shape(0:4, rr)        # "positive association"
```

That is the behaviour you want from a shape: a dose-response is about movement
with dose, and a curve that sits high but never moves has no dose-response
whatever its level.

Two arguments adjust where the thresholds fall, for when the defaults disagree
with your reading of a curve: `tol`, the fraction of the total range below which
one step counts as flat, and `plateau_frac`, how small the closing steps must be
to earn `Emax`.

---

## 8. Build the book locally

```bash
quarto render book
```

The rendered site lands in `book/_book`. Rendering executes the R chunks, and
those chunks call this package, so install it first:

```bash
R CMD INSTALL .
```

The statistics chapter also needs `dosresmeta`.

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

The book's full bibliography is in [`book/references.bib`](book/references.bib).

---

## License

MIT for the code in `R/`, `man/` and `tests/`. See [`LICENSE`](LICENSE).

CC BY-NC-SA 4.0 for the book text, figures and examples in `book/`. See
[`LICENSE-text.md`](LICENSE-text.md).

---

## Related

Part of the MetaRank 4 catalog from the [Meta-Analysis
Academy](https://metahub.metaanalysisacademy.com): one design, taught end to
end, per repo. The book is free. The video course, templates and support are the
paid half.

Questions and corrections:
[issues](https://github.com/vanioljantunes/easyDR/issues).
