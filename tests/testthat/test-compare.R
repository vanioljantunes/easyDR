# Comparing shapes, and comparing knot placements.
#
# Both follow Crippa's published scripts, so what is tested is that they do
# what those scripts do: every candidate fitted on the same table, ranked by
# an information criterion that can compare models which are not nested, and
# every combination of the named percentiles tried.

test_that("every candidate shape is fitted on the same table", {
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  cs <- easydr_compare_shapes(f, shapes = c("linear", "quadratic", "rcs"),
                              plot = FALSE)
  expect_s3_class(cs, "easydr_comparison")
  expect_setequal(cs$table$Shape, c("linear", "quadratic", "rcs"))
  expect_equal(nrow(cs$table), 3L)
  for (fit in cs$fits) expect_s3_class(fit, "easydr_result")
})

test_that("the table is ranked by AIC", {
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  cs <- easydr_compare_shapes(f, shapes = c("linear", "quadratic", "rcs"),
                              plot = FALSE)
  expect_false(is.unsorted(cs$table$AIC))
})

test_that("the fit handed in is reused rather than refitted", {
  f <- easydr(aripiprazole, sm = "MD", shape = "quadratic",
              dose.unit = "mg/day")
  cs <- easydr_compare_shapes(f, shapes = c("linear", "quadratic"),
                              plot = FALSE)
  expect_equal(coef(cs$fits$quadratic$fit), coef(f$fit))
})

test_that("a shape that cannot be fitted is named, not dropped silently", {
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  cs <- easydr_compare_shapes(f, shapes = c("linear", "quadratic", "rcs"),
                              plot = FALSE)
  expect_true(is.character(cs$failed))
})

test_that("comparing needs at least two shapes, and one fit", {
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  expect_error(easydr_compare_shapes(f, shapes = "rcs"), "at least two")
  m <- easydr(aripiprazole, sm = "MD", shape = c("linear", "rcs"),
              dose.unit = "mg/day")
  expect_error(easydr_compare_shapes(m), "holds 2")
})

test_that("knot sensitivity tries every combination of the percentiles", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  ks <- easydr_knot_sensitivity(f, plot = FALSE)
  expect_s3_class(ks, "easydr_knots")
  # Five percentiles taken three at a time is ten combinations, which is what
  # Crippa's script produces.
  expect_equal(nrow(ks$table), 10L)
  expect_equal(length(ks$knots), 10L)
  for (k in ks$knots) expect_length(k, 3L)
})

test_that("the knots come from the dose distribution", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  ks <- easydr_knot_sensitivity(f, plot = FALSE)
  want <- sort(unique(round(stats::quantile(f$frame$dose,
                                            c(.1, .25, .5, .75, .9),
                                            names = FALSE), 2)))
  expect_true(all(unlist(ks$knots) %in% want))
})

test_that("the percentiles and the knot count can both be changed", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  # Four percentiles taken three at a time is four combinations.
  ks <- easydr_knot_sensitivity(f, percentiles = c(.1, .4, .6, .9),
                                nknots = 3, plot = FALSE)
  expect_equal(nrow(ks$table), 4L)
})

test_that("one combination is refused as too few, not as a failure to fit", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  expect_error(
    easydr_knot_sensitivity(f, percentiles = c(.1, .5, .9), nknots = 3),
    "only one set")
})

test_that("knot sensitivity refuses what it cannot do", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  expect_error(easydr_knot_sensitivity(f, percentiles = c(0.5, 0.9)),
               "at least `nknots`")
  expect_error(easydr_knot_sensitivity(f, percentiles = c(0, 0.5, 0.9)),
               "between 0 and 1")
})

test_that("both workflows write a figure when asked", {
  dir <- withr::local_tempdir()
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  easydr_compare_shapes(f, shapes = c("linear", "quadratic"),
                        file.name = "cmp", path = dir)
  v <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  easydr_knot_sensitivity(v, percentiles = c(.1, .4, .6, .9), nknots = 3,
                          file.name = "kn", path = dir)
  expect_true(file.exists(file.path(dir, "cmp.png")))
  expect_true(file.exists(file.path(dir, "kn.png")))
})
