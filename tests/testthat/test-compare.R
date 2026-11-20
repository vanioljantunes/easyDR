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

test_that("AIC compares shapes on one-stage likelihoods, whatever the dose unit", {
  # A two-stage likelihood is of each shape's own study coefficients, so it
  # moves with the dose unit and reverses rankings. The one-stage refit is of
  # the reported log relative risks, so rescaling the dose leaves it alone.
  shapes <- c("linear", "quadratic", "rcs")
  f1 <- easydr(alcohol_cvd, shape = "rcs", dose.unit = "g/day")
  d10 <- alcohol_cvd
  d10$dose <- d10$dose / 10
  f2 <- easydr(d10, shape = "rcs", dose.unit = "10 g/day")
  a1 <- easydr_compare_shapes(f1, shapes = shapes, plot = FALSE)$table
  a2 <- easydr_compare_shapes(f2, shapes = shapes, plot = FALSE)$table
  expect_equal(a1$AIC[match(shapes, a1$Shape)],
               a2$AIC[match(shapes, a2$Shape)], tolerance = 0.2)
  # Crippa and colleagues' one-stage linear fit of these data: AIC 32.67.
  expect_equal(a1$AIC[a1$Shape == "linear"], 32.7, tolerance = 0.05)
})

test_that("a two-stage shape leaves out studies with too few doses", {
  d <- alcohol_cvd
  d <- d[!(d$studlab == "Bianchi" & d$dose > 10), ]
  expect_message(q <- easydr(d, shape = "quadratic", dose.unit = "g/day"),
                 "left out 1 study")
  expect_equal(q$counts$studies, 5L)
  expect_identical(q$provenance$dropped$studies, "Bianchi")
  expect_output(print(q), "Left out of this shape")
  l <- easydr(d, shape = "linear", dose.unit = "g/day")
  expect_equal(l$counts$studies, 6L)
  cs <- suppressMessages(easydr_compare_shapes(
    l, shapes = c("linear", "quadratic"), plot = FALSE))$table
  expect_equal(cs$Studies[match(c("linear", "quadratic"), cs$Shape)], c(6L, 5L))
})

test_that("each test column can be switched on or off", {
  tab <- data.frame(Shape = "rcs", Studies = 6L, `p nonlinearity` = "0.06",
                    Deviance = 39.4, `p deviance` = "0.002", `R-squared` = 0.16,
                    `R-squared adj` = 0.06, AIC = 28.8, `Delta AIC` = 0.9,
                    BIC = 33.5, check.names = FALSE)
  expect_named(dr_test_columns(tab),
               c("Shape", "Studies", "p nonlinearity", "Deviance", "p deviance",
                 "Delta AIC"))
  all_on <- dr_test_columns(tab, r2 = T, r2adj = T, aic = T, bic = T)
  expect_named(all_on, names(tab))
  expect_named(dr_test_columns(tab, wald = F, deviance = F, delta.aic = F),
               c("Shape", "Studies"))
  expect_equal(dr_shape_title(c("rcs", "fp", "hyperbolic_emax", "other")),
               c("Restricted cubic spline", "Fractional polynomial",
                 "Hyperbolic Emax", "other"))
})

test_that("the curve's footer carries one line per test switched on", {
  f <- easydr(alcohol_cvd, shape = "rcs", dose.unit = "g/day")
  expect_length(dr_test_notes(f), 1L)
  n <- dr_test_notes(f, wald = FALSE, deviance = TRUE, r2 = TRUE, r2adj = TRUE)
  expect_length(n, 2L)
  expect_match(n[1], "^Deviance: D = ")
  expect_match(n[2], "R-squared = .*adjusted R-squared = ")
})

test_that("knot sensitivity takes explicit knot sets", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  ks <- easydr_knot_sensitivity(f, sets = list(c(5, 20, 50), c(10, 25, 60)),
                                plot = FALSE)
  expect_equal(nrow(ks$table), 2L)
  expect_equal(ks$table$Knots, c("5, 20, 50", "10, 25, 60"))
  expect_error(easydr_knot_sensitivity(f, sets = list(c(5, 20, 50))),
               "at least two")
})
