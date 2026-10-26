# easydr() is the one public entry point. These tests pin the three things a
# reader is entitled to rely on: that it reproduces the published analysis, that
# it refuses to choose a shape, and that it folds the two checks in rather than
# leaving them to be remembered.

test_that("easydr reproduces the published alcohol and colorectal cancer fit", {
  # Orsini et al. 2012, eight prospective cohorts. The anchors are the pooled
  # slope on the log scale and the rate ratios at intakes people recognise.
  fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
                dose.unit = "g/day", at = c(10, 25, 50))

  expect_s3_class(fit, "easydr_result")
  expect_equal(fit$slope$estimate, 0.006437300, tolerance = 1e-7)
  expect_equal(fit$slope$se, 0.001352076, tolerance = 1e-7)

  at <- fit$predicted[!fit$predicted$is_ref, ]
  expect_equal(at$dose, c(10, 25, 50))
  expect_equal(round(at$est, 4), c(1.0665, 1.1746, 1.3797))
  expect_equal(round(at$ci.lb, 4), c(1.0386, 1.0993, 1.2085))
  expect_equal(round(at$ci.ub, 4), c(1.0951, 1.2551, 1.5752))

  expect_equal(fit$heterogeneity$Q, 4.78, tolerance = 0.01)
  expect_equal(fit$heterogeneity$df, 7L)
  expect_equal(fit$heterogeneity$i2, 0)
})

test_that("easydr gives the same answer as the worker it delegates to", {
  # The adapter must add the two checks and nothing else. If these ever differ,
  # easydr() has started doing arithmetic of its own.
  a <- easydr(alcohol_cancer, sm = "IRR", shape = "linear", dose.unit = "g/day")
  b <- easydrbin(alcohol_cancer, sm = "IRR", shape = "linear",
                 dose.unit = "g/day")
  expect_equal(a$slope$estimate, b$slope$estimate)
  expect_equal(a$slope$se, b$slope$se)
  expect_equal(a$predicted, b$predicted)
  expect_equal(a$curve, b$curve)
  expect_equal(a$heterogeneity, b$heterogeneity)
})

test_that("the shape is never chosen for the caller", {
  expect_error(easydr(alcohol_cancer, sm = "IRR"),
               "name the `shape` you want fitted")
  expect_error(easydr(alcohol_cancer, sm = "IRR"), "There is no default")
  # The error lists what you may name, so the fix is in the message.
  for (shape in c("linear", "rcs", "fp", "quadratic", "log")) {
    expect_error(easydr(alcohol_cancer, sm = "IRR"), shape)
  }
})

test_that("sm selects the outcome kind, and the wrong columns are refused", {
  cont <- easydr(aripiprazole, sm = "SMD", shape = "linear",
                 dose.unit = "mg/day")
  expect_s3_class(cont, "easydr_result")
  expect_equal(cont$entry, "easydrcont")
  expect_false(cont$ratio)

  bin <- easydr(alcohol_cancer, sm = "IRR", shape = "linear")
  expect_equal(bin$entry, "easydrbin")
  expect_true(bin$ratio)

  # A column belonging to the other kind of outcome is an error, not a name
  # that is quietly ignored while the fit looks like it used it.
  expect_error(easydr(alcohol_cancer, sm = "IRR", shape = "linear",
                      sd.e = "sd"),
               "sd.e is not read")
  expect_error(easydr(aripiprazole, sm = "SMD", shape = "linear",
                      TE = "TE"),
               "TE is not read")
})

test_that("every shape reaches easydr through the same argument", {
  for (shape in c("linear", "quadratic", "log", "rcs", "fp")) {
    fit <- easydr(alcohol_cancer, sm = "IRR", shape = shape,
                  dose.unit = "g/day")
    expect_s3_class(fit, "easydr_result")
    expect_equal(fit$shape$name, shape)
  }
})

test_that("the viability check travels on the result and is printed", {
  fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear")
  expect_s3_class(fit$viability, "dr_viability")
  expect_true(fit$viability$ok)
  expect_equal(nrow(fit$viability$summary), 8L)
  expect_true(all(fit$viability$summary$two_stage_ok))
  # The reference comes from the fit, not from a guess at a zero dose.
  expect_true(all(fit$viability$summary$has_reference))

  out <- paste(capture.output(print(fit)), collapse = " ")
  expect_match(out, "Viability: 8 of 8 studies carry 3 or more dose levels")
})

test_that("min.levels moves the floor the check applies and the line reports it", {
  fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear", min.levels = 7L)
  expect_false(fit$viability$ok)
  out <- paste(capture.output(print(fit)), collapse = " ")
  expect_match(out, "7 or more dose levels")
  expect_match(out, "See \\$viability")
})

test_that("the curve is named, and named on the model scale", {
  fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
                dose.unit = "g/day")
  expect_equal(fit$curve_label, "monotonic increasing")
  expect_match(paste(capture.output(print(fit)), collapse = " "),
               "Curve shape: monotonic increasing")

  # The label describes how the effect changes, so taking the log of a
  # monotone curve cannot change it.
  expect_equal(curve_shape(fit$curve$dose, log(fit$curve$est)),
               curve_shape(fit$curve$dose, fit$curve$est))
})

test_that("the workers leave the two new slots empty, so refits print as before", {
  # A sensitivity analysis refits through easydrbin(). If those refits grew a
  # viability block the printed layout of features 001 to 003 would change.
  w <- easydrbin(alcohol_cancer, sm = "IRR", shape = "linear")
  expect_null(w$viability)
  expect_null(w$curve_label)
  out <- paste(capture.output(print(w)), collapse = " ")
  expect_no_match(out, "Viability:")
  expect_no_match(out, "Curve shape:")
})

test_that("an easydr fit can still be handed to the sensitivity analyses", {
  fit <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
  expect_s3_class(easydr_loo(fit, quiet = TRUE), "easydr_sensitivity")
  expect_s3_class(easydr_metareg(fit, ~ design), "easydr_moderator")
})
