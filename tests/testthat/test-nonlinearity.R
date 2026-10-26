# The departure-from-linearity test (FR-007 to FR-010, FR-024, SC-002).
#
# Reference values come from an engine run on dosresmeta 2.2.0 under R 4.6.0 on
# 2026-10-23, using the engine's own alcohol_crc data. They are engine runs, not
# figures read from a publication, and are labelled as such.

test_that("the rcs fit reproduces the captured Wald test", {
  fit <- easydrbin(crc_table(), sm = "IRR", shape = "rcs",
                   knots = c(5, 25, 50))
  expect_true(fit$nonlinearity$estimable)
  expect_equal(fit$nonlinearity$statistic, 3.7784, tolerance = 1e-3)
  expect_equal(fit$nonlinearity$df, 1)
  expect_equal(fit$nonlinearity$p, 0.05192, tolerance = 1e-3)
  expect_match(fit$nonlinearity$comparison, "against the linear fit")
})

test_that("the basis parameterization is what the test depends on", {
  # With the dose in the first basis column the test asks whether the nonlinear
  # part is needed, and here the answer is borderline: p = 0.052. A basis that
  # mixes the dose into every column gives p = 3.8e-06 on the same data, which
  # is a different question wearing the same label (research decision D1). If a
  # future change swaps the basis to drop a dependency, this test fails.
  fit <- easydrbin(crc_table(), sm = "IRR", shape = "rcs",
                   knots = c(5, 25, 50))
  expect_gt(fit$nonlinearity$p, 0.01)
  expect_lt(fit$nonlinearity$p, 0.10)
  B <- fit$shape$basis(crc_table()$dose)
  expect_equal(as.numeric(B[, 1]), crc_table()$dose)
})

test_that("the quadratic fit reproduces its captured Wald test", {
  fit <- easydrbin(crc_table(), sm = "IRR", shape = "quadratic")
  expect_equal(fit$nonlinearity$statistic, 3.9929, tolerance = 1e-3)
  expect_equal(fit$nonlinearity$p, 0.04569, tolerance = 1e-3)
})

test_that("every multi-term shape reports a test, and one-term shapes do not", {
  # SC-002.
  for (s in c("rcs", "fp", "quadratic")) {
    fit <- easydrbin(crc_table(), sm = "IRR", shape = s)
    expect_true(fit$nonlinearity$estimable, info = s)
    out <- paste(capture.output(print(fit)), collapse = " ")
    expect_match(out, "Departure from linearity: X2 =", info = s)
    expect_false(grepl("nothing to test it against", out))
  }
  for (s in c("linear", "log")) {
    fit <- easydrbin(crc_table(), sm = "IRR", shape = s)
    expect_false(fit$nonlinearity$estimable, info = s)
    out <- paste(capture.output(print(fit)), collapse = " ")
    expect_match(out, "nothing to test it against", info = s)
  }
})

test_that("goodness of fit comes from the engine and is reported", {
  fit <- easydrbin(crc_table(), sm = "IRR", shape = "rcs",
                   knots = c(5, 25, 50))
  expect_true(fit$gof$available)
  expect_equal(fit$gof$deviance, 44.0347, tolerance = 1e-3)
  expect_equal(fit$gof$df, 38)
  expect_equal(fit$gof$p, 0.2314, tolerance = 1e-3)
  expect_equal(fit$gof$r2, 0.3787, tolerance = 1e-3)
  expect_match(paste(capture.output(print(fit)), collapse = " "),
               "Goodness of fit: deviance D =")
})

test_that("the test never changes the model that is reported", {
  # FR-010: a non-significant test is information, not a trigger to refit.
  fit <- easydrbin(crc_table(), sm = "IRR", shape = "rcs",
                   knots = c(5, 25, 50))
  expect_gt(fit$nonlinearity$p, 0.05)
  expect_equal(fit$shape$name, "rcs")
  expect_null(fit$slope)
})
