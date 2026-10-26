# Validation against published results, continuous outcomes (FR-025).
#
# Estimation paths reachable from easydrcont, and the check covering each:
#   standardised mean difference, two-stage -> SMD block
#   mean difference, two-stage              -> MD block
#   standardised mean difference, one-stage -> one-stage block
#
# Source: the aripiprazole example distributed with dosresmeta 2.2.0 and
# analysed in Crippa and Orsini (2016) BMC Med Res Methodol 16:91. Expected
# values captured from a direct engine run on R 4.6.0 on 2026-10-23:
#   covariance = "smd": slope 0.01532420243, SE 0.002596301563
#   covariance = "md" : slope 0.3632491684,  SE 0.06035296574
#   qtest: Q = 2.851, df = 4, p = 0.583

test_that("aripiprazole reproduces the published standardised slope", {
  fit <- easydrcont(ari_table(), sm = "SMD", dose.unit = "mg/day",
                    at = c(5, 10, 20))

  expect_equal(fit$slope$estimate, 0.01532420243, tolerance = 1e-6)
  expect_equal(fit$slope$se, 0.002596301563, tolerance = 1e-6)

  expect_true(fit$heterogeneity$estimable)
  expect_equal(fit$heterogeneity$Q, 2.851, tolerance = 1e-3)
  expect_equal(fit$heterogeneity$df, 4)
  expect_equal(fit$heterogeneity$p, 0.583, tolerance = 1e-2)

  expect_equal(fit$provenance$covariance, "smd")
  expect_equal(fit$counts$studies, 5)

  # A standardised mean difference is not back-transformed: the predicted value
  # is the slope times the dose contrast, on the scale it was estimated on.
  pred <- fit$predicted
  expect_equal(pred$est[pred$dose == 10], 10 * 0.01532420243, tolerance = 1e-6)
})

test_that("aripiprazole reproduces the published mean-difference slope", {
  fit <- easydrcont(ari_table(), sm = "MD", dose.unit = "mg/day")
  expect_equal(fit$slope$estimate, 0.3632491684, tolerance = 1e-6)
  expect_equal(fit$slope$se, 0.06035296574, tolerance = 1e-6)
  expect_equal(fit$provenance$covariance, "md")
})

test_that("one-stage pooling works for continuous outcomes too", {
  fit <- easydrcont(ari_table(), sm = "SMD", proc = "1stage")
  expect_equal(fit$provenance$proc, "1stage")
  expect_gt(fit$slope$estimate, 0)
})

test_that("the standardisation method is honoured", {
  cohens <- easydrcont(ari_table(), sm = "SMD", method.smd = "cohens")
  hedges <- easydrcont(ari_table(), sm = "SMD", method.smd = "hedges")
  expect_false(isTRUE(all.equal(cohens$slope$estimate, hedges$slope$estimate,
                                tolerance = 1e-10)))
})

# Feature 002: nonlinear shapes on the continuous path.
#
# Captured from easydr on dosresmeta 2.2.0 under R 4.6.0 on 2026-10-23, using
# the engine's own ari data. These are engine runs, not figures read from a
# publication, and are labelled as such. Where a paper reports the same analysis
# with the same knots, its figure should replace these.

test_that("the quadratic shape on ari reproduces its captured values", {
  fit <- easydrcont(ari_table(), sm = "SMD", shape = "quadratic")
  expect_equal(as.numeric(coef(fit$fit)),
               c(0.046070605, -0.0011202281), tolerance = 1e-6)
  expect_true(fit$nonlinearity$estimable)
  expect_equal(fit$nonlinearity$statistic, 14.795224, tolerance = 1e-5)
  expect_equal(fit$nonlinearity$df, 1)
  # expect_equal compares on a relative scale, so a p value of order 1e-4 needs
  # its captured digits rather than a rounded display value.
  expect_equal(fit$nonlinearity$p, 0.0001198386, tolerance = 1e-5)
  expect_equal(fit$gof$deviance, 6.691022, tolerance = 1e-5)
  expect_equal(fit$gof$p, 0.8235248, tolerance = 1e-5)
  expect_equal(fit$gof$df, 11)
  expect_null(fit$slope)
})

test_that("the spline shape on ari reproduces its captured values", {
  fit <- easydrcont(ari_table(), sm = "SMD", shape = "rcs")
  expect_equal(as.numeric(coef(fit$fit)),
               c(0.035448972, -0.035913366), tolerance = 1e-6)
  expect_equal(fit$nonlinearity$statistic, 14.069064, tolerance = 1e-5)
  expect_equal(fit$nonlinearity$p, 0.0001762185, tolerance = 1e-5)
  expect_equal(fit$gof$deviance, 7.417995, tolerance = 1e-5)
  expect_equal(fit$gof$p, 0.7642941, tolerance = 1e-5)
  # Unlike the alcohol data, here the departure from linearity is clear.
  expect_lt(fit$nonlinearity$p, 0.001)
})

test_that("the fractional polynomial on ari fits and reports its shift", {
  fit <- easydrcont(ari_table(), sm = "SMD", shape = "fp", powers = c(-1, 0))
  expect_equal(fit$shape$powers, c(-1, 0))
  # A dose of zero is present, so the shift must be reported.
  expect_match(fit$shape$label, "shifted by")
  expect_true(fit$gof$available)
})

test_that("a continuous curve prints the same blocks as a binary curve", {
  # FR-026: the layout does not change with the outcome type.
  blocks <- function(x) {
    out <- paste(capture.output(print(x)), collapse = "
")
    c(grepl("Dose-response meta-analysis", out),
      grepl("Pooling:", out),
      grepl("Shape:", out),
      grepl("No single slope describes this shape", out),
      grepl("Reference dose:", out),
      grepl("Departure from linearity: X2 =", out),
      grepl("Goodness of fit: deviance D =", out),
      grepl("Heterogeneity:", out),
      grepl("Predicted", out),
      grepl("dose categories", out),
      grepl("Defaults applied", out))
  }
  expect_equal(blocks(easydrcont(ari_table(), sm = "SMD", shape = "rcs")),
               blocks(easydrbin(crc_table(), sm = "IRR", shape = "rcs")))
  expect_true(all(blocks(easydrcont(ari_table(), sm = "SMD",
                                    shape = "quadratic"))))
})

test_that("the summary view carries the goodness of fit", {
  s <- summary(easydrcont(ari_table(), sm = "SMD", shape = "rcs"))
  out <- paste(capture.output(print(s)), collapse = "
")
  expect_match(out, "Goodness of fit: deviance D =")
  expect_match(out, "R-squared")
  expect_match(out, "Per-study contribution")
})
