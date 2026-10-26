# Validation against published results, binary outcomes (FR-025, FR-026).
#
# Estimation paths reachable from easydrbin, and the check covering each:
#   reported effect + counts, Greenland-Longnecker, two-stage -> first block
#   reported effect + counts, Greenland-Longnecker, one-stage -> one-stage block
#   reported effect + counts, Hamling                         -> Hamling block
#   counts only, effect computed through escalc                -> counts block
#   reported effect only, no denominators, independent         -> indep block
#
# Source for the first block: alcohol and colorectal cancer, the eight-cohort
# example distributed with dosresmeta 2.2.0 and analysed in Orsini et al.
# (2012) AJE 175:66-73. Values captured from a direct engine run on R 4.6.0 on
# 2026-10-23, recorded in specs/001-core-dr-fit/research.md decision D8.

test_that("alcohol and colorectal cancer reproduces the published slope", {
  fit <- easydrbin(crc_table("effect_counts"), sm = "IRR",
                   dose.unit = "g/day", at = c(10, 25, 50))

  # Slope on the log scale, per 1 g/day. Tolerance 1e-4, as D8 states.
  expect_equal(fit$slope$estimate, 0.006437538, tolerance = 1e-4)
  expect_equal(fit$slope$se, 0.001352, tolerance = 1e-3)

  # Back-transformed rate ratios. Tolerance 1e-3, as D8 states.
  pred <- fit$predicted
  expect_equal(pred$est[pred$dose == 10], 1.066493, tolerance = 1e-3)
  expect_equal(pred$est[pred$dose == 25], 1.174613, tolerance = 1e-3)
  expect_equal(pred$est[pred$dose == 50], 1.379715, tolerance = 1e-3)

  expect_equal(pred$ci.lb[pred$dose == 10], 1.038602, tolerance = 1e-3)
  expect_equal(pred$ci.ub[pred$dose == 10], 1.095132, tolerance = 1e-3)
  expect_equal(pred$ci.lb[pred$dose == 50], 1.208499, tolerance = 1e-3)
  expect_equal(pred$ci.ub[pred$dose == 50], 1.575188, tolerance = 1e-3)

  # Heterogeneity is the engine's own, not a competing statistic.
  expect_true(fit$heterogeneity$estimable)
  expect_equal(fit$heterogeneity$Q, 4.779728, tolerance = 1e-4)
  expect_equal(fit$heterogeneity$df, 7)
  expect_equal(fit$heterogeneity$p, 0.6868233, tolerance = 1e-4)
  expect_equal(fit$heterogeneity$i2, 0, tolerance = 1e-6)

  expect_equal(fit$counts$studies, 8)
  expect_equal(fit$provenance$covariance, "gl")
})

test_that("one-stage pooling reaches the engine", {
  one <- easydrbin(crc_table("effect_counts"), sm = "IRR", proc = "1stage")
  two <- easydrbin(crc_table("effect_counts"), sm = "IRR", proc = "2stage")

  expect_equal(one$provenance$proc, "1stage")
  expect_match(paste(capture.output(print(one)), collapse = " "), "one-stage")

  # On this dataset the between-study variance is at the numerical floor, so
  # the two procedures agree to within 1e-8. That is a fact about the data, not
  # a defect: the check is that the argument reaches the engine and that the
  # estimate stays sane, not that the two differ.
  expect_equal(one$slope$estimate, two$slope$estimate, tolerance = 1e-6)
  expect_gt(one$slope$estimate, 0)
})

test_that("Hamling reconstruction is available and agrees closely with GL", {
  # Engine run, R 4.6.0, 2026-10-23: covariance = "h" on the same data gives
  # 0.00640759662759221, against 0.006437538 for Greenland-Longnecker.
  ham <- easydrbin(crc_table("effect_counts"), sm = "IRR", covariance = "h")
  expect_equal(ham$provenance$covariance, "h")
  expect_equal(ham$slope$estimate, 0.0064076, tolerance = 1e-4)
})

test_that("counts-only intake computes the effect and says so", {
  # No published anchor exists for a crude effect recomputed from these counts,
  # so this check fixes the behaviour rather than a published number: the
  # intake is recorded, the fit runs, and the slope keeps its direction.
  fit <- easydrbin(crc_table("counts"), sm = "IRR", dose.unit = "g/day")
  expect_equal(fit$provenance$intake, "counts")
  expect_equal(fit$provenance$covariance, "gl")
  expect_gt(fit$slope$estimate, 0)
})

test_that("effect without denominators falls back to independent contrasts", {
  fit <- easydrbin(crc_table("effect"), sm = "IRR", dose.unit = "g/day")
  expect_equal(fit$provenance$intake, "effect")
  expect_equal(fit$provenance$covariance, "indep")
  # The printed output must disclose the fallback, since it understates the
  # correlation between a study's categories.
  expect_match(paste(capture.output(print(fit)), collapse = " "),
               "independent contrasts")
})

test_that("the two intakes of the same data agree", {
  gl    <- easydrbin(crc_table("effect_counts"), sm = "IRR")
  indep <- easydrbin(crc_table("effect"), sm = "IRR")
  expect_equal(gl$slope$estimate, indep$slope$estimate, tolerance = 0.2)
})
