# The printed answer. One print must carry the whole result, labelled so no
# number can be misread (FR-011, FR-013, FR-014, SC-002, SC-004, SC-007).

out_of <- function(x) paste(capture.output(print(x)), collapse = "\n")

test_that("one print carries every block of the contract layout", {
  fit <- easydrbin(crc_table(), sm = "IRR", dose.unit = "g/day",
                   at = c(10, 25, 50))
  out <- out_of(fit)

  # SC-002: no follow-up call needed.
  expect_match(out, "Dose-response meta-analysis")
  expect_match(out, "incidence rate ratio")
  expect_match(out, "Pooling: two-stage random effects")
  expect_match(out, "Covariance: Greenland and Longnecker")
  expect_match(out, "Per 1 g/day increase")
  expect_match(out, "Reference dose: 0 g/day")
  expect_match(out, "Departure from linearity")
  expect_match(out, "Heterogeneity: Q =")
  expect_match(out, "I-squared")
  expect_match(out, "Predicted IRR")
  expect_match(out, "95% CI")
  expect_match(out, "95% PI")
  expect_match(out, "8 studies, 48 dose categories")
  expect_match(out, "Defaults applied")
})

test_that("the slope is labelled by easydr and never as an intercept", {
  # The engine prints the dose slope as "(Intercept)" because the model is
  # fitted without one. A reviewer copying that would report a slope as an
  # intercept, so easydr builds its own label (research D6).
  fit <- easydrbin(crc_table(), sm = "IRR", dose.unit = "g/day")
  out <- out_of(fit)
  expect_match(out, "Per 1 g/day increase")
  expect_false(grepl("Intercept", out, fixed = TRUE))
})

test_that("every number carries its scale and its reference dose", {
  # SC-004: a log-scale slope must not read as a ratio.
  fit <- easydrbin(crc_table(), sm = "IRR", dose.unit = "g/day")
  out <- out_of(fit)
  expect_match(out, "IRR")
  expect_match(out, "Reference dose")
  # The reported slope is back-transformed, so it exceeds 1 for a positive log
  # slope rather than printing the log value.
  expect_gt(exp(fit$slope$estimate), 1)
  expect_match(out, "IRR 1[.]006")
})

test_that("the dose unit falls back to a readable phrase", {
  fit <- easydrbin(crc_table(), sm = "IRR")
  expect_match(out_of(fit), "Per 1 dose unit increase")
})

test_that("prediction and confidence intervals are separate labelled columns", {
  fit <- easydrbin(crc_table(), sm = "IRR", at = c(10, 25))
  out <- out_of(fit)
  expect_match(out, "95% CI")
  expect_match(out, "95% PI")
  expect_false(identical(fit$predicted$ci.lb, fit$predicted$pi.lb))
})

test_that("the prediction interval is wider than the confidence interval", {
  # alcohol_cvd has Psi = 1.0057e-4 and Q = 14.146 on 5 df, so there is real
  # between-study variation for the prediction interval to describe.
  fit <- easydrbin(cvd_table(), sm = "RR", dose.unit = "g/day",
                   at = c(10, 25, 50))
  pred <- fit$predicted[!fit$predicted$is_ref, ]
  expect_true(fit$prediction$estimable)
  expect_gt(as.numeric(fit$Psi), 0)
  expect_true(all(pred$pi.lb < pred$ci.lb))
  expect_true(all(pred$pi.ub > pred$ci.ub))
  expect_gt(fit$heterogeneity$i2, 0)
})

test_that("a prediction interval with no between-study variance matches the CI", {
  # alcohol_crc has Psi at the numerical floor, so the two intervals coincide
  # apart from the t against z critical value. The interval is still printed as
  # a prediction interval, never dropped.
  fit <- easydrbin(crc_table(), sm = "IRR", at = c(10, 50))
  pred <- fit$predicted[!fit$predicted$is_ref, ]
  expect_true(fit$prediction$estimable)
  expect_true(all(pred$pi.lb <= pred$ci.lb))
  expect_true(all(pred$pi.ub >= pred$ci.ub))
})

test_that("defaults applied are printed, so a methods paragraph can be written", {
  # SC-007.
  fit <- easydrbin(crc_table(), sm = "IRR")
  out <- out_of(fit)
  expect_match(out, "intake: reported effect with limits")
  expect_match(out, "reference category inferred in 8 studies")
  expect_match(out, "two-stage random effects, reml")
})

test_that("the continuous layout matches the binary layout block for block", {
  # User Story 2 acceptance 2: same shape, only the measure and covariance
  # lines differ.
  blocks <- function(x) {
    out <- out_of(x)
    c(grepl("Dose-response meta-analysis", out),
      grepl("Pooling:", out),
      grepl("increase", out),
      grepl("Reference dose:", out),
      grepl("Departure from linearity", out),
      grepl("Heterogeneity:", out),
      grepl("Predicted", out),
      grepl("dose categories", out),
      grepl("Defaults applied", out))
  }
  expect_equal(blocks(easydrbin(crc_table(), sm = "IRR")),
               blocks(easydrcont(ari_table(), sm = "SMD")))
})

test_that("summary adds the per-study contribution table", {
  s <- summary(easydrbin(crc_table(), sm = "IRR", dose.unit = "g/day"))
  expect_s3_class(s, "summary.easydr_result")
  expect_equal(nrow(s$per_study), 8)
  expect_true(all(c("study", "levels", "dose_range", "design", "slope", "se")
                  %in% names(s$per_study)))
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "Per-study contribution")
  expect_match(out, "slope [(]log scale[)]")
})

test_that("print returns its input invisibly", {
  fit <- easydrbin(crc_table(), sm = "IRR")
  expect_invisible(print(fit))
})
