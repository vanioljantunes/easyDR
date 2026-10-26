# Shape selection (FR-016 to FR-023, SC-004, SC-005).
#
# A small candidate set keeps these tests quick. The full set of about forty is
# exercised once, in the quickstart, not on every run.

small <- list(
  list(shape = "linear"),
  list(shape = "quadratic"),
  list(shape = "fp", powers = c(-1, 0)),
  list(shape = "rcs", knots = 3)
)

test_that("every candidate comes back, successes and failures alike", {
  # SC-004: the row count equals the candidate count.
  sel <- easydr_select(crc_table(), sm = "IRR", candidates = small, quiet = TRUE)
  expect_s3_class(sel, "easydr_selection")
  expect_equal(nrow(sel$ranked), length(small))
  expect_equal(sel$tried, length(small))
  expect_equal(sel$fitted + sel$failed, length(small))
  expect_true(all(c("shape", "terms", "AIC", "BIC", "deviance_p", "status",
                    "reason", "rank") %in% names(sel$ranked)))
})

test_that("the straight line is among the candidates", {
  # FR-017: without it the ranking cannot answer whether a curve is warranted.
  sel <- easydr_select(crc_table(), sm = "IRR", candidates = small, quiet = TRUE)
  expect_true("linear" %in% sel$ranked$shape)
  full <- dr_candidates()
  expect_true(any(vapply(full, function(s) identical(s$shape, "linear"),
                         logical(1))))
  # The documented set is bounded and reportable: 3 named shapes, 36 power
  # pairs, 3 spline knot counts.
  expect_equal(length(full), 42)
})

test_that("a failing candidate is recorded and does not stop the loop", {
  # FR-019, rule V18. Knots outside the dose range are refused by V16, so this
  # candidate must fail while its neighbours still fit.
  cands <- c(small, list(list(shape = "rcs", knots = c(500, 600, 700))))
  sel <- easydr_select(crc_table(), sm = "IRR", candidates = cands, quiet = TRUE)
  expect_equal(nrow(sel$ranked), length(cands))
  expect_equal(sel$failed, 1L)
  failed <- sel$ranked[sel$ranked$status == "failed", ]
  expect_match(failed$reason, "outside the observed dose range")
  expect_true(is.na(failed$rank))
  expect_gt(sel$fitted, 0L)
})

test_that("the ranking statistic is recorded and respected", {
  # FR-018: never implicit.
  on_aic <- easydr_select(crc_table(), sm = "IRR", candidates = small,
                      rank_on = "AIC", quiet = TRUE)
  on_bic <- easydr_select(crc_table(), sm = "IRR", candidates = small,
                      rank_on = "BIC", quiet = TRUE)
  expect_equal(on_aic$ranked_on, "AIC")
  expect_equal(on_bic$ranked_on, "BIC")
  fitted <- on_aic$ranked[on_aic$ranked$status == "fitted", ]
  expect_equal(fitted$AIC, sort(fitted$AIC))
  expect_match(paste(capture.output(print(on_bic)), collapse = " "),
               "Ranked on BIC")
})

test_that("candidates are refitted one-stage by maximum likelihood", {
  # A two-stage information criterion moves with the scale of the dose terms:
  # multiplying the dose by 1000 changes no model but moves the two-stage AIC
  # from -64.5 to -175.1. One-stage keeps the same outcome for every candidate.
  sel <- easydr_select(crc_table(), sm = "IRR", candidates = small, quiet = TRUE)
  expect_true(sel$forced_ml)
  expect_true(sel$forced_1stage)
  expect_equal(sel$best$provenance$proc, "1stage")
  expect_equal(sel$best$provenance$method, "ml")
  out <- paste(capture.output(print(sel)), collapse = " ")
  expect_match(out, "Refitted one-stage by maximum likelihood")
  expect_match(out, "moves with the scale")
})

test_that("asking for two-stage ranking warns rather than quietly misleading", {
  expect_warning(
    easydr_select(crc_table(), sm = "IRR", candidates = small, proc = "2stage",
              quiet = TRUE),
    "rewards a change of units")
})

test_that("ties are shown rather than broken silently", {
  # FR-022. On this data several candidates sit within 2 AIC of the best, and
  # the straight line is among them, which is the finding a reviewer needs.
  sel <- easydr_select(crc_table(), sm = "IRR", candidates = small, quiet = TRUE)
  expect_true(length(sel$ties) >= 2L)
  expect_match(paste(capture.output(print(sel)), collapse = " "),
               "of the best:")
})

test_that("the best-ranked fit comes back ready to use", {
  # FR-020.
  sel <- easydr_select(crc_table(), sm = "IRR", candidates = small, quiet = TRUE)
  expect_s3_class(sel$best, "easydr_result")
  expect_match(sel$ranked$shape[1], sel$best$shape$name)
})

test_that("the printed output calls the ranking a proposal and says how to adopt it", {
  # SC-005: a reviewer must be able to tell this from the output alone.
  sel <- easydr_select(crc_table(), sm = "IRR", candidates = small, quiet = TRUE)
  out <- paste(capture.output(print(sel)), collapse = " ")
  expect_match(out, "a proposal, not a decision")
  expect_match(out, "state it yourself")
  expect_match(out, "easydr[(]data")
  expect_match(out, "candidates tried")
  expect_match(out, "seconds")
})

test_that("selection never touches a fit the reviewer made themselves", {
  # FR-021. The guarantee is structural: easydr_select holds no reference to any
  # other fit. This test pins the behaviour anyway.
  stated <- easydrbin(crc_table(), sm = "IRR", shape = "rcs",
                      knots = c(5, 25, 50))
  before <- capture.output(print(stated))
  invisible(easydr_select(crc_table(), sm = "IRR", candidates = small,
                      quiet = TRUE))
  expect_equal(capture.output(print(stated)), before)
  expect_equal(stated$shape$knots, c(5, 25, 50))
  expect_equal(stated$provenance$proc, "2stage")
})

test_that("the shape arguments cannot be passed to the loop that varies them", {
  expect_error(easydr_select(crc_table(), sm = "IRR", shape = "rcs", quiet = TRUE),
               "varies the shape itself")
  expect_error(easydr_select(crc_table(), sm = "IRR", powers = c(1, 2),
                         quiet = TRUE),
               "varies the shape itself")
})

test_that("selection picks the continuous worker from sm, with no fn argument", {
  # The caller never names a fitting function: sm decides, exactly as it does
  # in easydr().
  sel <- easydr_select(ari_table(), sm = "SMD", candidates = small[1:3],
                       quiet = TRUE)
  expect_equal(nrow(sel$ranked), 3)
  expect_gt(sel$fitted, 0L)
  expect_match(paste(capture.output(print(sel)), collapse = " "),
               "easydr[(]data")
})
