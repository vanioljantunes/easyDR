# Meta-regression. The reference values in helper-anchors.R are captured engine
# runs on the bundled alcohol_cvd table, not published figures.
#
# T152 looked for a published anchor for this path, as Principle VI asks, and
# found none. The dosresmeta help page for gof() documents a meta-regression on
# milk_ov by study type, but Rd examples carry no expected values; the JSS paper
# (Crippa and Orsini 2016) cites meta-regression methods without printing a
# worked dose-response meta-regression; and alecri.github.io links its analyses
# out to rpubs rather than showing numbers. So this path is anchored to a
# captured run, and the gap is stated here rather than left implicit. If a
# publication reporting a dose-response meta-regression with its coefficients
# turns up, its figures should replace these.

test_that("a moderator reproduces the captured coefficients and tests", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  m   <- easyDR:::easydr_metareg(fit, ~ design)

  expect_s3_class(m, "easydr_moderator")
  expect_equal(nrow(m$terms), 2)
  expect_equal(m$terms$estimate[1], cvd_mod_slope_ref, tolerance = 1e-5)
  expect_equal(m$terms$estimate[2], cvd_mod_slope_change, tolerance = 1e-5)

  expect_true(m$test$estimable)
  expect_equal(m$test$chi2, cvd_mod_wald_chi2, tolerance = 1e-5)
  expect_identical(m$test$df, cvd_mod_wald_df)
  expect_equal(m$test$p, cvd_mod_wald_p, tolerance = 1e-5)
})

test_that("the residual heterogeneity is reported beside the heterogeneity without", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  m   <- easyDR:::easydr_metareg(fit, ~ design)

  expect_equal(m$heterogeneity$Q, cvd_mod_Q, tolerance = 1e-5)
  expect_identical(as.integer(m$heterogeneity$df), cvd_mod_df)
  expect_equal(m$heterogeneity$p, cvd_mod_p, tolerance = 1e-4)

  # The comparison is the fit the moderator was added to, not a refit.
  expect_equal(m$heterogeneity_before$Q, cvd_full_Q, tolerance = 1e-6)
  expect_identical(as.integer(m$heterogeneity_before$df), cvd_full_df)
})

test_that("easydr labels the coefficients and never leaves them as the engine does", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  m   <- easyDR:::easydr_metareg(fit, ~ design)

  # What the engine calls them. Both would be misread: the first looks like a
  # baseline risk and the second like a risk difference. Both are slopes.
  expect_identical(m$terms$coefficient,
                   c("(Intercept)", "designcumulative incidence"))
  expect_match(m$terms$label[1], "^slope at design = case-control$")
  expect_match(m$terms$label[2], "^change in slope, cumulative incidence$")

  out <- paste(capture.output(print(m)), collapse = "\n")
  expect_false(grepl("(Intercept)", out, fixed = TRUE))
  expect_false(grepl("designcumulative incidence", out, fixed = TRUE))
})

test_that("the reference level is named, because a coefficient without it cannot be read", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  m   <- easyDR:::easydr_metareg(fit, ~ design)

  expect_identical(m$reference_levels$design, cvd_mod_reference)
  out <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(out, "reference level case-control")
})

test_that("the printed output says what a coefficient is and that it is observational", {
  # Both sentences are required by the spec. Feature 002 lost a print
  # guarantee to a later change and only a pinned assertion caught it, so
  # these get a test rather than the memory of a reviewer.
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  out <- paste(capture.output(print(easyDR:::easydr_metareg(fit, ~ design))),
               collapse = "\n")

  expect_match(out,
    "A moderator coefficient is an effect on the dose-response, not on the\noutcome")
  expect_match(out,
    "observational association across studies, not a randomised\ncomparison")
  expect_match(out, "has not been shown to\ncause it")
  expect_match(out, "Wald test on the 1 moderator column")
  expect_match(out, "without the moderator: Q =")
})

test_that("a multi-term shape gets one moderator coefficient per shape term", {
  d <- alcohol_cancer
  d$year <- c(1990, 1995, 1998, 2001, 2004, 2007, 2010, 2013)[
    match(d$studlab, unique(d$studlab))]
  m <- easyDR:::easydr_metareg(easydrbin(d, sm = "IRR", shape = "quadratic"), ~ year)

  expect_equal(nrow(m$terms), 4)
  expect_equal(sum(m$terms$moderator), 2)
  expect_identical(m$test$df, 2L)
  expect_match(m$terms$label[1], "the dose term at the mean of year")
  expect_match(m$terms$label[3], "change in the dose term per unit of year")
  expect_false(any(grepl("(Intercept)", m$terms$label, fixed = TRUE)))
})

test_that("a continuous moderator is centred, which moves only the reference coefficient", {
  # Shifting a covariate cannot change its slope, so the moderator effect and
  # its test must be identical centred or not. What changes is the reference
  # coefficient, from the dose-response at year zero, an extrapolation two
  # thousand years out, to the dose-response at the mean year.
  d <- alcohol_cancer
  yrs <- c(1990, 1995, 1998, 2001, 2004, 2007, 2010, 2013)
  d$year <- yrs[match(d$studlab, unique(d$studlab))]
  dc <- d
  dc$year <- dc$year - mean(yrs)

  m  <- easyDR:::easydr_metareg(easydrbin(d,  sm = "IRR"), ~ year)
  mc <- easyDR:::easydr_metareg(easydrbin(dc, sm = "IRR"), ~ year)

  expect_equal(m$terms$estimate[2], mc$terms$estimate[2], tolerance = 1e-10)
  expect_equal(m$test$chi2, mc$test$chi2, tolerance = 1e-10)
  expect_equal(m$centres$year, mean(yrs))
  expect_equal(m$terms$estimate[1], mc$terms$estimate[1], tolerance = 1e-10)

  out <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(out, "was centred at its mean across studies")
  expect_match(out, "The moderator effect and its test are unchanged")
})

test_that("V19 stops a covariate that varies inside a study", {
  bad <- alcohol_cvd
  bad$site <- seq_len(nrow(bad))
  fit <- easydrbin(alcohol_cvd, sm = "RR")
  err <- tryCatch(easyDR:::easydr_metareg(fit, ~ site, data = bad),
                  error = conditionMessage)

  expect_match(err, "takes 5 different values inside study Bianchi")
  expect_match(err, "one value per study")
})

test_that("V20 stops a covariate that is the same in every study", {
  con <- alcohol_cvd
  con$one <- "same"
  fit <- easydrbin(alcohol_cvd, sm = "RR")
  err <- tryCatch(easyDR:::easydr_metareg(fit, ~ one, data = con),
                  error = conditionMessage)

  expect_match(err, "takes the same value")
  expect_match(err, "not estimable")
})

test_that("V21 stops when the model needs as many coefficients as there are studies", {
  four <- alcohol_cvd[alcohol_cvd$studlab %in%
    c("Bianchi", "Bobak", "Vliegenthart-wine", "Vliegenthart-beer"), ]
  fit <- easydrbin(four, sm = "RR", shape = "quadratic")
  err <- tryCatch(easyDR:::easydr_metareg(fit, ~ design), error = conditionMessage)

  expect_match(err, "needs 4 coefficients")
  expect_match(err, "only 4 studies")
})

test_that("V23 flags a continuous moderator whose range rests on one study", {
  d <- alcohol_cancer
  d$year <- c(1990, 1995, 1998, 2001, 2004, 2007, 2010, 2013)[
    match(d$studlab, unique(d$studlab))]
  d$year[d$studlab == unique(d$studlab)[1]] <- 1890
  fit <- easydrbin(d, sm = "IRR")

  expect_warning(m <- easyDR:::easydr_metareg(fit, ~ year), "rests on one study")
  expect_identical(m$narrow, "year")
  out <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(out, "rests on one study")
})

test_that("a moderator must be a one-sided formula naming columns that exist", {
  fit <- easydrbin(alcohol_cvd, sm = "RR")
  expect_error(easyDR:::easydr_metareg(fit, ~ nope), "not in the data")
  expect_error(easyDR:::easydr_metareg(fit, y ~ design), "must be one-sided")
  expect_error(easyDR:::easydr_metareg(fit, "design"), "one-sided formula")
})

test_that("the moderator fit inherits the model it was handed", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day",
                   proc = "1stage", method.tau = "ML")
  m   <- easyDR:::easydr_metareg(fit, ~ design)

  expect_identical(m$inherited$proc, "1stage")
  expect_identical(m$inherited$method, fit$provenance$method)
  expect_identical(m$inherited$sm, fit$sm)
  out <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(out, "one-stage")
})

test_that("two moderators at once are labelled and tested together", {
  d <- alcohol_cancer
  u <- unique(d$studlab)
  d$region <- ifelse(d$studlab %in% u[1:4], "europe", "america")
  d$year   <- c(1990, 1995, 1998, 2001, 2004, 2007, 2010, 2013)[
    match(d$studlab, u)]
  m <- easyDR:::easydr_metareg(easydrbin(d, sm = "IRR"), ~ region + year)

  expect_equal(nrow(m$terms), 3)
  expect_equal(sum(m$terms$moderator), 2)
  expect_identical(m$test$df, 2L)
  expect_identical(m$reference_levels$region, "america")
  # The reference coefficient sits at a level for one and at a mean for the
  # other, and the label has to say both.
  expect_match(m$terms$label[1], "region = america")
  expect_match(m$terms$label[1], "the mean of year")
})
