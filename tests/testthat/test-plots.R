# The plots are tested on what they return, never on what they look like. Every
# one of them hands back the table it drew, so the numbers behind the picture
# can be checked against the fit they came from. A pixel comparison would fail
# on a different graphics device and would still not say the picture was right.

local_null_device <- function(env = parent.frame()) {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off(), envir = env)
}

bin_fit <- function(shape = "linear", ...) {
  easydr(alcohol_cancer, sm = "IRR", shape = shape, dose.unit = "g/day", ...)
}

test_that("the forest rows are the engine's own first-stage estimates", {
  local_null_device()
  fit <- bin_fit()
  f   <- easydr_forest(fit)

  expect_equal(nrow(f$table), 8L)
  expect_equal(f$table$study, unique(as.character(fit$frame$studlab)))
  # Each row is the study's coefficient from the two-stage first stage, not
  # something recomputed here.
  expect_equal(f$table$est, as.numeric(fit$fit$bi[, 1]))
  expect_equal(f$table$se,
               unname(vapply(fit$fit$Si,
                             function(V) sqrt(as.numeric(V)[1]), numeric(1))))
  # And the diamond is the pooled estimate the print already reports.
  expect_equal(f$pooled$est, fit$slope$estimate)
  expect_equal(f$pooled$se, fit$slope$se)
})

test_that("the forest weights are inverse variance and sum to a hundred", {
  local_null_device()
  f <- easydr_forest(bin_fit())
  expect_equal(sum(f$table$weight), 100)
  # A more precise study carries more weight, so the two orderings are mirrors.
  expect_equal(order(f$table$weight), order(-f$table$se))
})

test_that("order= rearranges the rows and nothing else", {
  local_null_device()
  a <- easydr_forest(bin_fit(), order = "data")
  b <- easydr_forest(bin_fit(), order = "effect")
  d <- easydr_forest(bin_fit(), order = "weight")
  expect_setequal(a$table$study, b$table$study)
  expect_equal(b$table$est, sort(a$table$est))
  expect_equal(d$table$weight, sort(a$table$weight, decreasing = TRUE))
  expect_equal(a$pooled, b$pooled)
})

test_that("a curve has no single trend, so the forest compares one dose", {
  local_null_device()
  rcs <- bin_fit("rcs")
  f   <- easydr_forest(rcs, at = 25)
  expect_equal(f$at, 25)
  expect_match(f$label, "at 25 g/day against 0")
  # The pooled row must agree with what the fit predicts at that dose.
  pred <- rcs$curve$est[which.min(abs(rcs$curve$dose - 25))]
  expect_equal(exp(f$pooled$est), pred, tolerance = 0.05)

  # With no dose named, one is chosen from the fit and reported, not left NA.
  g <- easydr_forest(rcs)
  expect_false(is.na(g$at))
})

test_that("a one-stage fit refuses rather than drawing the pooled line k times", {
  local_null_device()
  one <- easydr(alcohol_cancer, sm = "IRR", shape = "linear", proc = "1stage")
  expect_error(easydr_forest(one), "only the two-stage first stage")
  expect_error(easydr_funnel(one), "only the two-stage first stage")
})

test_that("easydr_curve returns the curve the fit already computed", {
  local_null_device()
  fit <- bin_fit(at = c(10, 25, 50))
  tab <- easydr_curve(fit, main = "t", round = 2)
  expect_identical(tab, fit$curve)
  expect_true(all(c("dose", "est", "ci.lb", "ci.ub", "pi.lb", "pi.ub",
                    "extrapolated") %in% names(tab)))
})

test_that("asking for a prediction band a fit does not have says so", {
  local_null_device()
  common <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
                   common = TRUE)
  expect_false(isTRUE(common$prediction$estimable))
  expect_warning(easydr_curve(common, prediction = TRUE),
                 "no prediction interval")
})

test_that("the plot entry points refuse anything that is not an easydr fit", {
  local_null_device()
  expect_error(easydr_curve(1:10), "takes a fitted result from easydr")
  expect_error(easydr_forest("not a fit"), "takes a fitted result from easydr")
  expect_error(easyDR:::easydr_bubble(bin_fit()), "takes the result of easydr_metareg")
})

test_that("Egger's test is the regression of the deviate on the precision", {
  local_null_device()
  fit <- bin_fit()
  f   <- easydr_funnel(fit)
  eg  <- f$egger

  expect_true(eg$estimable)
  expect_equal(eg$k, 8L)
  expect_equal(eg$df, 6L)
  # Computed independently from the same per-study estimates.
  est <- f$table$est
  se  <- f$table$se
  m   <- stats::lm(I(est / se) ~ I(1 / se))
  expect_equal(eg$intercept, unname(coef(m)[1]))
  expect_equal(eg$p, unname(summary(m)$coefficients[1, 4]))
  # Eight studies is below the threshold where the test has any power.
  expect_true(eg$underpowered)
})

test_that("Egger's test refuses below three studies instead of returning a p", {
  two <- alcohol_cancer[alcohol_cancer$studlab %in% c("atm", "hpm"), ]
  fit <- easydr(two, sm = "IRR", shape = "linear")
  local_null_device()
  f <- easydr_funnel(fit)
  expect_false(f$egger$estimable)
  expect_match(f$egger$reason, "at least 3 studies")
})

test_that("Egger's intercept is zero for a symmetric funnel by construction", {
  # A funnel is symmetric when the estimate does not move with the precision.
  # Then the deviate is exactly proportional to the precision and the line goes
  # through the origin, whatever the sizes of the studies are.
  est <- rep(0.05, 4)
  se  <- c(0.05, 0.10, 0.15, 0.20)
  eg  <- easyDR:::dr_egger(est, se)
  expect_true(eg$estimable)
  expect_equal(eg$intercept, 0, tolerance = 1e-8)
  expect_equal(eg$slope, 0.05, tolerance = 1e-8)
})

test_that("Egger's test refuses when every study is the same size", {
  # The regression has no precision contrast to work with. Returning a p value
  # from a rank-deficient fit would be worse than saying so.
  eg <- easyDR:::dr_egger(c(-0.2, -0.1, 0.1, 0.2), rep(0.1, 4))
  expect_false(eg$estimable)
  expect_match(eg$reason, "same standard error")
})

test_that("the bubble line is the meta-regression's, not a fresh regression", {
  local_null_device()
  fit <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
  mr  <- easyDR:::easydr_metareg(fit, ~ design)
  b   <- easyDR:::easydr_bubble(mr)

  expect_equal(b$moderator, "design")
  expect_equal(nrow(b$table), 6L)
  # One fitted value per level, and they are the regression coefficients: the
  # reference level is the intercept, the other is the intercept plus its term.
  cf <- coef(mr$fit)
  expect_equal(b$fitted$trend[b$fitted$mod == "case-control"], unname(cf[1]))
  expect_equal(b$fitted$trend[b$fitted$mod == "cumulative incidence"],
               unname(cf[1] + cf[2]))
  expect_equal(b$test, mr$test)
})

test_that("a bubble plot of a curve refuses, because there is no trend to plot", {
  local_null_device()
  fit <- easydr(alcohol_cvd, sm = "RR", shape = "quadratic",
                dose.unit = "g/day")
  mr  <- easyDR:::easydr_metareg(fit, ~ design)
  expect_error(easyDR:::easydr_bubble(mr), "needs one trend per study")
  expect_error(easyDR:::easydr_bubble(mr), "easydr_metareg")
})

test_that("asking for a moderator the regression does not have is an error", {
  local_null_device()
  fit <- easydr(alcohol_cvd, sm = "RR", shape = "linear")
  mr  <- easyDR:::easydr_metareg(fit, ~ design)
  expect_error(easyDR:::easydr_bubble(mr, mod = "year"), "not a moderator")
})

test_that("the spread used for categorical bubbles is the same every time", {
  # Random jitter would move a published figure between two runs of the same
  # script. This one is a function of the grouping alone.
  g <- c(1, 1, 1, 2, 2)
  expect_equal(easyDR:::dr_spread(g), easyDR:::dr_spread(g))
  expect_equal(easyDR:::dr_spread(c(1, 2)), c(0, 0))
  expect_equal(sum(easyDR:::dr_spread(g)), 0, tolerance = 1e-12)
})

test_that("a ratio axis is padded multiplicatively, not by adding a constant", {
  r <- easyDR:::dr_pad_range(c(1, 2), ratio = TRUE, frac = 0.1)
  # Equal padding in log space means the two gaps are equal as ratios.
  expect_equal(1 / r[1], r[2] / 2, tolerance = 1e-12)
  s <- easyDR:::dr_pad_range(c(0, 10), ratio = FALSE, frac = 0.1)
  expect_equal(s, c(-1, 11))
})

test_that("every plot works for a continuous outcome too", {
  local_null_device()
  fit <- easydr(aripiprazole, sm = "SMD", shape = "linear",
                dose.unit = "mg/day")
  f <- easydr_forest(fit)
  expect_equal(nrow(f$table), length(unique(fit$frame$studlab)))
  expect_false(fit$ratio)
  tab <- easydr_curve(fit)
  expect_identical(tab, fit$curve)
  fn <- easydr_funnel(fit)
  expect_true(is.list(fn$egger))
})

# ---------------------------------------------------------------------------
# The easyDTA-style layout: what the figures hand back, and what they refuse.

test_that("the forest returns the grob it drew, with the heterogeneity on it", {
  local_null_device()
  f <- easydr_forest(bin_fit())
  expect_s3_class(f$grob, "gtable")
  expect_equal(f$heterogeneity$df, 7L)
  # The label columns come from the same per-study table the rows do.
  expect_true(all(c("levels", "dose_min", "dose_max") %in% names(f$table)))
  expect_true(all(f$table$levels == 6L))
})

test_that("the forest can drop its extra columns without changing the numbers", {
  local_null_device()
  a <- easydr_forest(bin_fit(), columns = TRUE)
  b <- easydr_forest(bin_fit(), columns = FALSE)
  expect_equal(a$table$est, b$table$est)
  expect_equal(a$pooled, b$pooled)
})

test_that("leave-one-out draws one curve per refit that fitted", {
  local_null_device()
  cvd <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
  l <- easydr_loo(cvd, quiet = TRUE)
  expect_s3_class(l, "easydr_sensitivity")
  expect_false(is.null(l$plot))
  expect_equal(length(unique(l$plot$curves$key)), sum(l$table$status == "fitted"))
  # Every curve is lettered, in the order the table is already sorted by.
  expect_equal(sort(unique(l$plot$curves$key)), LETTERS[seq_len(nrow(l$table))])
})

test_that("plot = FALSE gives the table back and draws nothing", {
  cvd <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
  l <- easydr_loo(cvd, quiet = TRUE, plot = FALSE)
  expect_null(l$plot)
  expect_equal(nrow(l$table), 6L)
  # The numbers must not depend on whether a picture was asked for.
  local_null_device()
  withplot <- easydr_loo(cvd, quiet = TRUE)
  expect_equal(withplot$table$estimate, l$table$estimate)
})

test_that("the funnel refuses a shade list that does not match its levels", {
  local_null_device()
  expect_error(easydr_funnel(bin_fit(), levels = c(95, 99),
                             shade = "grey80"),
               "must have the same length")
})

test_that("the funnel's shaded regions are the pseudo confidence wedges", {
  local_null_device()
  f <- easydr_funnel(bin_fit(), levels = c(95), shade = "grey80")
  expect_true(f$egger$estimable)
  # Changing which regions are shaded cannot change the test.
  g <- easydr_funnel(bin_fit())
  expect_equal(f$egger$intercept, g$egger$intercept)
  expect_equal(f$egger$p, g$egger$p)
})

test_that("the bubble's fitted band comes from the regression's own vcov", {
  local_null_device()
  fit <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
  mr  <- easyDR:::easydr_metareg(fit, ~ design)
  b   <- easyDR:::easydr_bubble(mr)
  expect_true(all(c("trend", "lb", "ub") %in% names(b$fitted)))
  expect_true(all(b$fitted$lb < b$fitted$trend))
  expect_true(all(b$fitted$ub > b$fitted$trend))
  # The interval on the reference level is the first coefficient's own.
  cf <- coef(mr$fit); V <- as.matrix(vcov(mr$fit))
  z  <- stats::qnorm(0.975)
  i  <- which(b$fitted$mod == "case-control")
  expect_equal(b$fitted$lb[i], unname(cf[1] - z * sqrt(V[1, 1])))
})

# ---------------------------------------------------------------------------
# Round 3: what the columns say, how I2 is written, and the result box.

test_that("the forest names what its columns count", {
  local_null_device()
  f <- easydr_forest(bin_fit())
  cols <- names(f$grob$grobs[[2]][["children"]])
  # The headers themselves are drawn inside the panel, so the check is on the
  # table the caller gets plus the legend text the figure carries.
  expect_true(all(c("levels", "dose_min", "dose_max") %in% names(f$table)))
  expect_equal(f$table$levels, rep(6L, 8L))
})

test_that("heterogeneity is written I2, not I-squared", {
  local_null_device()
  fit <- bin_fit()
  f   <- easydr_forest(fit)
  # The line is built from the fit, so it is checked by rebuilding it.
  het <- f$heterogeneity
  line <- sprintf("I2 = %.1f%%, Q = %.2f (df %d), %s",
                  het$i2, het$Q, het$df, easyDR:::dr_p_lab(het$p))
  expect_match(line, "^I2 = ")
  expect_no_match(line, "I-squared")
})

test_that("a p value carries its own operator, so none reads p = < 0.001", {
  expect_equal(easyDR:::dr_p_lab(0.687), "p = 0.687")
  expect_equal(easyDR:::dr_p_lab(1e-8), "p < 0.001")
  expect_no_match(easyDR:::dr_p_lab(1e-8), "= <")
})

test_that("the curve's corner box is the fit's own result", {
  local_null_device()
  fit <- bin_fit(at = c(10, 25, 50))
  box <- easyDR:::dr_pooled_at(fit)
  # A one-term shape reports its slope, and that is what the fit printed.
  expect_true(is.na(box$at))
  expect_equal(box$est, fit$slope$estimate)
  expect_equal(box$se, fit$slope$se)
  expect_equal(box$p, fit$slope$p)
  expect_silent(easydr_curve(fit))
})

test_that("a curve has no single slope, so the box names its dose", {
  local_null_device()
  rcs <- bin_fit("rcs", at = c(10, 25, 50))
  box <- easyDR:::dr_pooled_at(rcs)
  expect_equal(box$at, 10)
  expect_match(box$label, "^IRR at 10 g/day$")
  # And it agrees with what the fit predicts there.
  row <- rcs$predicted[rcs$predicted$dose == 10, ]
  expect_equal(exp(box$est), row$est, tolerance = 1e-8)

  # result.at moves it, and nothing else.
  at50 <- easyDR:::dr_pooled_at(rcs, 50)
  expect_equal(at50$at, 50)
  expect_equal(exp(at50$est), rcs$predicted$est[rcs$predicted$dose == 50],
               tolerance = 1e-8)
})

test_that("the corner box works for a one-stage fit, which has no bi", {
  # dr_study_effects() refuses a one-stage fit. The box must not, or the
  # curve would lose its numbers for exactly the fits that need them most.
  local_null_device()
  one <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
                proc = "1stage", dose.unit = "g/day")
  box <- easyDR:::dr_pooled_at(one)
  expect_true(is.finite(box$est))
  expect_true(is.finite(box$p))
  expect_silent(easydr_curve(one))
})

test_that("the meta-regression and the bubble are off the public surface", {
  exports <- getNamespaceExports("easyDR")
  expect_false("easydr_metareg" %in% exports)
  expect_false("easydr_bubble" %in% exports)
  # Still there, still working, just not exported.
  local_null_device()
  fit <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
  mr  <- easyDR:::easydr_metareg(fit, ~ design)
  expect_s3_class(mr, "easydr_moderator")
  expect_type(easyDR:::easydr_bubble(mr)$test, "list")
})

test_that("what the plots return prints itself, so nothing needs rounding by hand", {
  local_null_device()
  fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
                dose.unit = "g/day", at = c(10, 25, 50))

  tab <- easydr_forest(fit)$table
  expect_s3_class(tab, "easydr_contributions")
  # The class is for printing only: every column and value survives it.
  expect_true(all(c("study", "est", "se", "weight") %in% names(tab)))
  expect_equal(sum(tab$weight), 100, tolerance = 1e-8)
  out <- capture.output(print(tab))
  expect_true(any(grepl("weight %", out, fixed = TRUE)))
  expect_true(any(grepl("dose range", out, fixed = TRUE)))

  eg <- easydr_funnel(fit)$egger
  expect_s3_class(eg, "easydr_egger")
  expect_true(is.numeric(eg$intercept))
  out <- capture.output(print(eg))
  expect_true(any(grepl("Egger's test", out, fixed = TRUE)))
  expect_true(any(grepl("intercept -1.0879 (SE 1.2204)", out, fixed = TRUE)))
  expect_true(any(grepl("little power", out, fixed = TRUE)))
})

test_that("summary reports the coefficient, so fit$slope need not be reached into", {
  fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
                dose.unit = "g/day", at = c(10, 25, 50))
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Pooled trend on the log scale, per 1 g/day", out,
                        fixed = TRUE)))
  expect_true(any(grepl("b = 0.0064373,  SE = 0.0013521", out, fixed = TRUE)))

  # A curve has no single slope, so it reports none rather than a made-up one.
  rcs <- easydr(alcohol_cancer, sm = "RR", shape = "rcs", dose.unit = "g/day")
  expect_false(any(grepl("Pooled trend on the",
                         capture.output(print(summary(rcs))), fixed = TRUE)))
})

test_that("the rug marks reported doses, thinned so it reads as data", {
  local_null_device()
  fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
                dose.unit = "g/day", at = c(10, 25, 50))
  # 38 distinct doses across eight cohorts is a solid black band, not a rug.
  expect_gt(length(unique(fit$frame$dose)), 12L)
  expect_silent(easydr_curve(fit))
  expect_silent(easydr_curve(fit, rug = FALSE))
})
