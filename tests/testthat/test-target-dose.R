# Target doses.
#
# The point estimate has a closed form for a hyperbolic Emax with a known
# ED50, so it is checked against the algebra rather than against itself. The
# interval has no closed form, so what is checked there is that it carries
# the uncertainty it claims to and refuses when it cannot.

emax_fit <- function(ed50 = NULL) {
  easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax", ed50 = ed50,
         dose.unit = "mg/day")
}

test_that("the target dose matches the algebra for a known ED50", {
  f <- emax_fit(ed50 = 5)
  got <- ed_gamma(f, gamma = c(0.5, 0.8))
  # E(d) is proportional to d / (d + 5), rising throughout, so the largest
  # effect inside the range is at its top and the dose reaching a fraction
  # g of it solves d / (d + 5) = g * hi / (hi + 5).
  hi <- f$counts$dose_max
  want <- vapply(c(0.5, 0.8), function(g) {
    k <- g * hi / (hi + 5)
    5 * k / (1 - k)
  }, numeric(1))
  expect_equal(unname(got$ed), want, tolerance = 1e-3)
  expect_equal(got$ed_max, hi)
  expect_named(got$ed, c("ED50", "ED80"))
})

test_that("a maximum at the end of the range is flagged, not reported", {
  # A curve that is still climbing has no maximum inside the range, only an
  # edge. Reporting the edge as ED max would make a range limit look like a
  # finding.
  f <- emax_fit(ed50 = 5)
  expect_true(ed_gamma(f)$at_edge)
})

test_that("the direction is read from the curve", {
  # Aripiprazole lowers the score, so the measured effect runs one way; the
  # alcohol tables raise the risk and run the other. Both must work without
  # being told which.
  up <- easydr(alcohol_cancer, sm = "RR", shape = "rcs", dose.unit = "g/day")
  expect_gt(ed_gamma(up, 0.5)$ed, 0)
  expect_gt(abs(ed_gamma(emax_fit(5), 0.5)$ed), 0)
})

test_that("gamma is checked", {
  f <- emax_fit(5)
  for (bad in list(0, -1, 1.5, NA_real_, "half", numeric(0))) {
    expect_error(ed_gamma(f, bad), "fraction of the maximum")
  }
  expect_silent(ed_gamma(f, 1))
})

test_that("a multi fit is refused by name", {
  m <- easydr(aripiprazole, sm = "MD",
              shape = c("linear", "hyperbolic_emax"), dose.unit = "mg/day")
  expect_error(ed_gamma(m), "holds 2")
  expect_error(easydr_target_dose(m), "holds 2")
})

test_that("the bootstrap carries where the curve bends", {
  skip_on_cran()
  # With ED50 supplied there is nothing left to vary: this shape has one
  # coefficient, and a target dose depends on the curve rather than on its
  # height, so every draw gives the same answer. With ED50 profiled the
  # draws disagree, and that disagreement is the interval.
  fixed <- easydr_target_dose(emax_fit(ed50 = 5), B = 60, seed = 9)
  expect_true(fixed$degenerate)
  expect_false(fixed$reprofiled)

  prof <- easydr_target_dose(emax_fit(), B = 80, seed = 9)
  expect_false(prof$degenerate)
  expect_true(prof$reprofiled)
  ed50 <- prof$table[prof$table$Quantity == "ED50", ]
  expect_lt(ed50$Lower, ed50$Dose)
  expect_gt(ed50$Upper, ed50$Dose)
})

test_that("the reported table names what it reports", {
  skip_on_cran()
  td <- easydr_target_dose(emax_fit(), gamma = c(0.5, 0.9), B = 60, seed = 3)
  expect_equal(td$table$Quantity, c("ED50", "ED90", "ED max"))
  expect_equal(nrow(td$table), 3L)
  # A percentile interval is read off where the draws fell, so it need not
  # contain the point estimate. When it does not, that is reported rather
  # than hidden.
  contains <- all(td$table$Lower <= td$table$Dose + 1e-8) &&
    all(td$table$Upper >= td$table$Dose - 1e-8)
  expect_equal(td$excludes_point, !contains)
})

test_that("too few draws is refused", {
  expect_error(easydr_target_dose(emax_fit(5), B = 10), "at least")
})

test_that("easydr_table carries the goodness of fit when asked", {
  f <- emax_fit(5)
  plain <- easydr_table(f, doses = c(5, 10))
  expect_null(attr(plain, "gof"))
  with_gof <- easydr_table(f, doses = c(5, 10), gof = TRUE)
  g <- attr(with_gof, "gof")
  expect_s3_class(g, "data.frame")
  expect_equal(g$Statistic,
               c("Deviance (df), p", "R-squared",
                 "Departure from linearity"))
})

test_that("easydr_table carries target doses when asked", {
  skip_on_cran()
  f <- emax_fit(5)
  tab <- easydr_table(f, doses = c(5, 10), target = TRUE, B = 60, seed = 1)
  td <- attr(tab, "target")
  expect_s3_class(td, "data.frame")
  expect_equal(td$Quantity, c("ED50", "ED80", "ED max"))
  # The effects themselves are untouched by asking for the extras.
  expect_equal(tab[[2L]], easydr_table(f, doses = c(5, 10))[[2L]])
})
