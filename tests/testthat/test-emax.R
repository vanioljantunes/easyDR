# The Emax family.
#
# Emax is not linear in its parameters, so easyDR profiles the ones that are
# not: it fits the GLS once per candidate and keeps the best. What is worth
# testing is that the profile actually searches, that a parameter the caller
# names is used instead of being searched for, that the information criteria
# charge for the parameters the search spent, and that a table the family
# cannot be fitted to is refused rather than fitted to something else.

test_that("both members of the family fit and name what they fitted", {
  h <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              dose.unit = "mg/day")
  e <- easydr(aripiprazole, sm = "MD", shape = "emax", dose.unit = "mg/day")
  expect_s3_class(h, "easydr_result")
  expect_s3_class(e, "easydr_result")
  expect_match(h$shape$label, "hyperbolic Emax, ED50")
  expect_match(e$shape$label, "Emax, ED50 .* steepness")
  # One coefficient, Emax. The reference dose absorbs the usual E0 term.
  expect_equal(h$shape$terms, 1L)
  expect_equal(e$shape$terms, 1L)
})

test_that("the hyperbolic form holds its steepness at one", {
  h <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              dose.unit = "mg/day")
  expect_equal(h$shape$steepness, 1)
  expect_equal(h$shape$profiled_params, 1L)
})

test_that("a free steepness costs a second profiled parameter", {
  e <- easydr(aripiprazole, sm = "MD", shape = "emax", dose.unit = "mg/day")
  expect_equal(e$shape$profiled_params, 2L)
  # Two grids multiplied, so the search is larger than the hyperbolic one.
  h <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              dose.unit = "mg/day")
  expect_gt(e$shape$profiled_over, h$shape$profiled_over)
})

test_that("a named ED50 is used rather than searched for", {
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              ed50 = 7, dose.unit = "mg/day")
  expect_equal(f$shape$ed50, 7)
  expect_null(f$shape$profiled_over)
  expect_match(f$shape$label, "ED50 = 7")
})

test_that("the profile searches rather than returning its starting value", {
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              dose.unit = "mg/day")
  obs <- aripiprazole$dose[aripiprazole$dose > 0]
  grid <- stats::quantile(obs, probs = seq(0.05, 0.95, length.out = 20),
                          names = FALSE)
  expect_true(f$shape$ed50 %in% sort(unique(grid[grid > 0])))
  expect_equal(f$shape$profiled_over, length(sort(unique(grid[grid > 0]))))
})

test_that("the basis is the Emax transform and is zero at a dose of zero", {
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              ed50 = 10, dose.unit = "mg/day")
  b <- f$shape$basis(c(0, 10, 1e6))
  expect_equal(unname(b[1, 1]), 0)
  # At the ED50 the curve is half of Emax, which is what ED50 means.
  expect_equal(unname(b[2, 1]), 0.5)
  # And it saturates.
  expect_gt(unname(b[3, 1]), 0.999)
})

test_that("a table with no dose above zero cannot carry an Emax curve", {
  # The shape itself refuses it. A whole table of zero doses is stopped
  # earlier than this, by the data preparation, which is why the guard is
  # checked on the shape rather than through easydr().
  expect_error(easyDR:::dr_shape_hyperbolic_emax(c(0, 0, 0)),
               "doses above zero")
  expect_error(easyDR:::dr_shape_emax(c(0, 0, 0)), "doses above zero")
})

test_that("a bad parameter is refused rather than silently replaced", {
  expect_error(easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
                      ed50 = -5), "above zero")
  expect_error(easydr(aripiprazole, sm = "MD", shape = "emax",
                      ed50 = 10, steepness = 0), "above zero")
})

test_that("the information criteria charge for the profiled parameters", {
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              ed50 = 7, dose.unit = "mg/day")
  ml <- easyDR:::dr_refit_ml(f)
  plain <- stats::AIC(ml$fit)
  # One profiled parameter costs two AIC, which is the whole point: without
  # it the family would be compared at a discount no other shape gets.
  expect_equal(easyDR:::dr_ic(f, "AIC", ml), plain + 2)
  f$shape$profiled_params <- NULL
  expect_equal(easyDR:::dr_ic(f, "AIC", ml), plain)
})

test_that("the comparison refits under ML, because REML would not compare", {
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  expect_equal(f$provenance$method, "reml")
  ml <- easyDR:::dr_refit_ml(f)
  expect_equal(ml$provenance$method, "ml")
  # The same curve, estimated a different way.
  expect_equal(ml$shape$name, f$shape$name)
  # An ML fit is handed straight back rather than refitted again.
  expect_identical(easyDR:::dr_refit_ml(ml), ml)
})

test_that("the family takes its place beside the other shapes", {
  m <- easydr(aripiprazole, sm = "MD",
              shape = c("linear", "hyperbolic_emax"), dose.unit = "mg/day")
  tab <- easyDR:::dr_multi_table(m)
  expect_equal(tab$Shape, c("linear", "hyperbolic_emax"))
  expect_true(all(is.finite(tab$AIC)))
  # Emax saturates, so on a dose-finding table it follows the doses more
  # closely than a straight line does.
  expect_lt(tab$Deviance[tab$Shape == "hyperbolic_emax"],
            tab$Deviance[tab$Shape == "linear"])
})
