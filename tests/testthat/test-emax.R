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

test_that("the profile searches a grid that reaches past the observed doses", {
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              dose.unit = "mg/day")
  obs <- aripiprazole$dose[aripiprazole$dose > 0]
  grid <- exp(seq(log(min(obs) / 10), log(max(obs) * 2), length.out = 25L))
  expect_true(any(abs(f$shape$ed50 - grid) < 1e-8))
  expect_equal(f$shape$profiled_over, 25L)
  # A grid that only spanned the observed doses would pin ED50 to its own
  # edge whenever the curve saturates early, and report the edge as an
  # estimate. It reaches below the smallest dose and above the largest.
  expect_lt(min(grid), min(obs))
  expect_gt(max(grid), max(obs))
})

test_that("a profile that stops at the end of its grid says so", {
  # Forced by handing the profile a grid of one value, which is then both
  # ends of it.
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              dose.unit = "mg/day")
  expect_null(f$shape$ed50_at_boundary)
  shp <- easyDR:::dr_shape_hyperbolic_emax(aripiprazole$dose)
  shp$candidates <- data.frame(ed50 = c(0.001, 0.002), steepness = 1)
  got <- easyDR:::dr_profile_shape(
    shp, f$frame, sm = "MD", intake = f$provenance$intake,
    proc = f$provenance$proc, method = f$provenance$method,
    covariance = f$provenance$covariance,
    method.smd = f$provenance$method.smd)
  expect_true(isTRUE(got$ed50_at_boundary))
  expect_match(got$label, "edge of the search range")
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

# --- the bootstrap ---------------------------------------------------------

test_that("a draw is centred on the model, not on the observed contrasts", {
  # The check that caught two real errors. On a straight line nothing is
  # profiled, so the bootstrap has a closed form to be compared against: the
  # centres must agree. They did not when the draw was taken around the
  # observed values, nor when the study's own baseline was dropped from it.
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  doses <- c(5, 10, 20, 30)
  bs <- easyDR:::dr_boot_curve(f, doses, B = 300, reprofile = FALSE, seed = 11)
  lim <- easyDR:::dr_boot_limits(bs$draws)
  ana <- easyDR:::dr_predicted(f$fit, sort(unique(c(0, doses))), 0,
                               f$provenance$level, f$sm,
                               f$provenance$single_study, f$shape)$table
  ana <- ana[match(doses, ana$dose), ]
  mid_boot <- (lim$lb + lim$ub) / 2
  mid_ana  <- (ana$ci.lb + ana$ci.ub) / 2
  expect_equal(mid_boot, mid_ana, tolerance = 0.05)
})

test_that("every draw is usable on a table the shape suits", {
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  bs <- easyDR:::dr_boot_curve(f, c(10, 20), B = 60, reprofile = FALSE,
                               seed = 3)
  expect_equal(bs$n, 60L)
  expect_equal(dim(bs$draws), c(60L, 2L))
})

test_that("a seed makes the interval reproducible", {
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  a <- easyDR:::dr_boot_curve(f, 10, B = 60, reprofile = FALSE, seed = 5)
  b <- easyDR:::dr_boot_curve(f, 10, B = 60, reprofile = FALSE, seed = 5)
  expect_equal(a$draws, b$draws)
})

test_that("too few draws is refused rather than reported", {
  f <- easydr(aripiprazole, sm = "MD", shape = "linear", dose.unit = "mg/day")
  expect_error(easyDR:::dr_boot_curve(f, 10, B = 10), "at least")
  expect_error(easyDR:::dr_boot_curve(f, 10, B = "many"), "at least")
})

test_that("re-profiling carries the uncertainty in ED50", {
  # Near the bend the position of ED50 is most of what is unknown, so a
  # bootstrap that re-profiles must be wider there than one that holds ED50
  # at the value the real data chose. Out on the plateau it makes little
  # difference, which is the shape of the thing rather than a failure.
  skip_on_cran()
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax",
              dose.unit = "mg/day")
  fixed <- easyDR:::dr_boot_curve(f, 5, B = 200, reprofile = FALSE, seed = 21)
  free  <- easyDR:::dr_boot_curve(f, 5, B = 200, reprofile = TRUE,  seed = 21)
  wf <- diff(as.numeric(unlist(easyDR:::dr_boot_limits(fixed$draws))))
  wr <- diff(as.numeric(unlist(easyDR:::dr_boot_limits(free$draws))))
  expect_gt(wr, wf)
})
