# Which model gets fitted, and whether the output says so truthfully.
#
# Three things were wrong here before 2026-10-25, and all three were the same
# kind of wrong: the code stated something that was not true. A fourth, about
# the shape-selection loop claiming maximum likelihood it had not used, went
# with that loop when it was deleted.
#
#   1. `random` was accepted, documented as affecting the output, stored for
#      refits, and never read. `random = FALSE` returned a bit-identical
#      random-effects fit.
#   2. The printed pooling line said "random effects" whatever was fitted, so
#      a common-effect fit read "two-stage random effects, fixed".
#   3. The prediction interval reported that it needed more studies than
#      coefficients, on a fit with six studies and one coefficient. The real
#      reason was that a common-effect fit has no between-study variance.
#
# None of them had a test, which is why all of them shipped. These are
# those tests.

cvd_fit <- function(...) {
  easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day", ...)
}

test_that("random = FALSE actually changes the model", {
  r <- cvd_fit()
  f <- cvd_fit(random = FALSE)

  expect_false(identical(r$slope$estimate, f$slope$estimate))
  expect_identical(r$provenance$method, "reml")
  expect_identical(f$provenance$method, "fixed")
})

test_that("every route to a common-effect fit gives the same fit", {
  # common = TRUE, random = FALSE and method.tau = "fixed" are three ways of
  # asking for one thing, so they must not disagree.
  a  <- cvd_fit(common = TRUE)
  b  <- cvd_fit(random = FALSE)
  c3 <- cvd_fit(method.tau = "fixed")

  expect_equal(a$slope$estimate, b$slope$estimate, tolerance = 1e-12)
  expect_equal(a$slope$estimate, c3$slope$estimate, tolerance = 1e-12)
  expect_identical(a$provenance$method, "fixed")
  expect_identical(b$provenance$method, "fixed")
  expect_identical(c3$provenance$method, "fixed")
})

test_that("asking for both leaves the common effect, and says which was fitted", {
  # easydr fits one model, unlike meta, which can display two rows from one
  # object. Asking for the common effect wins, and the print names it.
  both <- cvd_fit(common = TRUE, random = TRUE)
  expect_identical(both$provenance$method, "fixed")
  out <- paste(capture.output(print(both)), collapse = "\n")
  expect_match(out, "two-stage common effect")
})

test_that("the printed pooling line names the model that was fitted", {
  rand <- paste(capture.output(print(cvd_fit())), collapse = "\n")
  expect_match(rand, "Pooling: two-stage random effects, reml")
  expect_false(grepl("common effect", rand, fixed = TRUE))

  for (args in list(list(common = TRUE), list(random = FALSE),
                    list(method.tau = "fixed"))) {
    out <- paste(capture.output(print(do.call(cvd_fit, args))),
                 collapse = "\n")
    expect_match(out, "Pooling: two-stage common effect\\.")
    # The contradiction that used to be printed.
    expect_false(grepl("random effects, fixed", out, fixed = TRUE))
    # No tau estimator is named, because none was estimated.
    expect_false(grepl("common effect, ", out, fixed = TRUE))
  }
})

test_that("a common-effect fit gives the real reason for having no prediction interval", {
  f <- cvd_fit(common = TRUE)
  expect_false(f$prediction$estimable)
  expect_match(f$prediction$reason, "common-effect fit")
  expect_match(f$prediction$reason, "no between-study variation")
  # It has six studies and one coefficient, so the old reason was false.
  expect_false(grepl("needs more studies than fitted coefficients",
                     f$prediction$reason, fixed = TRUE))
  expect_equal(f$counts$studies, 6)
})

test_that("a single-study fit still gives its own reason", {
  one <- alcohol_cvd[alcohol_cvd$studlab == "Bianchi", ]
  f <- easydrbin(one, sm = "RR")
  expect_false(f$prediction$estimable)
  expect_match(f$prediction$reason, "one study only")
})

test_that("random and common must be TRUE or FALSE", {
  expect_error(cvd_fit(random = "no"), "`random` must be TRUE or FALSE")
  expect_error(cvd_fit(common = NA), "`common` must be TRUE or FALSE")
  expect_error(cvd_fit(random = c(TRUE, FALSE)),
               "`random` must be TRUE or FALSE")
})

test_that("the continuous entry point behaves the same way", {
  r <- easydrcont(aripiprazole, sm = "SMD", dose.unit = "mg/day")
  f <- easydrcont(aripiprazole, sm = "SMD", dose.unit = "mg/day",
                  random = FALSE)

  expect_identical(r$provenance$method, "reml")
  expect_identical(f$provenance$method, "fixed")
  out <- paste(capture.output(print(f)), collapse = "\n")
  expect_match(out, "two-stage common effect")
  expect_error(easydrcont(aripiprazole, sm = "SMD", common = "yes"),
               "`common` must be TRUE or FALSE")
})

test_that("a sensitivity analysis inherits the common-effect choice", {
  # The model is inherited, so a refit of a common-effect fit must stay one.
  f <- cvd_fit(common = TRUE)
  l <- easydr_loo(f, quiet = TRUE)
  expect_identical(l$inherited$method, "fixed")
  for (r in l$refits) {
    expect_identical(r$provenance$method, "fixed")
  }
})

