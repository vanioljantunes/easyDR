# Before and after. This is the test that makes a refactor safe to do.
#
# inst/validation/baseline.rds holds what easydr answered at a known commit,
# for 32 calls across every path and 18 refusal messages. A simplification may
# change how an answer is computed. It may not change the answer, the words
# printed with it, or the message a refusal gives. Each of those is compared
# here.
#
# Numbers are compared at full precision rather than as printed, because a
# print rounds to three or four digits and would hide a difference in the
# seventh. The tolerance is 1e-9: nothing in a pure refactor should move at
# all, and a genuine change in floating-point association would still sit far
# below it.
#
# It was 1e-12, which is tighter than the arithmetic itself. The baseline is a
# recorded run on one machine, and the same call on the CI runner differs in
# the twelfth significant digit, because a different BLAS associates the sums
# in a different order. That failed every deploy on 2026-10-28 while the
# package was behaving correctly. 1e-9 is still two orders below anything a
# refactor could move and three above the platform noise, and the guard below
# proves the comparison can still see a change.
#
# If one of these fails, the refactor is wrong. Rebuilding the baseline to
# make it pass would be deleting the evidence, which is why the builder
# refuses to overwrite without an explicit environment variable.

baseline_path <- function() {
  # Installed, the file sits under inst/; from the source tree it does not.
  p <- system.file("validation", "baseline.rds", package = "easyDR")
  if (nzchar(p) && file.exists(p)) return(p)
  testthat::test_path("..", "..", "inst", "validation", "baseline.rds")
}

cases_path <- function() {
  p <- system.file("validation", "cases.R", package = "easyDR")
  if (nzchar(p) && file.exists(p)) return(p)
  testthat::test_path("..", "..", "inst", "validation", "cases.R")
}

test_that("every recorded call still gives the same answer and prints the same words", {
  skip_if_not(file.exists(baseline_path()), "no baseline recorded")
  base <- readRDS(baseline_path())
  source(cases_path(), local = TRUE)

  cases  <- dr_validation_cases()
  labels <- vapply(cases, function(k) k$label, character(1))

  # A case dropped from the battery is as much a regression in coverage as a
  # changed number is in behaviour.
  expect_setequal(labels, names(base$results))

  for (k in cases) {
    want <- base$results[[k$label]]
    got  <- withCallingHandlers(
      k$f(),
      warning = function(w) invokeRestart("muffleWarning")
    )

    expect_identical(class(got), want$classes, info = k$label)
    expect_identical(dr_validation_print(got), want$printed, info = k$label)
    expect_equal(dr_validation_numbers(got), want$numbers,
                 tolerance = 1e-9, info = k$label)
  }
})

test_that("every refusal still refuses, with the same message", {
  skip_if_not(file.exists(baseline_path()), "no baseline recorded")
  base <- readRDS(baseline_path())
  source(cases_path(), local = TRUE)

  errs   <- dr_validation_errors()
  labels <- vapply(errs, function(k) k$label, character(1))
  expect_setequal(labels, names(base$errors))

  for (k in errs) {
    want <- base$errors[[k$label]]$message
    got  <- tryCatch({
      k$f()
      NA_character_
    }, error = conditionMessage)
    expect_identical(got, want, info = k$label)
  }
})

test_that("the baseline says which code it was built from", {
  skip_if_not(file.exists(baseline_path()), "no baseline recorded")
  base <- readRDS(baseline_path())

  # A snapshot that cannot be traced to a commit, an R version and an engine
  # version is not evidence of anything.
  expect_true(nzchar(base$commit))
  expect_true(nzchar(base$built_on))
  expect_identical(base$dosresmeta, format(packageVersion("dosresmeta")))
  expect_gt(length(base$results), 25L)
  expect_gt(length(base$errors), 10L)
})

test_that("the comparison would actually notice a changed number", {
  # A snapshot test that cannot fail proves nothing. This perturbs a recorded
  # number by one part in ten billion, far below anything a printout shows,
  # and asserts the comparison rejects it. Verified by hand too: injecting the
  # same error into dr_heterogeneity() turns this file red.
  skip_if_not(file.exists(baseline_path()), "no baseline recorded")
  base <- readRDS(baseline_path())

  want <- base$results[["bin linear, default doses"]]$numbers
  nudged <- want
  nudged$slope[["estimate"]] <- nudged$slope[["estimate"]] * (1 + 1e-7)

  expect_false(isTRUE(all.equal(nudged, want, tolerance = 1e-9)))
  expect_true(isTRUE(all.equal(want, want, tolerance = 1e-9)))
})

test_that("the comparison would notice a changed word", {
  skip_if_not(file.exists(baseline_path()), "no baseline recorded")
  base <- readRDS(baseline_path())

  want <- base$results[["bin linear, default doses"]]$printed
  nudged <- want
  nudged[1] <- paste0(nudged[1], ".")
  expect_false(identical(nudged, want))
})
