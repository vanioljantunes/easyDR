# SC-008, asserted rather than assumed.
#
# Features 002 and 003 both changed the result object and the printed output.
# Neither was allowed to change what a feature 001 call prints, and the only
# way to know that is to compare against output captured before any of it
# landed. The fixtures were captured by feature 002 and, until this test, were
# never compared with anything.
#
# If one of these fails, the fix belongs in the code. Re-capturing the fixture
# to make it pass would be deleting the evidence.

compare_fixture <- function(fixture, fit) {
  path <- testthat::test_path("fixtures", fixture)
  skip_if_not(file.exists(path), paste("missing fixture", fixture))
  want <- readLines(path, warn = FALSE)
  got  <- capture.output(print(fit))
  expect_identical(got, want)
}

test_that("the feature 001 binary linear output is byte for byte unchanged", {
  compare_fixture(
    "feature001-linear-bin.txt",
    # The call the fixture was captured from, with the default doses. Changing
    # the call here would compare two different things and prove nothing.
    easydrbin(alcohol_cancer, sm = "IRR", dose.unit = "g/day")
  )
})

test_that("the feature 001 continuous linear output is byte for byte unchanged", {
  compare_fixture(
    "feature001-linear-cont.txt",
    easydrcont(aripiprazole, sm = "SMD", dose.unit = "mg/day")
  )
})

test_that("the slots features 001 and 002 readers use are still there", {
  # The object gained slots in feature 003. Gaining them is safe; losing or
  # renaming one would break a script written against an earlier release.
  fit <- easydrbin(alcohol_cancer, sm = "IRR", dose.unit = "g/day")
  for (slot in c("fit", "slope", "sm", "sm_label", "ratio", "dose_unit",
                 "reference", "heterogeneity", "Psi", "shape", "nonlinearity",
                 "gof", "curve", "prediction", "predicted", "counts",
                 "provenance", "frame", "call")) {
    expect_true(slot %in% names(fit), info = slot)
  }
  expect_false(is.null(fit$slope))
  expect_null(easydrbin(alcohol_cancer, sm = "IRR",
                        shape = "quadratic")$slope)
})
