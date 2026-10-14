test_that("a study with two exposure levels is flagged for two-stage", {
  d <- data.frame(
    study = c("A", "A", "A", "B", "B"),
    dose  = c(0, 10, 20, 0, 15),
    cases = c(10, 20, 31, 5, 9),
    n     = c(100, 100, 100, 50, 50)
  )
  v <- dr_viability(d)
  expect_false(v$ok)
  expect_equal(v$summary$levels, c(3L, 2L))
  expect_equal(v$summary$two_stage_ok, c(TRUE, FALSE))
  expect_match(v$problems[1], "one-stage")
})

test_that("a clean three-level dataset passes", {
  d <- data.frame(
    study = rep(c("A", "B"), each = 3),
    dose  = c(0, 5, 10, 0, 8, 16),
    cases = c(5, 9, 14, 7, 11, 18),
    n     = rep(100, 6)
  )
  expect_true(dr_viability(d)$ok)
})

test_that("a missing zero-dose reference is reported", {
  d <- data.frame(
    study = rep("A", 3),
    dose  = c(2, 5, 10),
    cases = c(5, 9, 14),
    n     = rep(100, 3)
  )
  v <- dr_viability(d)
  expect_false(v$ok)
  expect_true(any(grepl("reference", v$problems)))
})

test_that("missing counts block covariance reconstruction", {
  d <- data.frame(
    study = rep("A", 3),
    dose  = c(0, 5, 10),
    cases = c(5, NA, 14),
    n     = rep(100, 3)
  )
  expect_true(any(grepl("covariance", dr_viability(d)$problems)))
})

test_that("a missing column is an error, not a silent pass", {
  expect_error(dr_viability(data.frame(x = 1)), "not found")
})
