test_that("a monotonic rise is a positive association", {
  expect_equal(curve_shape(0:4, c(0, 0.2, 0.4, 0.6, 0.8)), "positive association")
})

test_that("a rise that flattens is Emax", {
  expect_equal(curve_shape(0:4, c(0, 0.5, 0.8, 0.85, 0.86)), "Emax")
})

test_that("a rise then fall is non-monotonic", {
  expect_equal(curve_shape(0:4, c(0, 0.5, 0.9, 0.4, -0.1)), "non-monotonic")
})

test_that("a monotonic fall is an inverse association", {
  expect_equal(curve_shape(0:4, c(0, -0.2, -0.5, -0.7, -0.9)), "inverse association")
})

test_that("a flat curve is no association", {
  expect_equal(curve_shape(0:4, rep(0, 5)), "no association")
})

test_that("input order does not change the answer", {
  up <- curve_shape(c(4, 1, 0, 3, 2), c(0.8, 0.2, 0, 0.6, 0.4))
  expect_equal(up, "positive association")
})

test_that("too few points is an error", {
  expect_error(curve_shape(0:1, c(0, 1)))
})
