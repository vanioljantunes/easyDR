# The six labels, and what separates each from its neighbour.
#
# The names are the published ones: a dose-response is reported as plateauing
# after an initial rise, as still rising at the maximum dose, or as bell-shaped
# (Leucht et al. 2020). Direction is made explicit for the monotone case and
# the mirror of a bell is named, because a risk that falls and returns means
# the opposite of one that rises and falls.

test_that("a curve still rising at the top dose is monotonic increasing", {
  expect_equal(curve_shape(0:4, c(0, 0.2, 0.4, 0.6, 0.8)),
               "monotonic increasing")
  # Accelerating counts too: what matters is that it has not flattened.
  expect_equal(curve_shape(0:4, c(0, 0.1, 0.3, 0.6, 1.0)),
               "monotonic increasing")
})

test_that("a curve still falling at the top dose is monotonic decreasing", {
  expect_equal(curve_shape(0:4, c(0, -0.2, -0.4, -0.6, -0.8)),
               "monotonic decreasing")
})

test_that("a curve that flattens is a plateau, in either direction", {
  expect_equal(curve_shape(0:4, c(0, 0.5, 0.8, 0.85, 0.86)), "plateau")
  expect_equal(curve_shape(0:4, c(0, -0.5, -0.8, -0.85, -0.86)), "plateau")
})

test_that("plateau is judged against the curve's own steepest step", {
  # A shallow curve and a steep one with the same relative tail earn the same
  # name, which is what makes the label mean anything across outcomes.
  steep   <- c(0, 0.50, 0.80, 0.85, 0.86)
  shallow <- steep / 10
  expect_equal(curve_shape(0:4, steep), curve_shape(0:4, shallow))
  expect_equal(curve_shape(0:4, shallow), "plateau")
})

test_that("the flat floor is absolute, and deliberately so", {
  # The plateau judgement is scale free, but "nothing happened" cannot be: a
  # curve spanning 0.0086 on the log scale is a rate ratio of 1.009 across the
  # whole dose range, and calling that a plateau would dress up a null result.
  tiny <- c(0, 0.50, 0.80, 0.85, 0.86) / 100
  expect_equal(curve_shape(0:4, tiny), "flat")
})

test_that("a rise then a fall is bell-shaped, and the mirror is U-shaped", {
  expect_equal(curve_shape(0:4, c(0, 0.5, 0.9, 0.4, -0.1)), "bell-shaped")
  expect_equal(curve_shape(0:4, c(0, -0.5, -0.9, -0.4, 0.1)), "U-shaped")
})

test_that("bell-shaped and U-shaped are not collapsed into one name", {
  # They mean opposite things. A bell says the effect is lost at high dose; a
  # U says risk returns at high dose. One label for both would be useless.
  bell <- curve_shape(0:4, c(0, 0.5, 0.9, 0.4, -0.1))
  u    <- curve_shape(0:4, c(0, -0.5, -0.9, -0.4, 0.1))
  expect_false(identical(bell, u))
})

test_that("a curve that does nothing is flat", {
  expect_equal(curve_shape(0:4, rep(0, 5)), "flat")
  expect_equal(curve_shape(0:4, c(0, 0.001, -0.001, 0.0005, 0)), "flat")
})

test_that("no label from the old vocabulary survives", {
  # Emax was a model name borrowed as a shape label, and the association
  # wording was not what the dose-response literature publishes.
  grid <- 0:4
  curves <- list(c(0, 0.2, 0.4, 0.6, 0.8), c(0, 0.5, 0.8, 0.85, 0.86),
                 c(0, 0.5, 0.9, 0.4, -0.1), c(0, -0.2, -0.4, -0.6, -0.8),
                 rep(0, 5))
  got <- vapply(curves, function(y) curve_shape(grid, y), character(1))
  expect_false(any(got %in% c("Emax", "positive association",
                              "inverse association", "non-monotonic",
                              "no association")))
  expect_true(all(got %in% c("monotonic increasing", "monotonic decreasing",
                             "plateau", "bell-shaped", "U-shaped", "flat")))
})

test_that("the label describes change, so the scale does not move it", {
  # A monotone curve on the log scale and the same curve exponentiated must
  # get the same name, which is what lets easydr() pass either.
  y <- c(0, 0.2, 0.45, 0.7, 0.95)
  expect_equal(curve_shape(0:4, y), curve_shape(0:4, exp(y)))
})

test_that("the dose order does not matter", {
  d <- c(3, 0, 4, 1, 2)
  y <- c(0.6, 0, 0.8, 0.2, 0.4)
  expect_equal(curve_shape(d, y), "monotonic increasing")
})

test_that("fewer than three doses cannot be described", {
  expect_error(curve_shape(0:1, c(0, 1)))
})
