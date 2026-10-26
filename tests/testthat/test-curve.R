# The curve as numbers (FR-011 to FR-015, SC-003).

test_that("a curve has no single slope and says so", {
  # SC-003.
  for (s in c("rcs", "fp", "quadratic")) {
    fit <- easydrbin(crc_table(), sm = "IRR", shape = s, dose.unit = "g/day")
    expect_null(fit$slope, info = s)
    out <- paste(capture.output(print(fit)), collapse = " ")
    expect_match(out, "No single slope describes this shape", info = s)
    expect_match(out, "Shape:", info = s)
  }
  line <- easydrbin(crc_table(), sm = "IRR", dose.unit = "g/day")
  expect_false(is.null(line$slope))
  expect_match(paste(capture.output(print(line)), collapse = " "),
               "Per 1 g/day increase")
})

test_that("the curve grid carries both intervals on the reported scale", {
  fit <- easydrbin(crc_table(), sm = "IRR", shape = "rcs", dose.unit = "g/day")
  cv <- fit$curve
  expect_true(all(c("dose", "est", "ci.lb", "ci.ub", "pi.lb", "pi.ub",
                    "is_ref", "extrapolated") %in% names(cv)))
  expect_gte(nrow(cv), 20)
  body <- cv[!cv$is_ref, ]
  expect_true(all(body$ci.lb < body$est))
  expect_true(all(body$ci.ub > body$est))
  expect_equal(cv$est[cv$is_ref], 1)
})

test_that("the prediction interval is at least as wide as the confidence interval", {
  fit <- easydrbin(cvd_table(), sm = "RR", shape = "quadratic")
  cv <- fit$curve[!fit$curve$is_ref, ]
  expect_true(fit$prediction$estimable)
  expect_true(all(cv$pi.lb <= cv$ci.lb))
  expect_true(all(cv$pi.ub >= cv$ci.ub))
})

test_that("the generalised interval agrees with the engine on a quadratic", {
  # Research D5 obtained these limits twice, from predict() and by hand: at
  # dose 25 the quadratic fit gives 1.065058, 0.9445542 to 1.200936.
  fit <- easydrbin(crc_table(), sm = "IRR", shape = "quadratic",
                   at = c(10, 25, 50))
  p <- fit$predicted
  expect_equal(p$est[p$dose == 25], 1.065058, tolerance = 1e-4)
  expect_equal(p$ci.lb[p$dose == 25], 0.9445542, tolerance = 1e-4)
  expect_equal(p$ci.ub[p$dose == 25], 1.200936, tolerance = 1e-4)
})

test_that("a dose beyond the data is marked as extrapolation", {
  fit <- easydrbin(crc_table(), sm = "IRR", shape = "rcs",
                   dose.unit = "g/day", grid = c(0, 10, 20, 40, 200))
  cv <- fit$curve
  expect_true(cv$extrapolated[cv$dose == 200])
  expect_false(any(cv$extrapolated[cv$dose <= 40]))
})

test_that("a grid point at the observed maximum is not called extrapolation", {
  # Rounding a grid endpoint to four significant digits can push it a hair past
  # the observed maximum without leaving the data, which once produced a false
  # extrapolation flag in the default grid.
  fit <- easydrbin(crc_table(), sm = "IRR", shape = "rcs")
  expect_false(any(fit$curve$extrapolated))
})

test_that("a single study reports the prediction interval as not estimable", {
  tab <- crc_table()
  one <- tab[tab$studlab == "atm", ]
  fit <- easydrbin(one, sm = "IRR")
  expect_false(fit$prediction$estimable)
  expect_true(all(is.na(fit$curve$pi.lb)))
})
