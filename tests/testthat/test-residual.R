test_that("easydr_residual returns one decorrelated residual per non-reference dose", {
  f <- suppressMessages(easydr(alcohol_cvd, sm = "RR", shape = "rcs"))
  tmp <- tempdir()
  r <- easydr_residual(f, file.name = "res-one", path = tmp)
  expect_named(r, c("shape", "studlab", "dose", "residual"))
  expect_equal(nrow(r), sum(!f$frame$is_ref))
  expect_true(all(is.finite(r$residual)))
  expect_true(file.exists(file.path(tmp, "res-one.png")))
})

test_that("easydr_residual draws one panel per shape", {
  f <- suppressMessages(easydr(alcohol_cvd, sm = "RR", shape = "rcs"))
  r <- easydr_residual(f, shapes = c("rcs", "linear"),
                       file.name = "res-two", path = tempdir())
  expect_setequal(unique(r$shape), c("rcs", "linear"))
  expect_equal(dr_data_name(f), "alcohol_cvd")
})

test_that("easydr_residual rejects an empty shapes argument", {
  f <- suppressMessages(easydr(alcohol_cvd, sm = "RR", shape = "rcs"))
  expect_error(easydr_residual(f, shapes = character(0)), "shapes")
})
