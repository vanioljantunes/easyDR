# Saving a figure, and the table that goes beside it.
#
# The size a figure is saved at is worked out from the figure itself, so the
# thing worth testing is that the file exists, that its dimensions follow the
# figure rather than a constant, and that a caller who names a size gets that
# size instead. What the picture looks like inside is checked by eye, not here.

cvd_fit <- function(shape = "rcs") {
  easydr(alcohol_cvd, sm = "OR", shape = shape, dose.unit = "g/day",
         at = c(10, 25, 50))
}

png_size <- function(file) {
  # A PNG carries its width and height as big-endian integers at bytes 17-24.
  raw <- readBin(file, "raw", n = 33L)
  int <- function(i) sum(as.integer(raw[i:(i + 3L)]) * 256L^(3:0))
  c(width = int(17L), height = int(21L))
}

test_that("every plot writes the format it was asked for", {
  dir <- withr::local_tempdir()
  fit <- cvd_fit()

  easydr_curve(fit,  file.name = "curve",  path = dir)
  easydr_forest(fit, file.name = "forest", path = dir, filetype = "pdf")
  easydr_funnel(fit, file.name = "funnel", path = dir, filetype = "tiff")
  invisible(easydr_loo(fit, quiet = TRUE, file.name = "loo", path = dir))

  expect_true(file.exists(file.path(dir, "curve.png")))
  expect_true(file.exists(file.path(dir, "forest.pdf")))
  expect_true(file.exists(file.path(dir, "funnel.tiff")))
  expect_true(file.exists(file.path(dir, "loo.png")))
})

test_that("the extension is added once and only once", {
  dir <- withr::local_tempdir()
  fit <- cvd_fit("linear")
  easydr_curve(fit, file.name = "already.png", path = dir)
  expect_true(file.exists(file.path(dir, "already.png")))
  expect_false(file.exists(file.path(dir, "already.png.png")))
})

test_that("a named size is used instead of the measured one", {
  dir <- withr::local_tempdir()
  fit <- cvd_fit("linear")
  easydr_curve(fit, file.name = "fixed", path = dir, width = 6, height = 4)
  # The device is opened at 300 dpi, so six inches is 1800 pixels.
  expect_equal(unname(png_size(file.path(dir, "fixed.png"))), c(1800, 1200))
})

test_that("the measured size follows the figure rather than a constant", {
  dir <- withr::local_tempdir()
  six <- cvd_fit("linear")
  easydr_forest(six, file.name = "six", path = dir)
  eight <- easydr(alcohol_cancer, sm = "RR", shape = "linear")
  easydr_forest(eight, file.name = "eight", path = dir)
  # Eight studies need two more rows than six, so the file is taller. A fixed
  # default would make these identical, which is the mistake this guards.
  expect_gt(png_size(file.path(dir, "eight.png"))[["height"]],
            png_size(file.path(dir, "six.png"))[["height"]])
})

test_that("bad file arguments are refused before anything is drawn", {
  dir <- withr::local_tempdir()
  fit <- cvd_fit("linear")
  expect_error(easydr_curve(fit, file.name = "x", path = file.path(dir, "no")),
               "does not exist")
  expect_error(easydr_curve(fit, file.name = "x", path = dir, width = -2),
               "positive number")
  expect_error(easydr_curve(fit, file.name = c("a", "b"), path = dir),
               "single file name")
  expect_error(easydr_curve(fit, file.name = "x", path = dir,
                            filetype = "svg"))
})

test_that("nothing is written when no file name is given", {
  dir <- withr::local_tempdir()
  fit <- cvd_fit("linear")
  easydr_curve(fit)
  expect_length(list.files(dir), 0L)
})

test_that("easydr_table reports the doses it was given, plus the reference", {
  fit <- cvd_fit()
  tab <- easydr_table(fit, doses = c(10, 25, 50))
  expect_s3_class(tab, "data.frame")
  expect_equal(tab[[1L]], c(0, 10, 25, 50))
  expect_equal(names(tab)[1:2], c("Dose (g/day)", "OR"))
  # The reference row is the comparison every other row is measured against,
  # so it is exactly one on a ratio scale.
  expect_equal(tab[["OR"]][1L], 1)
})

test_that("easydr_table agrees with the fit it came from", {
  fit <- cvd_fit()
  tab <- easydr_table(fit, doses = c(10, 25, 50), digits = 6)
  from_fit <- fit$predicted[fit$predicted$dose == 25, "est"]
  expect_equal(tab[[2L]][tab[[1L]] == 25], round(from_fit, 6))
})

test_that("an extrapolated dose is reported and flagged", {
  fit <- cvd_fit()
  beyond <- max(alcohol_cvd$dose) + 20
  tab <- easydr_table(fit, doses = c(10, beyond))
  expect_true("Extrapolated" %in% names(tab))
  expect_true(tab[["Extrapolated"]][tab[[1L]] == beyond])
  expect_false(tab[["Extrapolated"]][tab[[1L]] == 10])
  # A table with nothing beyond the data does not carry a column of blanks.
  expect_false("Extrapolated" %in% names(easydr_table(fit, doses = 10)))
})

test_that("easydr_table refuses what it cannot report", {
  fit <- cvd_fit()
  expect_error(easydr_table(alcohol_cvd), "fitted result")
  expect_error(easydr_table(fit, doses = "ten"), "numeric vector")
  expect_error(easydr_table(fit, doses = c(10, NA)), "numeric vector")
})

test_that("easydr_table writes a spreadsheet and a PDF", {
  skip_if_not_installed("writexl")
  dir <- withr::local_tempdir()
  fit <- cvd_fit()
  easydr_table(fit, doses = c(10, 25, 50), main = "Table 2",
               filetype = "xlsx", file.name = "t2", path = dir)
  easydr_table(fit, doses = c(10, 25, 50), main = "Table 2",
               filetype = "pdf", file.name = "t2", path = dir)
  expect_true(file.exists(file.path(dir, "t2.xlsx")))
  expect_true(file.exists(file.path(dir, "t2.pdf")))
})
