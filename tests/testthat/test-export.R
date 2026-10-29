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

# --- several shapes, one figure --------------------------------------------

multi_fit <- function(shapes = c("linear", "quadratic", "rcs")) {
  easydr(alcohol_cancer, sm = "RR", shape = shapes, dose.unit = "g/day")
}

test_that("a shape vector gives one fit per shape, in the order named", {
  m <- multi_fit(c("rcs", "linear"))
  expect_s3_class(m, "easydr_multi")
  expect_equal(m$shapes, c("rcs", "linear"))
  expect_equal(names(m$fits), c("rcs", "linear"))
  for (f in m$fits) expect_s3_class(f, "easydr_result")
  # One shape is still a plain result, not a multi of length one.
  expect_s3_class(multi_fit("linear"), "easydr_result")
})

test_that("a multi fit matches fitting each shape on its own", {
  m <- multi_fit(c("linear", "quadratic"))
  one <- easydr(alcohol_cancer, sm = "RR", shape = "quadratic",
                dose.unit = "g/day")
  expect_equal(coef(m$fits$quadratic$fit), coef(one$fit))
})

test_that("a shape named twice is refused", {
  expect_error(multi_fit(c("linear", "linear")), "more than once")
})

test_that("the comparison table carries what a protocol compares on", {
  m <- multi_fit()
  tab <- easyDR:::dr_multi_table(m)
  expect_equal(tab$Shape, m$shapes)
  expect_true(all(c("Terms", "Deviance", "AIC", "BIC") %in% names(tab)))
  # AIC and BIC are what make non-nested shapes comparable, so they must be
  # real numbers rather than the NA a failed extraction would leave.
  expect_true(all(is.finite(tab$AIC)))
  expect_true(all(is.finite(tab$BIC)))
})

test_that("the single-fit functions refuse a multi fit by name", {
  m <- multi_fit()
  for (f in list(easydr_forest, easydr_funnel, easydr_table)) {
    expect_error(f(m), "holds 3")
  }
  expect_error(easydr_loo(m), "holds 3")
  # The message has to say how to get out of the situation.
  expect_error(easydr_forest(m), "x\\$fits\\$linear")
})

test_that("panels need mult.curves, and mult.curves needs panels", {
  m <- multi_fit()
  expect_error(easydr_curve(m), "mult.curves = TRUE")
  one <- easydr(alcohol_cancer, sm = "RR", shape = "linear")
  expect_error(easydr_curve(one, mult.curves = TRUE), "holds one")
})

test_that("curve.panels chooses the shapes and their order", {
  dir <- withr::local_tempdir()
  m <- multi_fit()
  out <- easydr_curve(m, mult.curves = TRUE,
                      curve.panels = c("rcs", "linear"),
                      file.name = "two", path = dir)
  expect_equal(out$shapes, c("rcs", "linear"))
  expect_equal(out$layout$nrow, 1L)
  expect_equal(out$layout$ncol, 2L)
  expect_true(file.exists(file.path(dir, "two.png")))
  expect_error(easydr_curve(m, mult.curves = TRUE, curve.panels = "emax"),
               "not fitted here")
})

test_that("the layout is the one the author specified", {
  expect_equal(unlist(easyDR:::dr_panel_grid(2)[c("nrow", "ncol")]),
               c(nrow = 1L, ncol = 2L))
  expect_equal(unlist(easyDR:::dr_panel_grid(3)[c("nrow", "ncol")]),
               c(nrow = 1L, ncol = 3L))
  expect_equal(unlist(easyDR:::dr_panel_grid(4)[c("nrow", "ncol")]),
               c(nrow = 2L, ncol = 2L))
  # Five is three over two, so the last cell of the grid stays empty.
  five <- easyDR:::dr_panel_grid(5)
  expect_equal(five$layout[1, ], c(1, 2, 3))
  expect_equal(five$layout[2, ], c(4, 5, NA))
  expect_equal(unlist(easyDR:::dr_panel_grid(6)[c("nrow", "ncol")]),
               c(nrow = 2L, ncol = 3L))
  expect_error(easyDR:::dr_panel_grid(7), "more than this layout covers")
})

test_that("panels share one pair of axes", {
  m <- multi_fit()
  lims <- easyDR:::dr_panel_limits(m$fits, m$ratio)
  every <- range(unlist(lapply(m$fits, function(f)
    range(c(f$curve$ci.lb, f$curve$ci.ub), na.rm = TRUE))))
  # The shared band has to cover every panel's own band, or a curve would be
  # drawn outside the axis it is plotted on.
  expect_lte(lims$ylim[1], every[1])
  expect_gte(lims$ylim[2], every[2])
})

test_that("the table under the panels reports one quantity per shape", {
  m <- multi_fit()
  tab <- easyDR:::dr_panel_table(m$fits, m$shapes)
  expect_equal(names(tab), c(" ", m$shapes))
  expect_equal(nrow(tab), 6L)
  # Every column is the same quantity at the same dose, or the columns would
  # not be comparable.
  expect_match(tab[[1L]][1], "at")
})
