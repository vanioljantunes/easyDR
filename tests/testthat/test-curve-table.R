fit_rcs <- function() {
  easydr(alcohol_cancer, sm = "IRR", shape = "rcs", dose.unit = "g/day")
}

test_that("the dose table has a Dose column and one row per requested dose", {
  x  <- fit_rcs()
  dt <- easyDR:::dr_dose_table(x, c(10, 25, 50))
  expect_equal(names(dt)[1L], "Dose (g/day)")
  expect_equal(dt[[1L]], c("10", "25", "50"))
  expect_equal(names(dt)[2:4], c("IRR", "95% CI", "p value"))
})

test_that("each table row is the effect against the reference", {
  x  <- fit_rcs()
  dt <- easyDR:::dr_dose_table(x, 25, digits = 6)
  pr <- easyDR:::dr_predicted(x$fit, c(x$reference$curve, 25),
                              x$reference$curve, x$provenance$level, x$sm,
                              x$provenance$single_study, x$shape,
                              range(x$frame$dose))$table
  expect_equal(as.numeric(dt[[2L]]), round(pr$est[pr$dose == 25], 6))
})

test_that("the subtitle reports the tests as p values only", {
  x   <- fit_rcs()
  tok <- easyDR:::dr_curve_subtitle_tokens(x)
  txt <- paste(unlist(lapply(tok, function(g)
    vapply(g, `[[`, character(1), "text"))), collapse = " ")
  expect_match(txt, "Q test p")
  expect_match(txt, "Linearity p")
  expect_match(txt, "Deviance p")
  expect_match(txt, "Restricted cubic spline")
  expect_match(txt, "studies")
  expect_match(txt, "· ")
  expect_false(grepl("X2 =|deviance =|I-squared", txt))
})

test_that("a curve's headline per increment is the linear trend of its studies", {
  x   <- fit_rcs()
  hl  <- easyDR:::dr_curve_headline(x)
  lin <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
                dose.unit = "g/day")
  expect_true(hl$trend)
  expect_match(hl$label, "^Per 1 g/day increase \\(linear trend\\)$")
  expect_equal(hl$pooled$est, lin$slope$estimate, tolerance = 1e-6)
  # A straight line reports its own slope, with no trend label.
  hl_lin <- easyDR:::dr_curve_headline(lin)
  expect_false(hl_lin$trend)
  expect_match(hl_lin$label, "^Per 1 g/day increase$")
  # A dose contrast keeps the curve's own contrast wording.
  hl_pair <- easyDR:::dr_curve_headline(x, c(25, 0))
  expect_false(hl_pair$trend)
  expect_match(hl_pair$label, "^25 vs 0 g/day$")
})

test_that("results.table must be numeric", {
  x <- fit_rcs()
  expect_error(easydr_curve(x, results.table = "a"), "results.table")
})

test_that("a dose above the observed range draws without error", {
  x <- fit_rcs()
  f <- tempfile(fileext = ".png")
  suppressMessages(easydr_curve(x, results.table = c(20, 2 * max(x$frame$dose)),
                                file.name = basename(f), path = dirname(f)))
  expect_true(file.exists(f))
})
