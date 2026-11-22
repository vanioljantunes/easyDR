test_that("easydr(result.at) sets the effect the tables report", {
  f1 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  f2 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day",
               result.at = 20)
  expect_identical(f1$result_at, list(type = "increment", value = 1))
  expect_equal(f2$result_at, list(type = "increment", value = 20))
  expect_match(easyDR:::dr_pooled_at(f1)$label, "Per 1 g/day increase")
  expect_match(easyDR:::dr_pooled_at(f2)$label, "Per 20 g/day increase")
  # A plot can still override the fit's setting.
  expect_match(easyDR:::dr_pooled_at(f2, "per.unit")$label, "Per 1 g/day")
})

test_that("the forest stays per unit whatever easydr() was told", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day",
              result.at = 20)
  e <- easyDR:::dr_study_effects(f, "per.unit")
  expect_equal(e$label, "Per 1 g/day increase")
  expect_identical(formals(easydr_forest)$at, "per.unit")
  expect_null(formals(easydr_funnel)$at)
})

test_that("result.at is per.unit or one dose", {
  expect_error(easydr(alcohol_cvd, sm = "RR", shape = "linear",
                      result.at = "two"), "increment")
  expect_error(easydr(alcohol_cvd, sm = "RR", shape = "linear",
                      overall.per.unit = FALSE, result.at = 5), "numerator")
})

test_that("the funnel table lies across and is headed Summary", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
  out <- withr::with_tempdir(easydr_funnel(f, file.name = "fun", path = ".",
                                           filetype = "png"))
  expect_identical(formals(easydr_funnel)$table, "horizontal")
  tab <- data.frame(Statistic = c("Number of studies", "p value"),
                    Value = c("6", "0.5"))
  g <- easyDR:::dr_summary_grob(tab, "horizontal")
  expect_true("summary-title" %in% g$layout$name)
  expect_s3_class(out$grob, "grob")
})

test_that("leave-one-out follows easydr(result.at)", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  lo <- easydr_loo(f, plot = FALSE, quiet = TRUE)
  expect_match(lo$label, "Per 1 g/day increase")
  f2 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day",
               result.at = 20)
  lo2 <- easydr_loo(f2, plot = FALSE, quiet = TRUE)
  expect_match(lo2$label, "Per 20 g/day increase")
})

test_that("overall.per.unit scales a line by the increment and pairs two doses", {
  f1  <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
  f10 <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day",
                result.at = 10)
  expect_equal(easyDR:::dr_pooled_at(f10)$est, 10 * easyDR:::dr_pooled_at(f1)$est)
  expect_equal(easyDR:::dr_pooled_at(f10)$label, "Per 10 g/day increase")
  fp <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day",
               overall.per.unit = FALSE, result.at = c(30, 10))
  a <- easyDR:::dr_pooled_at(fp)
  expect_equal(a$label, "RR 30 vs 10 g/day")
  b30 <- easyDR:::dr_pooled_at(fp, c(30, 0))$est
  b10 <- easyDR:::dr_pooled_at(fp, c(10, 0))$est
  expect_equal(a$est, b30 - b10)
})

test_that("the forest names its curve and reports the deviance test", {
  expect_equal(easyDR:::dr_forest_deviance(list(available = TRUE, deviance = 24.13, p = 0.17)),
               "Deviance (deviance = 24.13; p = 0.170)")
  expect_match(easyDR:::dr_forest_deviance(list(available = FALSE)), "not available")
  expect_equal(easyDR:::dr_shape_title("rcs"), "Restricted cubic spline")
})
