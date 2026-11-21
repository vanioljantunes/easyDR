test_that("easydr(result.at) sets the effect the tables report", {
  f1 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
  f2 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day",
               result.at = 20)
  expect_identical(f1$result_at, "per.unit")
  expect_equal(f2$result_at, 20)
  expect_match(easyDR:::dr_pooled_at(f1)$label, " vs ")
  expect_match(easyDR:::dr_pooled_at(f2)$label, "at 20 g/day")
  expect_match(easyDR:::dr_summary_table(f2, NULL)[1, 1], "at 20")
  # A plot can still override the fit's setting.
  expect_match(easyDR:::dr_pooled_at(f2, "per.unit")$label, " vs ")
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
                      result.at = "two"), "per.unit")
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
  expect_match(lo$label, "1 vs 0")
  f2 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day",
               result.at = 20)
  lo2 <- easydr_loo(f2, plot = FALSE, quiet = TRUE)
  expect_match(lo2$label, "at 20")
})
