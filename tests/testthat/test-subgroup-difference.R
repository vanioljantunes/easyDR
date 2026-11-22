# The test for subgroup interaction is the moderator Wald test of Berlin,
# Longnecker and Greenland (1993), run on the overall fit.

red_design_fit <- function(method.tau) {
  data("red_bc", package = "dosresmeta", envir = environment())
  red <- data.frame(
    studlab = paste0(red_bc$author, ", ", red_bc$year, " (", red_bc$id, ")"),
    dose = red_bc$dose, event = red_bc$cases, n = red_bc$n,
    TE = red_bc$rr, lower = red_bc$lrr, upper = red_bc$urr,
    design = ifelse(red_bc$type == "cc", "Case-control", "Cohort"))
  red$lower[red$TE == 1 & red_bc$se == 0] <- NA
  red$upper[is.na(red$lower)] <- NA
  suppressWarnings(easydr(red, sm = "RR", shape = "linear",
                          method.tau = method.tau, subgroup = "design"))
}

test_that("the subgroup test is the moderator Wald test of the overall fit", {
  fit <- red_design_fit("reml")
  t   <- fit$subgroup$test
  expect_true(t$estimable)
  expect_equal(t$df, 1L)
  m <- easydr_metareg(fit, mods = ~ design)
  expect_equal(t$p, m$test$p)
  # Crippa, 2018 (10.1007/s00394-016-1356-0) reports P = 0.02 by design.
  expect_equal(round(t$p, 2), 0.02)
})

test_that("the forest and curve notes write the test", {
  fit <- red_design_fit("reml")
  expect_match(dr_subgroup_test_note(fit$subgroup),
               "^Test for subgroup interaction: X2 = [0-9.]+ \\(df = 1\\), p = ")
  expect_identical(dr_subgroup_test_note(list(test = list(estimable = FALSE))),
                   "Test for subgroup interaction: not available")
})
