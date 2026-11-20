# Crippa 2018, red meat and bladder cancer, Fig. 2: 13 studies, five of them
# without category counts, pooled per 100 g/day by DerSimonian and Laird.
red_bc_frame <- function() {
  e <- new.env()
  utils::data("red_bc", package = "dosresmeta", envir = e)
  r <- e$red_bc
  lab <- sub(" et al.", "", r$author, fixed = TRUE)
  lab[r$id == 10] <- "Michaud (10)"
  lab[r$id == 11] <- "Michaud (11)"
  d <- data.frame(studlab = lab, dose = r$dose, event = r$cases, n = r$n,
                  TE = r$rr, lower = r$lrr, upper = r$urr,
                  design = factor(ifelse(r$type == "cc", "Case-control", "Cohort"),
                                  c("Cohort", "Case-control")))
  d$lower[d$TE == 1 & r$se == 0] <- NA
  d$upper[is.na(d$lower)] <- NA
  d
}

test_that("studies without counts keep independent contrasts beside the others", {
  f <- easydr(red_bc_frame(), sm = "RR", shape = "linear", dose.unit = "g/day",
              method.tau = "mm")
  expect_match(f$provenance$covariance_label, "independent contrasts for the 5 studies")
  e <- easyDR:::dr_study_effects(f, increment = 100)
  rr <- round(exp(e$table$est), 2)
  names(rr) <- e$table$study
  expect_equal(unname(rr[c("Isa", "Lin", "Tavani", "Michaud (11)")]),
               c(1.94, 2.85, 2.13, 0.94))
  expect_equal(round(exp(c(e$pooled$est, e$pooled$lb, e$pooled$ub)), 2),
               c(1.22, 1.05, 1.41))
  expect_equal(e$label, "Per 100 g/day increase")
})

test_that("a moment estimate of zero heterogeneity is really zero", {
  d <- red_bc_frame()
  f <- easydr(d[d$design == "Cohort", ], sm = "RR", shape = "linear",
              dose.unit = "g/day", method.tau = "mm")
  q <- easyDR:::dr_pooled_at(f, increment = 100)
  expect_equal(round(exp(c(q$est, q$lb, q$ub)), 2), c(1.01, 0.97, 1.06))
})

test_that("increment is one positive number", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")
  expect_error(easyDR:::dr_pooled_at(f, increment = 0), "positive number")
  a <- easyDR:::dr_pooled_at(f)
  b <- easyDR:::dr_pooled_at(f, increment = 10)
  expect_equal(b$est, 10 * a$est)
})
