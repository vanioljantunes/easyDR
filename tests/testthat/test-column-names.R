# The dosresmeta spellings are read without naming them: a `type` column with
# cc / ci / ir, and plain event, n, mean, sd.

test_that("type, event and n are found without arguments", {
  d <- data.frame(studlab = rep(c("A", "B", "C"), each = 3),
                  dose = rep(c(0, 10, 20), 3),
                  event = c(20, 25, 31, 15, 18, 24, 30, 33, 41),
                  n = rep(c(200, 210, 190), 3),
                  type = rep(c("ci", "cc", "ci"), each = 3),
                  ref = rep(c(TRUE, FALSE, FALSE), 3))
  roles <- easyDR:::dr_roles(d, list(), "bin")
  expect_identical(roles$event.e, "event")
  expect_identical(roles$n.e, "n")
  expect_identical(roles$design, "type")
})

test_that("mean, sd and n are found for continuous data", {
  d <- data.frame(studlab = rep("A", 3), dose = c(0, 5, 10),
                  mean = c(1, 2, 3), sd = c(1, 1, 1), n = c(30, 30, 30))
  roles <- easyDR:::dr_roles(d, list(), "cont")
  expect_identical(c(roles$mean.e, roles$sd.e, roles$n.e), c("mean", "sd", "n"))
})

test_that("the .e spelling still wins when both are present", {
  d <- data.frame(studlab = "A", dose = 0, event.e = 1, event = 2, n.e = 10,
                  n = 20)
  roles <- easyDR:::dr_roles(d, list(), "bin")
  expect_identical(c(roles$event.e, roles$n.e), c("event.e", "n.e"))
})

test_that("logrr and se are read as the reported effect, dosresmeta style", {
  # Loaded into its own environment: a fixture of the same name is used by
  # other test files, and data() would otherwise shadow it for them.
  e <- new.env()
  utils::data("alcohol_cvd", package = "dosresmeta", envir = e)
  d <- e$alcohol_cvd
  names(d)[names(d) == "author"] <- "studlab"
  names(d)[names(d) == "cases"] <- "event"
  d$studlab <- paste(d$studlab, d$id)
  d$ref <- d$dose == 0
  te <- d
  te$TE <- exp(d$logrr)
  te$lower <- exp(d$logrr - qnorm(0.975) * d$se)
  te$upper <- exp(d$logrr + qnorm(0.975) * d$se)
  te$logrr <- NULL; te$se <- NULL
  a <- suppressWarnings(easydr(d, sm = "OR", shape = "linear"))
  b <- suppressWarnings(easydr(te, sm = "OR", shape = "linear"))
  expect_equal(coef(a$fit), coef(b$fit), tolerance = 1e-10)
})

test_that("a subgroup on cc / ci codes is labelled with the words", {
  e <- new.env()
  utils::data("alcohol_cvd", package = "dosresmeta", envir = e)
  d <- e$alcohol_cvd
  names(d)[names(d) == "author"] <- "studlab"
  names(d)[names(d) == "cases"] <- "event"
  d$studlab <- paste(d$studlab, d$id)
  d$ref <- d$dose == 0
  r <- suppressWarnings(easydr(d, sm = "OR", shape = "linear", subgroup = "type"))
  expect_setequal(r$subgroup$levels, c("case-control", "cumulative incidence"))
})
