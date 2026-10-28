# Shape resolution and the rules that refuse a shape the data cannot support
# (FR-001 to FR-006, SC-001, SC-007).

test_that("every shape resolves and records its own parameters", {
  d <- crc_table()$dose
  s <- crc_table()$studlab

  lin <- dr_shape("linear", d, studies = s)
  expect_equal(lin$terms, 1L)
  expect_null(lin$knots)

  rcs <- dr_shape("rcs", d, studies = s)
  expect_equal(rcs$terms, 2L)
  expect_length(rcs$knots, 3)
  expect_match(rcs$label, "knots at")

  fp <- dr_shape("fp", d, studies = s)
  expect_equal(fp$powers, c(-1, 0))
  expect_match(fp$label, "powers -1, 0")

  quad <- dr_shape("quadratic", d, studies = s)
  expect_equal(quad$terms, 2L)

  lg <- dr_shape("log", d, studies = s)
  expect_equal(lg$terms, 1L)
  # A dose of zero is present, so the shift must be reported, not assumed.
  expect_match(lg$label, "shifted by")
})

test_that("the spline basis puts the dose in its first column", {
  # This is what makes the Wald test a test of departure from linearity. A
  # basis that mixes the dose into every column answers a different question:
  # on this data the correct test gives p = 0.052 and that substitution gives
  # p = 3.8e-06 (research decision D1).
  d <- crc_table()$dose
  shp <- dr_shape("rcs", d, knots = c(5, 25, 50), studies = crc_table()$studlab)
  B <- shp$basis(d)
  expect_equal(as.numeric(B[, 1]), d)
  expect_equal(colnames(B)[1], "dose")
  expect_equal(ncol(B), 2L)
})

test_that("the quadratic basis also puts the dose first", {
  d <- c(0, 1, 2, 3)
  B <- dr_shape("quadratic", d)$basis(d)
  expect_equal(as.numeric(B[, 1]), d)
  expect_equal(as.numeric(B[, 2]), d^2)
})

test_that("the basis comes from the fitted parameters, not from the grid", {
  # Rebuilding the basis from a prediction grid would silently change the model
  # between fitting and prediction (plan risk 2).
  d <- crc_table()$dose
  shp <- dr_shape("rcs", d, studies = crc_table()$studlab)

  narrow <- shp$basis(c(0, 5, 10))
  wide   <- shp$basis(c(0, 100, 200))
  expect_equal(as.numeric(shp$basis(c(0, 5, 10))[2, ]),
               as.numeric(shp$basis(c(5, 999))[1, ]))
  expect_equal(as.numeric(narrow[1, ]), as.numeric(wide[1, ]))
})

test_that("the fractional-polynomial shift is fixed at resolution", {
  # Left to the engine, the shift and scale are derived from whatever vector is
  # handed over, so fitting and prediction would use two different transforms.
  d <- crc_table()$dose
  shp <- dr_shape("fp", d, powers = c(-1, 0), studies = crc_table()$studlab)
  expect_match(shp$label, "shifted by")
  expect_equal(as.numeric(shp$basis(c(10, 20))[1, ]),
               as.numeric(shp$basis(c(10, 999))[1, ]))
  # A single dose must work, although the engine needs two values of its own.
  expect_equal(nrow(shp$basis(10)), 1L)
})

test_that("V14 refuses a shape no study can inform", {
  tab <- crc_table()
  keep <- unlist(lapply(split(seq_len(nrow(tab)), tab$studlab),
                        function(i) i[1:2]))
  two <- tab[sort(keep), ]
  err <- tryCatch(easydrbin(two, sm = "IRR", shape = "rcs"),
                  error = conditionMessage)
  expect_match(err, "needs a study with at least")
  expect_match(err, "The best any study offers")
  expect_match(err, "does not substitute one shape for another")
})

test_that("V15 refuses more terms than the studies can support", {
  tab <- crc_table()
  two_studies <- tab[tab$studlab %in% unique(tab$studlab)[1:2], ]
  err <- tryCatch(easydrbin(two_studies, sm = "IRR", shape = "rcs",
                            knots = c(5, 10, 25, 50)),
                  error = conditionMessage)
  expect_match(err, "degrees of freedom")
})

test_that("V16 refuses knots outside the range or repeated", {
  tab <- crc_table()
  err <- tryCatch(easydrbin(tab, sm = "IRR", shape = "rcs",
                            knots = c(200, 300, 400)),
                  error = conditionMessage)
  expect_match(err, "outside the observed dose range")

  err2 <- tryCatch(easydrbin(tab, sm = "IRR", shape = "rcs",
                             knots = c(5, 5, 50)),
                   error = conditionMessage)
  expect_match(err2, "repeated values")
})

test_that("a curve differs from a line by one argument only", {
  # SC-001.
  args <- list(data = crc_table(), sm = "IRR", dose.unit = "g/day")
  line  <- do.call(easydrbin, args)
  curve <- do.call(easydrbin, c(args, list(shape = "rcs")))
  expect_s3_class(line, "easydr_result")
  expect_s3_class(curve, "easydr_result")
  expect_equal(line$shape$name, "linear")
  expect_equal(curve$shape$name, "rcs")
})

# ------------------------------------- the cubic and the one-power FP ------

test_that("the cubic is a third-degree polynomial, not the spline", {
  fit <- easydr(alcohol_cancer, sm = "RR", shape = "cubic",
                dose.unit = "g/day")
  expect_equal(fit$shape$name, "cubic")
  expect_equal(fit$shape$terms, 3L)
  expect_match(fit$shape$label, "^cubic in dose")

  # A different model from the restricted cubic spline, which is the only
  # thing the package called "cubic" before.
  rcs <- easydr(alcohol_cancer, sm = "RR", shape = "rcs", dose.unit = "g/day")
  expect_false(identical(coef(fit$fit), coef(rcs$fit)))
  expect_equal(fit$shape$min_levels, 4L)

  # The dose comes first, so the departure-from-linearity test is a Wald test
  # on the two remaining terms.
  expect_equal(colnames(fit$shape$basis(c(1, 2, 3))),
               c("dose", "dose_sq", "dose_cu"))
  expect_true(isTRUE(fit$nonlinearity$estimable))
  expect_equal(fit$nonlinearity$df, 2L)

  # A cubed dose of 60 is near 2e5, so the basis is scaled the way the
  # fractional polynomial already is, and says so.
  expect_equal(fit$shape$scaling$scale, 100)
  expect_match(fit$shape$label, "well conditioned")
  expect_true(all(is.finite(coef(fit$fit))))
})

test_that("a fractional polynomial takes one power or two", {
  one <- easydr(alcohol_cancer, sm = "RR", shape = "fp", powers = -0.5,
                dose.unit = "g/day")
  expect_equal(one$shape$terms, 1L)
  expect_equal(one$shape$powers, -0.5)
  expect_match(one$shape$label, "power -0.5")

  two <- easydr(alcohol_cancer, sm = "RR", shape = "fp", dose.unit = "g/day")
  expect_equal(two$shape$terms, 2L)
  expect_equal(two$shape$powers, c(-1, 0))
  expect_match(two$shape$label, "powers -1, 0")

  # Three is a model this package does not fit, and must not be truncated.
  expect_error(
    easydr(alcohol_cancer, sm = "RR", shape = "fp", powers = c(-1, 0, 1)),
    "3 values were given")

  # A power of zero is the logarithm, following Royston and Altman.
  z <- easydr(alcohol_cancer, sm = "RR", shape = "fp", powers = 0,
              dose.unit = "g/day")
  expect_equal(z$shape$terms, 1L)
  expect_match(z$shape$label, "power 0")
})

test_that("a one-term shape reports its coefficient per unit of its own term", {
  # Only the straight line has a term that is the dose itself. A log fit
  # reporting "Per 1 g/day increase" was describing a step in log dose.
  lin <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
                dose.unit = "g/day")
  lg  <- easydr(alcohol_cancer, sm = "RR", shape = "log", dose.unit = "g/day")
  fp1 <- easydr(alcohol_cancer, sm = "RR", shape = "fp", powers = -0.5,
                dose.unit = "g/day")

  expect_equal(lin$slope$label, "Per 1 g/day increase")
  expect_equal(lg$slope$label, "Per 1 unit increase in log dose")
  expect_equal(fp1$slope$label, "Per 1 unit increase in dose^-0.5")

  # The estimate itself is untouched by the relabelling.
  expect_true(is.finite(lg$slope$estimate))
  expect_false(is.null(lg$slope$se))

  # A shape with more than one term still reports no slope at all.
  expect_null(easydr(alcohol_cancer, sm = "RR", shape = "cubic",
                     dose.unit = "g/day")$slope)
})
