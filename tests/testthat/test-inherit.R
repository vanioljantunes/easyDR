# A refit must differ from the original fit in exactly one way: which studies
# are in it. Everything here exists because the opposite failure is silent. A
# re-derived spline basis, an auto-detected column resolved differently on a
# subset, a pooling stage quietly switched for comparability: each of those
# produces a table that prints perfectly and answers a different question on
# every row.

test_that("a refit keeps the knots the full fit placed", {
  full <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day", shape = "rcs")
  inh  <- dr_inherited(full)

  for (s in inh$studies) {
    r <- dr_refit(inh, setdiff(inh$studies, s))
    expect_identical(r$shape$knots, full$shape$knots)
    expect_identical(r$shape$name, full$shape$name)
    expect_identical(r$shape$terms, full$shape$terms)
  }
})

test_that("a refit keeps the measure, stage, covariance and reference dose", {
  full <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day",
                    proc = "1stage", method.tau = "ML")
  inh  <- dr_inherited(full)
  r    <- dr_refit(inh, inh$studies[-1])

  expect_identical(r$sm, full$sm)
  expect_identical(r$provenance$proc, full$provenance$proc)
  expect_identical(r$provenance$method, full$provenance$method)
  expect_identical(r$provenance$covariance, full$provenance$covariance)
  expect_identical(r$provenance$level, full$provenance$level)
  expect_identical(r$reference$curve, full$reference$curve)
  expect_identical(r$dose_unit, full$dose_unit)
})

test_that("a refit keeps the fractional-polynomial shift and scale", {
  # The powers alone do not define an fp basis: the shift and the conditioning
  # scale do too, and both were fixed from the full dose range. Re-deriving
  # them on a subset would move the design without moving the powers, which is
  # the kind of change no printed label would show.
  full <- easydrbin(alcohol_cvd, sm = "RR", shape = "fp", powers = c(-1, 0))
  inh  <- dr_inherited(full)
  r    <- dr_refit(inh, inh$studies[-1])

  expect_identical(r$shape$powers, full$shape$powers)
  expect_identical(r$shape$scaling, full$shape$scaling)
  expect_identical(r$shape$label, full$shape$label)
})

test_that("a refit resolves the same columns, not whatever the subset suggests", {
  full <- easydrbin(alcohol_cvd, sm = "RR")
  inh  <- dr_inherited(full)
  r    <- dr_refit(inh, inh$studies[-1])

  expect_identical(r$provenance$intake, full$provenance$intake)
  expect_identical(r$args$studlab, full$args$studlab)
  expect_identical(r$args$dose, full$args$dose)
  expect_identical(r$args$TE, full$args$TE)
})

test_that("a refit repeats the fitted shape and never re-derives it", {
  # A curve's knots were placed on the full set of doses. Re-deriving them on
  # a subset would answer a different question on every row, so the shape
  # travels as an object and the refit uses the fitted basis.
  fit <- easydr(alcohol_cvd, sm = "RR", shape = "quadratic",
                dose.unit = "g/day")
  inh <- dr_inherited(fit)
  r   <- dr_refit(inh, inh$studies[-1])

  expect_identical(r$shape$name, fit$shape$name)
  expect_identical(r$shape$terms, fit$shape$terms)
  expect_identical(r$shape$label, fit$shape$label)
})

test_that("the subset is the only thing that changes", {
  full <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  inh  <- dr_inherited(full)
  keep <- inh$studies[-1]
  r    <- dr_refit(inh, keep)

  expect_equal(r$counts$studies, length(keep))
  expect_false(any(as.character(r$frame$studlab) == inh$studies[1]))
  expect_setequal(unique(as.character(r$frame$studlab)), keep)
})

test_that("inheritance works the same way for continuous outcomes", {
  full <- easydrcont(aripiprazole, sm = "SMD", dose.unit = "mg/day",
                     shape = "quadratic")
  inh  <- dr_inherited(full)
  r    <- dr_refit(inh, inh$studies[-1])

  expect_identical(r$sm, full$sm)
  expect_identical(r$shape$name, full$shape$name)
  expect_identical(r$provenance$method, full$provenance$method)
  expect_identical(r$args$method.smd, full$args$method.smd)
})

test_that("a refit of a shape the subset cannot support fails rather than simplifying", {
  # A shape needing more terms than the subset leaves room for. The refusal is
  # the point: a row reporting a simpler shape under the label of the fitted
  # one would be worse than a row reporting nothing.
  full <- easydrbin(alcohol_cvd, sm = "RR", shape = "quadratic")
  inh  <- dr_inherited(full)
  expect_error(dr_refit(inh, inh$studies[1:2]), "degrees of freedom")
})

test_that("a result without its table says so instead of guessing", {
  full <- easydrbin(alcohol_cvd, sm = "RR")
  full$data <- NULL
  expect_error(dr_inherited(full), "not available")
})

test_that("these analyses refuse anything that is not a fitted result", {
  expect_error(dr_inherited(alcohol_cvd), "fitted result")
})
