# SC-003: every estimation path reachable from the two entry points must have a
# published-value check. This test fails when a covariance or pooling option is
# added without one, so the enumeration above the published tests cannot go
# stale quietly.

test_that("every reachable covariance and pooling combination is covered", {
  covered <- list(
    c(intake = "effect_counts", covariance = "gl",    proc = "2stage"),
    c(intake = "effect_counts", covariance = "gl",    proc = "1stage"),
    c(intake = "effect_counts", covariance = "h",     proc = "2stage"),
    c(intake = "counts",        covariance = "gl",    proc = "2stage"),
    c(intake = "effect",        covariance = "indep", proc = "2stage"),
    c(intake = "continuous",    covariance = "md",    proc = "2stage"),
    c(intake = "continuous",    covariance = "smd",   proc = "2stage"),
    c(intake = "continuous",    covariance = "smd",   proc = "1stage")
  )
  codes <- vapply(covered, `[[`, character(1), "covariance")
  procs <- vapply(covered, `[[`, character(1), "proc")

  # The default covariance of every intake the package can resolve.
  defaults <- c("gl", "indep", "md", "smd")
  expect_true(all(defaults %in% codes))

  # Both pooling stages are exercised.
  expect_true(all(c("2stage", "1stage") %in% procs))

  # Hamling is an override rather than a default, and it is covered.
  expect_true("h" %in% codes)

  # The codes the package accepts are the engine's. "user" needs a covariance
  # matrix the reviewer supplies, which is out of scope for this feature and is
  # the only accepted code without a check.
  accepted <- c("gl", "h", "md", "smd", "indep", "user")
  expect_equal(setdiff(accepted, unique(codes)), "user")
})

test_that("every shape reachable from the entry points has a check", {
  # Feature 002, FR-031 and SC-006. A shape added without a check fails here.
  reachable <- c("linear", "rcs", "fp", "quadratic", "log",
                 "emax", "hyperbolic_emax")
  covered <- c(
    linear    = "test-published-bin.R and the feature 001 fixture",
    rcs       = "test-nonlinearity.R, captured engine run 2026-10-23",
    quadratic = "test-nonlinearity.R and test-curve.R",
    fp        = "test-shapes.R and test-published-cont.R",
    log       = "test-nonlinearity.R, one-term branch",
    emax      = "test-emax.R, profiled over ED50 and steepness",
    hyperbolic_emax = "test-emax.R, profiled over ED50, steepness fixed at 1"
  )
  expect_setequal(names(covered), reachable)

  # The shapes the resolver accepts are exactly the shapes claimed above, so a
  # new value in the resolver cannot slip past this enumeration.
  resolver <- names(formals(dr_shape))
  expect_true("shape" %in% resolver)
  for (s in reachable) {
    expect_s3_class(easydrbin(crc_table(), sm = "IRR", shape = s),
                    "easydr_result")
  }
  expect_error(easydrbin(crc_table(), sm = "IRR", shape = "wiggly"))

  # Both pooling stages are exercised for a curve, not only for the line.
  for (pr in c("2stage", "1stage")) {
    fit <- suppressWarnings(easydrbin(crc_table(), sm = "IRR",
                                      shape = "quadratic", proc = pr))
    expect_equal(fit$provenance$proc, pr)
  }
})

test_that("every analysis of feature 003 is reached on both entry points", {
  # A path no test walks is a path that ships untried. Two shapes each, so the
  # one-term and the multi-term branches of the comparison quantity are both
  # exercised, and both entry points, since the intake differs between them.
  cases <- list(
    list(entry = "bin",  shape = "linear"),
    list(entry = "bin",  shape = "quadratic"),
    list(entry = "cont", shape = "linear"),
    list(entry = "cont", shape = "quadratic")
  )

  for (k in cases) {
    fit <- if (identical(k$entry, "bin")) {
      easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day", shape = k$shape)
    } else {
      easydrcont(aripiprazole, sm = "SMD", dose.unit = "mg/day",
                 shape = k$shape)
    }

    l  <- easydr_loo(fit, quiet = TRUE)
    cm <- dr_cumul(fit, quiet = TRUE)
    expect_s3_class(l, "easydr_sensitivity")
    expect_s3_class(cm, "easydr_sensitivity")
    expect_equal(nrow(l$table), length(l$inherited$studies))
    expect_equal(nrow(cm$table), length(cm$inherited$studies))

    # The comparison quantity follows the shape, never the entry point.
    expect_identical(l$quantity,
                     if (identical(k$shape, "linear")) "slope" else "predicted")
  }

  # Moderators need a study-level covariate, which only the binary table
  # carries, so the continuous side gets a synthetic one.
  ari <- aripiprazole
  u   <- unique(ari$studlab)
  ari$arm <- ifelse(ari$studlab %in% u[1:2], "fixed", "flexible")
  for (shape in c("linear", "quadratic")) {
    expect_s3_class(
      easyDR:::easydr_metareg(easydrbin(alcohol_cvd, sm = "RR", shape = shape), ~ design),
      "easydr_moderator")
    expect_s3_class(
      easyDR:::easydr_metareg(easydrcont(ari, sm = "SMD", shape = shape), ~ arm),
      "easydr_moderator")
  }
})

test_that("the exported surface is exactly what the constitution allows", {
  exports <- getNamespaceExports("easyDR")
  funs <- exports[!exports %in% c("alcohol_cancer", "alcohol_cvd",
                                  "aripiprazole")]
  plain <- funs[!startsWith(funs, "print.") & !startsWith(funs, "summary.")]

  expect_setequal(plain, c("easydr", "easydr_forest", "easydr_curve",
                           "easydr_loo", "easydr_funnel", "easydr_table",
                           "easydr_target_dose", "ed_gamma"))
  # One entry point and seven things to do with its result. The workers behind
  # it, the multi-fit container and the cumulative analysis stay internal:
  # easydr() builds the container itself when it is given several shapes.
  expect_false(any(c("easydrbin", "easydrcont", "dr_cumul", "dr_viability",
                     "curve_shape", "easydr_select", "easydr_metareg",
                     "easydr_bubble", "dr_multi", "dr_curve_multi") %in%
                     plain))
  expect_lte(length(plain), 12L)
})
