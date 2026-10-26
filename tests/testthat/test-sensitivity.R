# Leave-one-out and cumulative analyses. The reference values come from
# helper-anchors.R and every one of them is a captured engine run on the
# bundled alcohol_cvd table, not a published figure. The helper says so too.

test_that("leave-one-out returns one row per study with the influential one first", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  l   <- easydr_loo(fit, quiet = TRUE)

  expect_s3_class(l, "easydr_sensitivity")
  expect_identical(l$kind, "loo")
  expect_equal(nrow(l$table), length(cvd_studies))
  expect_setequal(l$table$label, cvd_studies)
  expect_identical(l$table$label[1], "Bianchi")
  expect_true(all(diff(l$table$abs_shift) <= 0))
})

test_that("each leave-one-out slope matches the captured engine run", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  l   <- easydr_loo(fit, quiet = TRUE)

  expect_equal(l$full$raw, cvd_full_slope, tolerance = 1e-8)
  for (s in cvd_studies) {
    row <- l$table[l$table$label == s, ]
    expect_equal(row$estimate_raw, unname(cvd_loo_slopes[s]), tolerance = 1e-8)
    expect_equal(row$shift, unname(cvd_loo_slopes[s]) - cvd_full_slope,
                 tolerance = 1e-8)
  }
})

test_that("leave-one-out reports the heterogeneity of each refit", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  l   <- easydr_loo(fit, quiet = TRUE)

  expect_true(all(l$table$df == 4L))
  expect_true(all(l$table$Q >= 0))
  expect_true(all(l$table$i2 >= 0 & l$table$i2 <= 100))
  expect_equal(l$fitted, length(cvd_studies))
  expect_equal(l$failed, 0)
  expect_equal(l$tried, nrow(l$table))
})

test_that("the refits themselves come back so a row can be inspected", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  l   <- easydr_loo(fit, quiet = TRUE)

  expect_identical(names(l$refits), l$table$label)
  one <- l$refits[["Bianchi"]]
  expect_s3_class(one, "easydr_result")
  expect_equal(one$slope$estimate, unname(cvd_loo_slopes["Bianchi"]),
               tolerance = 1e-8)
  expect_equal(one$counts$studies, length(cvd_studies) - 1L)
})

test_that("a curve is compared at a named dose and never given a slope", {
  cv <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day",
                  shape = "quadratic", at = c(10, 25))
  l  <- easydr_loo(cv, at = 25, quiet = TRUE)

  expect_identical(l$quantity, "predicted")
  expect_equal(l$dose, 25)
  expect_equal(nrow(l$table), length(cvd_studies))
  expect_true(all(l$table$status == "fitted"))

  out <- paste(capture.output(print(l)), collapse = "\n")
  expect_match(out, "predicted effect at 25")
  expect_false(grepl("the slope, because this shape has one term", out))
})

test_that("a curve refit keeps the knots of the full fit on every row", {
  # The failure this guards against leaves no trace in the output: knots
  # re-placed on five studies instead of six give a different model per row,
  # and every row would still print cleanly.
  cv <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day", shape = "rcs")
  l  <- easydr_loo(cv, quiet = TRUE)

  for (r in l$refits) {
    expect_identical(r$shape$knots, cv$shape$knots)
  }
})

test_that("a refit that cannot be fitted keeps its row and its reason", {
  # A quadratic on three studies leaves no degrees of freedom once one is
  # removed, so every leave-one-out refit fails. The rows must still be there.
  three <- alcohol_cvd[alcohol_cvd$studlab %in% cvd_studies[1:3], ]
  fit   <- easydrbin(three, sm = "RR", shape = "quadratic")
  l     <- easydr_loo(fit, quiet = TRUE)

  expect_equal(nrow(l$table), 3)
  expect_equal(l$tried, 3)
  expect_equal(l$failed, 3)
  expect_true(all(l$table$status == "failed"))
  expect_true(all(nzchar(l$table$reason)))
  expect_match(l$table$reason[1], "degrees of freedom")

  out <- paste(capture.output(print(l)), collapse = "\n")
  expect_match(out, "Refits that could not be fitted")
})

test_that("a two-study fit reports what a one-study refit cannot estimate", {
  two <- alcohol_cvd[alcohol_cvd$studlab %in% cvd_studies[1:2], ]
  fit <- easydrbin(two, sm = "RR")
  l   <- easydr_loo(fit, quiet = TRUE)

  expect_equal(nrow(l$table), 2)
  expect_true(all(l$table$studies == 1))
  expect_true(all(is.na(l$table$Q)))
  expect_true(all(is.na(l$table$i2)))

  out <- paste(capture.output(print(l)), collapse = "\n")
  expect_match(out, "not estimable")
})

test_that("leave-one-out needs at least two studies", {
  one <- alcohol_cvd[alcohol_cvd$studlab == cvd_studies[1], ]
  fit <- easydrbin(one, sm = "RR")
  expect_error(easydr_loo(fit, quiet = TRUE), "at least two studies")
})

test_that("the printed leave-one-out output names what it is measured against", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  out <- paste(capture.output(print(easydr_loo(fit, quiet = TRUE))),
               collapse = "\n")

  expect_match(out, "Leave-one-out sensitivity \\(risk ratio\\)")
  expect_match(out, "Inherited: linear in dose; two-stage reml; reference 0 g/day")
  expect_match(out, "Ordered by the absolute shift")
  expect_match(out, "Full fit:")
  expect_match(out, "carries no p value")
})

test_that("the cumulative sequence accumulates in order and ends at the full fit", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  cm  <- dr_cumul(fit, quiet = TRUE)

  expect_identical(cm$kind, "cumulative")
  expect_identical(cm$table$label, cvd_studies)
  expect_identical(cm$table$studies, seq_along(cvd_studies))
  expect_false(cm$order_stated)

  last <- cm$table[nrow(cm$table), ]
  expect_equal(last$estimate_raw, cvd_full_slope, tolerance = 1e-8)
  expect_equal(last$shift, 0, tolerance = 1e-10)
  expect_equal(last$Q, cvd_full_Q, tolerance = 1e-6)
  expect_identical(last$df, cvd_full_df)
})

test_that("a stated cumulative order is used and reported as stated", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  ord <- rev(cvd_studies)
  cm  <- dr_cumul(fit, order = ord, quiet = TRUE)

  expect_true(cm$order_stated)
  expect_identical(cm$table$label, ord)
  expect_equal(cm$table$estimate_raw[nrow(cm$table)], cvd_full_slope,
               tolerance = 1e-8)

  out <- paste(capture.output(print(cm)), collapse = "\n")
  expect_match(out, "Order: as stated")
})

test_that("the default cumulative order says it was the order in the data", {
  fit <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day")
  out <- paste(capture.output(print(dr_cumul(fit, quiet = TRUE))),
               collapse = "\n")

  expect_match(out, "no order was stated")
  expect_match(out, "The final step reproduces the full fit")
  expect_match(out, "is not a trend")
})

test_that("V22 stops a stated order that does not match the data", {
  fit <- easydrbin(alcohol_cvd, sm = "RR")
  err <- tryCatch(dr_cumul(fit, order = c("Bianchi", "Hokusai")),
                  error = conditionMessage)

  expect_match(err, "not in the data: Hokusai")
  expect_match(err, "in the data but not in the order: Bobak")
  expect_error(dr_cumul(fit, order = c(cvd_studies, "Bianchi")), "repeats")
})

test_that("both analyses run on a continuous outcome too", {
  fit <- easydrcont(aripiprazole, sm = "SMD", dose.unit = "mg/day")
  l   <- easydr_loo(fit, quiet = TRUE)
  cm  <- dr_cumul(fit, quiet = TRUE)

  expect_equal(nrow(l$table), length(unique(aripiprazole$studlab)))
  expect_equal(nrow(cm$table), length(unique(aripiprazole$studlab)))
  expect_equal(cm$table$estimate_raw[nrow(cm$table)], fit$slope$estimate,
               tolerance = 1e-8)
  # Not a ratio measure, so the shift is on the reported scale.
  expect_false(l$ratio)
  out <- paste(capture.output(print(l)), collapse = "\n")
  expect_match(out, "the shift is on the reported scale")
})

test_that("a refit may not override anything about the model", {
  fit <- easydrbin(alcohol_cvd, sm = "RR")
  inh <- dr_inherited(fit)
  expect_error(dr_refit(inh, inh$studies, list(sm = "OR")),
               "may not override")
  expect_error(dr_refit(inh, inh$studies, list(proc = "1stage")),
               "may not override")
})
