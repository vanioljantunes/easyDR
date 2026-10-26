# One assertion per validation rule, V1 to V13. Every message must name the
# study and the rows involved, so a reviewer can find the problem in their own
# spreadsheet (FR-018, SC-006).

test_that("V1 stops on a non-numeric dose", {
  tab <- crc_table()
  tab$dose <- as.character(tab$dose)
  expect_error(easydrbin(tab, sm = "IRR"), "not numeric")
  expect_error(easydrbin(tab, sm = "IRR"), "category label")
})

test_that("V2 stops when a study has no reference row or several", {
  tab <- crc_table()
  tab$ref <- FALSE
  expect_error(easydrbin(tab, sm = "IRR"), "0 reference rows")

  tab$ref <- TRUE
  err <- tryCatch(easydrbin(tab, sm = "IRR"), error = conditionMessage)
  expect_match(err, "reference rows")
  expect_match(err, "atm")
  expect_match(err, "rows")
})

test_that("V3 stops when every row shares one dose", {
  tab <- crc_table()
  tab$dose <- 0
  expect_error(easydrbin(tab, sm = "IRR"), "single dose value")
})

test_that("V4 stops on duplicated study and dose pairs", {
  tab <- crc_table()
  tab <- rbind(tab, tab[1, ])
  err <- tryCatch(easydrbin(tab, sm = "IRR"), error = conditionMessage)
  expect_match(err, "duplicated study and dose")
  expect_match(err, "atm")
})

test_that("V5 stops on impossible counts", {
  tab <- crc_table("counts")
  tab$design <- "cumulative incidence"
  tab$event.e[3] <- tab$n.e[3] + 1
  err <- tryCatch(easydrbin(tab, sm = "RR"), error = conditionMessage)
  expect_match(err, "impossible counts")
  expect_match(err, "atm")
  expect_match(err, "3")
})

test_that("V6 warns and records the zero-cell correction", {
  tab <- crc_table("counts")
  tab$event.e[2] <- 0
  expect_warning(fit <- easydrbin(tab, sm = "IRR"), "zero cases")
  expect_true(length(fit$provenance$continuity) > 0L)
  expect_match(paste(capture.output(print(fit)), collapse = " "),
               "zero-cell correction")
})

test_that("V7 stops when a reported effect sits outside its own limits", {
  tab <- crc_table("effect")
  tab$TE[2] <- tab$upper[2] * 2
  err <- tryCatch(easydrbin(tab, sm = "IRR"), error = conditionMessage)
  expect_match(err, "outside its own")
  expect_match(err, "atm")
})

test_that("V8 stops when limits are present for some rows of a study only", {
  tab <- crc_table("effect")
  tab$lower[3] <- NA_real_
  err <- tryCatch(easydrbin(tab, sm = "IRR"), error = conditionMessage)
  expect_match(err, "some dose")
  expect_match(err, "atm")
})

test_that("V9 warns when per-study dose ranges differ wildly", {
  tab <- crc_table()
  tab$dose[tab$studlab == "hpm"] <- tab$dose[tab$studlab == "hpm"] / 1000
  expect_warning(easydrbin(tab, sm = "IRR"), "factor of 100")
  expect_warning(easydrbin(tab, sm = "IRR"), "does not convert units")
})

test_that("V10 reports heterogeneity and the prediction interval as not estimable", {
  tab <- crc_table()
  tab <- tab[tab$studlab == "atm", ]
  fit <- easydrbin(tab, sm = "IRR")
  expect_false(fit$heterogeneity$estimable)
  expect_false(fit$prediction$estimable)
  out <- paste(capture.output(print(fit)), collapse = " ")
  expect_match(out, "not estimable")
  # Not estimable, never reported as zero.
  expect_true(all(is.na(fit$predicted$pi.lb)))
})

test_that("V11 records how many studies could inform a nonlinearity test", {
  fit <- easydrbin(crc_table(), sm = "IRR")
  expect_equal(fit$nonlinearity$studies, 8)
  expect_equal(fit$nonlinearity$studies_with_three_levels, 8)

  two_levels <- crc_table()
  keep <- unlist(lapply(split(seq_len(nrow(two_levels)), two_levels$studlab),
                        function(i) i[1:2]))
  fit2 <- easydrbin(two_levels[sort(keep), ], sm = "IRR")
  expect_equal(fit2$nonlinearity$studies_with_three_levels, 0)
  expect_false(fit2$nonlinearity$estimable)
})

test_that("V12 lists the names it looked for and the names present", {
  tab <- crc_table()
  names(tab)[names(tab) == "dose"] <- "exposure"
  err <- tryCatch(easydrbin(tab, sm = "IRR"), error = conditionMessage)
  expect_match(err, "Looked for")
  expect_match(err, "exposure")
  expect_match(err, "Columns present")
})

test_that("V13 stops when easydr would compute a measure the counts cannot give", {
  # The rule applies where easydr computes the effect itself. From
  # case-control counts only an odds ratio is available, and from counts over
  # totals only a risk ratio.
  counts <- crc_table("counts")
  expect_error(easydrbin(counts, sm = "OR"), "case-control")

  cc <- crc_table("counts")
  cc$design <- "case-control"
  expect_error(easydrbin(cc, sm = "RR"), "only give an odds ratio")
})

test_that("V13 requires a person-time denominator for a rate ratio either way", {
  # A denominator is a fact about the data, not a reporting choice, so this
  # check does not depend on who computed the effect.
  cc <- crc_table("counts")
  cc$design <- "case-control"
  expect_error(easydrbin(cc, sm = "IRR"), "person-time")

  rep <- crc_table("effect")
  rep$design <- "cumulative incidence"
  expect_error(easydrbin(rep, sm = "IRR"), "person-time")
})

test_that("V13 allows a reported effect whose design differs from sm", {
  # Where the authors reported the effect, sm names the scale their ratios are
  # pooled on and the design only selects the covariance approximation.
  rep <- crc_table("effect_counts")
  rep$design <- "case-control"
  expect_s3_class(easydrbin(rep, sm = "RR"), "easydr_result")
})

test_that("a mixed-design table pools and says so", {
  # alcohol_cvd is four case-control studies and two cumulative-incidence
  # ones. Pooling that mixture is standard in this literature, so easydr fits
  # it and names the mixture rather than letting one sm imply a shared design.
  fit <- easydrbin(cvd_table(), sm = "RR", dose.unit = "g/day")
  expect_s3_class(fit, "easydr_result")
  out <- paste(capture.output(print(fit)), collapse = "
")
  expect_match(out, "designs differ across studies")
  expect_match(out, "4 case-control")
  expect_match(out, "2 cumulative incidence")
})
