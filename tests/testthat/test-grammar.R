# The argument grammar: auto-detection, explicit columns, the one-call promise,
# and no engine vocabulary in the exported signatures (Principle I, FR-002,
# FR-003, FR-006, SC-001).

test_that("column roles are detected from the standard names", {
  fit <- easydrbin(crc_table(), sm = "IRR")
  expect_s3_class(fit, "easydr_result")
  expect_equal(fit$counts$studies, 8)
})

test_that("explicit column names work when the table uses its own", {
  tab <- crc_table()
  names(tab) <- c("trial", "exposure", "kind", "rr", "lo", "hi", "k", "pt")
  fit <- easydrbin(tab, studlab = "trial", dose = "exposure", design = "kind",
                   TE = "rr", lower = "lo", upper = "hi",
                   event.e = "k", n.e = "pt", sm = "IRR")
  expect_equal(fit$counts$studies, 8)
})

test_that("no engine argument name appears in either exported signature", {
  engine_names <- c("logrr", "se", "type", "cases", "id", "Slist", "mod", "v")
  for (f in list(easydrbin, easydrcont)) {
    expect_length(intersect(names(formals(f)), engine_names), 0)
  }
})

test_that("the grammar is the one meta users already know", {
  expected <- c("data", "dose", "studlab", "sm", "method.tau", "level",
                "random", "common")
  expect_true(all(expected %in% names(formals(easydrbin))))
  expect_true(all(c("mean.e", "sd.e", "n.e") %in% names(formals(easydrcont))))
  expect_true(all(c("event.e", "n.e", "TE", "lower", "upper") %in%
                    names(formals(easydrbin))))
})

test_that("a standard table fits with data as the only supplied argument", {
  # SC-001: one call, no reshaping. sm keeps its default.
  tab <- crc_table("counts")
  fit <- easydrbin(tab)
  expect_s3_class(fit, "easydr_result")
  expect_equal(fit$sm, "RR")
})

test_that("no entry point requires more than four arguments", {
  # FR-006: every choice beyond the first four has a documented default.
  for (f in list(easydrbin, easydrcont)) {
    required <- vapply(formals(f), function(a) identical(a, quote(expr = )),
                       logical(1))
    expect_lte(sum(required), 4L)
  }
})

test_that("row order and dose order do not change the fit", {
  # FR-002: long form, no reshaping, sorting or aggregation required.
  tab <- crc_table()
  set.seed(42)
  shuffled <- tab[sample(nrow(tab)), ]
  a <- easydrbin(tab, sm = "IRR")
  b <- easydrbin(shuffled, sm = "IRR")
  expect_equal(a$slope$estimate, b$slope$estimate)
  expect_equal(a$heterogeneity$Q, b$heterogeneity$Q)
  expect_equal(a$counts$categories, b$counts$categories)
})

test_that("a reference flag is honoured and a missing one is reported", {
  tab <- crc_table()
  tab$ref <- tab$dose == 0
  flagged <- easydrbin(tab, sm = "IRR")
  expect_length(flagged$provenance$inferred_ref, 0)

  inferred <- easydrbin(crc_table(), sm = "IRR")
  expect_length(inferred$provenance$inferred_ref, 8)
  expect_match(paste(capture.output(print(inferred)), collapse = " "),
               "reference category inferred")
})

test_that("a reference category above the lowest dose is kept", {
  # FR-005: the reference is per study and is not reset to the minimum. The
  # counts intake is used here because every row carries counts; in the
  # reported-effect intake the published reference row has no limits, so moving
  # the reference elsewhere leaves it unusable and V8 rejects the table.
  tab <- crc_table("counts")
  tab$ref <- FALSE
  for (s in unique(tab$studlab)) {
    rows <- which(tab$studlab == s)
    tab$ref[rows[2]] <- TRUE
  }
  fit <- easydrbin(tab, sm = "IRR")
  expect_gt(fit$reference$curve, 0)
})

test_that("bad arguments are refused with a readable message", {
  expect_error(easydrbin(crc_table(), sm = "IRR", level = 95), "between 0 and 1")
  expect_error(easydrbin(crc_table(), sm = "IRR", proc = "3stage"), "2stage")
  expect_error(easydrbin(crc_table(), sm = "IRR", method.tau = "nonsense"),
               "method.tau")
  expect_error(easydrbin(crc_table(), sm = "IRR", covariance = "wrong"),
               "covariance must be")
  expect_error(easydrbin(as.list(crc_table()), sm = "IRR"), "data frame")
})
