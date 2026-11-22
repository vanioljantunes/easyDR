# Extracted from test-sensitivity.R:71

# test -------------------------------------------------------------------------
cv <- easydrbin(alcohol_cvd, sm = "RR", dose.unit = "g/day",
                  shape = "quadratic", at = c(10, 25))
l  <- easydr_loo(cv, at = 25, quiet = TRUE)
expect_identical(l$quantity, "predicted")
expect_equal(l$dose, 25)
expect_equal(nrow(l$table), length(cvd_studies))
expect_true(all(l$table$status == "fitted"))
out <- paste(capture.output(print(l)), collapse = "\n")
expect_match(out, "predicted RR at 25 g/day")
