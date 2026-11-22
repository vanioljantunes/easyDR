# Extracted from test-nonlinearity.R:43

# test -------------------------------------------------------------------------
for (s in c("rcs", "fp", "quadratic")) {
    fit <- easydrbin(crc_table(), sm = "IRR", shape = s)
    expect_true(fit$nonlinearity$estimable, info = s)
    out <- paste(capture.output(print(fit)), collapse = " ")
    expect_match(out, "Departure from linearity: X2 =", info = s)
    expect_false(grepl("nothing to test it against", out))
  }
