# Extracted from test-result-at.R:44

# test -------------------------------------------------------------------------
f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
lo <- easydr_loo(f, plot = FALSE, quiet = TRUE)
expect_match(lo$label, "Per 1 g/day increase")
