# Extracted from test-result-at.R:48

# test -------------------------------------------------------------------------
f <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
lo <- easydr_loo(f, plot = FALSE, quiet = TRUE)
expect_match(lo$label, "Per 1 g/day increase")
f2 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day",
               result.at = 20)
lo2 <- easydr_loo(f2, plot = FALSE, quiet = TRUE)
expect_match(lo2$label, "Per 20 g/day increase")
