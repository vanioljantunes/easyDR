# Extracted from test-result-at.R:7

# test -------------------------------------------------------------------------
f1 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day")
f2 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", dose.unit = "g/day",
               result.at = 20)
expect_identical(f1$result_at, list(type = "increment", value = 1))
expect_equal(f2$result_at, list(type = "increment", value = 20))
expect_match(easyDR:::dr_pooled_at(f1)$label, "Per 1 g/day increase")
