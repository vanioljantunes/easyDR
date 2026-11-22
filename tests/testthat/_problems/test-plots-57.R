# Extracted from test-plots.R:57

# prequel ----------------------------------------------------------------------
local_null_device <- function(env = parent.frame()) {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off(), envir = env)
}
bin_fit <- function(shape = "linear", ...) {
  easydr(alcohol_cancer, sm = "IRR", shape = shape, dose.unit = "g/day", ...)
}

# test -------------------------------------------------------------------------
local_null_device()
rcs <- bin_fit("rcs")
f   <- easydr_forest(rcs, at = 25)
expect_equal(f$at, 25)
expect_match(f$label, "Per 25 g/day increase")
