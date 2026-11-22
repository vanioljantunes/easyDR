# Extracted from test-plots.R:335

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
rcs <- bin_fit("rcs", at = c(1, 10, 25, 50))
box <- easyDR:::dr_pooled_at(rcs)
expect_equal(box$at, 1)
expect_match(box$label, "^Per 1 g/day increase$")
