# Builds every figure and every console block the tutorial pages use.
#
#   Rscript tutorial/figures.R
#
# Figures go to man/figures/, which the README and the HTML pages both read.
# Console output goes to tutorial/_out/, so the text inside the HTML pages is
# copied from a real run rather than typed from memory.

suppressMessages(devtools::load_all("."))

figs <- file.path("man", "figures")
outs <- file.path("tutorial", "_out")
dir.create(figs, showWarnings = FALSE, recursive = TRUE)
dir.create(outs, showWarnings = FALSE, recursive = TRUE)

px <- function(name, w = 1100, h = 760) {
  grDevices::png(file.path(figs, name), width = w, height = h, res = 130)
}

capture_to <- function(name, expr) {
  txt <- utils::capture.output(expr)
  writeLines(txt, file.path(outs, paste0(name, ".txt")))
  invisible(txt)
}

# ---------------------------------------------------------------- the fits --
fit <- easydr(alcohol_cancer, sm = "IRR", shape = "linear",
              dose.unit = "g/day", at = c(10, 25, 50))
rcs <- easydr(alcohol_cancer, sm = "IRR", shape = "rcs",
              dose.unit = "g/day", at = c(10, 25, 50))
cvd <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")

capture_to("easydr", print(fit))
capture_to("summary", print(summary(fit)))
capture_to("rcs", print(rcs))
capture_to("loo", print(easydr_loo(cvd, quiet = TRUE)))
capture_to("metareg", print(easydr_metareg(cvd, ~ design)))
capture_to("select", print(easydr_select(alcohol_cancer, sm = "IRR",
                                         dose.unit = "g/day", quiet = TRUE)))

# ------------------------------------------------------------- the figures --
px("fig-curve.png", 1100, 780)
easydr_curve(fit, main = "Alcohol and colorectal cancer", round = 2)
grDevices::dev.off()

px("fig-curve-rcs.png", 1100, 780)
easydr_curve(rcs, main = "The same studies, fitted as a spline",
             prediction = TRUE, round = 2)
grDevices::dev.off()

px("fig-forest.png", 1180, 800)
easydr_forest(fit)
grDevices::dev.off()

px("fig-funnel.png", 1100, 860)
easydr_funnel(fit, labels = TRUE)
grDevices::dev.off()

px("fig-bubble.png", 1100, 780)
easydr_bubble(easydr_metareg(cvd, ~ design), labels = TRUE)
grDevices::dev.off()

# The five labels easydr() uses to name a curve. Drawn by hand because the
# point is the vocabulary, not any one fit.
shapes <- list(
  "positive association" = function(d) 0.02 * d,
  "Emax"                 = function(d) 0.9 * d / (12 + d),
  "non-monotonic"        = function(d) 0.055 * d - 0.0011 * d^2,
  "inverse association"  = function(d) -0.018 * d,
  "no association"       = function(d) rep(0, length(d))
)
px("fig-shapes.png", 1500, 430)
op <- par(mfrow = c(1, 5), mar = c(4.4, 2.6, 3.4, 0.8), cex = 0.85,
          mgp = c(2.6, 0.8, 0))
d <- seq(0, 60, length.out = 200)
for (nm in names(shapes)) {
  y <- shapes[[nm]](d)
  plot(d, y, type = "l", lwd = 2.4, col = "#1F2A44", xlab = "dose",
       ylab = "", main = nm, ylim = c(-1.15, 1.15), las = 1, yaxt = "n")
  abline(h = 0, lty = 2, col = "#8A939B")
}
par(op)
grDevices::dev.off()

cat("figures written to", figs, "\n")
cat("console blocks written to", outs, "\n")
print(list.files(figs, pattern = "^fig-"))
