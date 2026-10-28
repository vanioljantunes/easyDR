# Builds every figure and every console block the tutorial pages use.
#
#   Rscript tutorial/figures.R
#
# Figures go to man/figures/, which the README and the HTML pages both read.
# Console output goes to tutorial/_out/, so the text inside the HTML pages is
# copied from a real run rather than typed from memory.
#
# Each figure is given the size its own layout needs, in inches. A forest has a
# fixed height per row and a fixed width per column, so a device that does not
# match it leaves a band of white round the figure.

suppressMessages(devtools::load_all("."))

figs <- file.path("man", "figures")
outs <- file.path("tutorial", "_out")
dir.create(figs, showWarnings = FALSE, recursive = TRUE)
dir.create(outs, showWarnings = FALSE, recursive = TRUE)

px <- function(name, w, h, res = 150) {
  grDevices::png(file.path(figs, name), width = w * res, height = h * res,
                 res = res)
}

capture_to <- function(name, expr) {
  txt <- utils::capture.output(expr)
  writeLines(txt, file.path(outs, paste0(name, ".txt")))
  invisible(txt)
}

# ---------------------------------------------------------------- the fits --
fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
              dose.unit = "g/day", at = c(10, 25, 50))
rcs <- easydr(alcohol_cancer, sm = "RR", shape = "rcs",
              dose.unit = "g/day", at = c(10, 25, 50))
cvd <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day")

capture_to("easydr", print(fit))
capture_to("summary", print(summary(fit)))
capture_to("rcs", print(rcs))
capture_to("loo", print(easydr_loo(cvd, quiet = TRUE, plot = FALSE)))
capture_to("shape-error",
           print(try(easydr(alcohol_cancer, sm = "RR"), silent = TRUE)))

# ------------------------------------------------------------- the figures --
px("fig-curve.png", 6.6, 4.4)
easydr_curve(fit, main = "Alcohol and colorectal cancer",
             xlab = "Alcohol intake (g/day)", ylab = "Risk ratio",
             round = 2)
grDevices::dev.off()

px("fig-curve-rcs.png", 6.6, 4.4)
easydr_curve(rcs, main = "The same studies, fitted as a spline",
             prediction = TRUE, round = 2)
grDevices::dev.off()

# The forest knows how wide it is, so it is measured once on a null device
# and the real one is opened at that size. Guessing clips a column.
grDevices::pdf(NULL)
fdim <- easydr_forest(fit)
grDevices::dev.off()
px("fig-forest.png", fdim$width + 0.2, fdim$height + 0.45)
easydr_forest(fit)
grDevices::dev.off()

px("fig-loo.png", 6.8, 6.2)
easydr_loo(cvd, quiet = TRUE)
grDevices::dev.off()

px("fig-funnel.png", 5.6, 6.0)
easydr_funnel(fit, labels = TRUE)
grDevices::dev.off()

# The six labels easydr() uses to name a curve. Drawn by hand because the
# point is the vocabulary, not any one fit.
shapes <- list(
  "monotonic increasing" = function(d) 0.02 * d,
  "monotonic decreasing" = function(d) -0.018 * d,
  "plateau"              = function(d) 0.9 * d / (12 + d),
  "bell-shaped"          = function(d) 0.055 * d - 0.0011 * d^2,
  "U-shaped"             = function(d) -0.055 * d + 0.0011 * d^2,
  "flat"                 = function(d) rep(0, length(d))
)
px("fig-shapes.png", 12, 2.5)
op <- par(mfrow = c(1, 6), mar = c(4.2, 2.4, 3.2, 0.8), cex = 0.9,
          mgp = c(2.5, 0.8, 0))
d <- seq(0, 60, length.out = 200)
for (nm in names(shapes)) {
  # The longest label is wider than one of six panels, so the title is set
  # below full size rather than clipped at both ends.
  plot(d, shapes[[nm]](d), type = "l", lwd = 2.6, col = "#1F2A44",
       xlab = "dose", ylab = "", main = "", ylim = c(-1.15, 1.15),
       las = 1, yaxt = "n")
  title(main = nm, cex.main = 0.82)
  abline(h = 0, lty = 2, col = "#8A939B")
}
par(op)
grDevices::dev.off()

# ------------------------------------------------------- subgroup figures --
# alcohol_cvd carries `design`: four case-control cohorts and two cumulative
# incidence ones. Both levels fit on their own, which is what makes it the
# example rather than alcohol_cancer, where every study shares one design.
sub <- easydr(alcohol_cvd, sm = "RR", shape = "linear", dose.unit = "g/day",
              at = c(10, 20, 30), subgroup = "design")

capture_to("subgroup", print(sub))

grDevices::pdf(NULL)
sdim <- easydr_forest(sub)
grDevices::dev.off()
px("fig-subgroup-forest.png", sdim$width, sdim$height)
easydr_forest(sub, main = "Alcohol and cardiovascular disease, by design")
grDevices::dev.off()

px("fig-subgroup-curve.png", 6.8, 5.6)
easydr_curve(sub, main = "Alcohol and cardiovascular disease, by design",
             xlab = "Alcohol intake (g/day)", ylab = "Risk ratio", round = 2)
grDevices::dev.off()

# ------------------------------------------------- the five shapes, fitted --
# One panel per shape on the same studies and the same axes, so the shapes can
# be compared rather than admired one at a time. result = FALSE drops the
# summary table: at panel size it would be unreadable.
shape_fits <- list(
  linear    = easydr(alcohol_cancer, sm = "RR", shape = "linear",
                     dose.unit = "g/day"),
  log       = easydr(alcohol_cancer, sm = "RR", shape = "log",
                     dose.unit = "g/day"),
  quadratic = easydr(alcohol_cancer, sm = "RR", shape = "quadratic",
                     dose.unit = "g/day"),
  rcs       = easydr(alcohol_cancer, sm = "RR", shape = "rcs",
                     dose.unit = "g/day"),
  fp        = easydr(alcohol_cancer, sm = "RR", shape = "fp",
                     dose.unit = "g/day"))

panels <- lapply(names(shape_fits), function(nm) {
  tabf <- shape_fits[[nm]]$curve
  ggplot2::ggplot(tabf, ggplot2::aes(x = .data$dose)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci.lb, ymax = .data$ci.ub),
                         fill = "grey80") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = "grey35", linewidth = 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = .data$est), colour = "#1F2A44",
                       linewidth = 0.9) +
    ggplot2::scale_y_log10(limits = c(0.9, 2.2)) +
    ggplot2::labs(title = paste0("shape = \"", nm, "\""),
                  x = "g/day", y = if (nm == "linear") "Risk ratio" else "") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold",
                                                      size = 9.5))
})
px("fig-shapes-fitted.png", 12, 2.6)
gridExtra::grid.arrange(grobs = panels, nrow = 1)
grDevices::dev.off()

cat("figures written to", figs, "\n")
print(list.files(figs, pattern = "^fig-"))
