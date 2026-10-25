# Figures for the README images.
#
#   Rscript tutorial/plot.R
#
# Writes PNGs into tutorial/, which tutorial/*.html then embeds. Base graphics
# only: the package itself has no plotting and no graphics dependency, and
# these are development artefacts rather than package code.

suppressMessages(library(easyDR))
suppressMessages(library(dosresmeta))
data(alcohol_crc)

navy <- "#1F2A44"
red  <- "#C0392B"
blue <- "#4E7CA8"
grey <- "#8A939B"
band <- "#EDEFF2"

px <- function(file, w, h) {
  png(file.path("tutorial", file), width = w, height = h, units = "px",
      res = 150, bg = "white")
}

# 1. Exposure levels per study: what dr_viability() counts, and the floor it
#    checks against.
px("fig-levels.png", 1500, 760)
par(mar = c(4.2, 6.4, 1.6, 1.2), family = "sans")
lv <- tapply(alcohol_crc$dose, alcohol_crc$id, function(x) length(unique(x)))
lv <- sort(lv)
bp <- barplot(lv, horiz = TRUE, las = 1, col = band, border = navy,
              xlab = "distinct exposure levels", xlim = c(0, 7),
              cex.names = 1.1, cex.axis = 1.1, cex.lab = 1.15)
abline(v = 3, col = red, lwd = 3, lty = 2)
text(3, max(bp) + 0.6, "two-stage floor", col = red, pos = 4, cex = 1.05,
     font = 2, xpd = NA)
text(lv - 0.18, bp, lv, col = navy, font = 2, cex = 1.1)
dev.off()

# 2. The fitted dose-response, which is the thing curve_shape() then names.
px("fig-curve.png", 1500, 980)
par(mar = c(4.4, 4.6, 1.2, 1.2), family = "sans")
fit <- dosresmeta(formula = logrr ~ dose, type = type, id = id, se = se,
                  cases = cases, n = peryears, data = alcohol_crc)
grid <- data.frame(dose = seq(0, 60, length.out = 120))
pr <- predict(fit, newdata = grid, expo = TRUE)
plot(grid$dose, pr$pred, type = "n", ylim = range(pr$ci.lb, pr$ci.ub),
     xlab = "alcohol, g/day", ylab = "rate ratio", cex.lab = 1.2,
     cex.axis = 1.1, panel.first = grid(col = band, lty = 1))
polygon(c(grid$dose, rev(grid$dose)), c(pr$ci.lb, rev(pr$ci.ub)),
        col = adjustcolor(blue, alpha.f = 0.18), border = NA)
lines(grid$dose, pr$pred, col = navy, lwd = 3)
abline(h = 1, col = grey, lty = 2, lwd = 2)
at <- c(10, 25, 50)
# The reference dose has to be in newdata. Without the 0 row, predict() takes
# the first row as the reference and every marker comes out relative to 10
# g/day instead of to no intake.
pa <- predict(fit, newdata = data.frame(dose = c(0, at)), expo = TRUE)
pa <- pa[pa$dose %in% at, ]
points(pa$dose, pa$pred, pch = 21, bg = red, col = "white", cex = 1.9,
       lwd = 2)
text(pa$dose, pa$pred, sprintf("%.3f", pa$pred), pos = 2, offset = 0.9,
     col = navy, font = 2, cex = 1.05)
legend("topleft", bty = "n", cex = 1.1,
       legend = c("pooled rate ratio", "95% confidence band", "no effect"),
       col = c(navy, adjustcolor(blue, alpha.f = 0.45), grey),
       lwd = c(3, 10, 2), lty = c(1, 1, 2))
dev.off()

# 3. The five labels curve_shape() returns, one panel each.
px("fig-shapes.png", 1900, 470)
par(mfrow = c(1, 5), mar = c(2.6, 2.2, 2.6, 0.8), family = "sans")
shapes <- list(
  "positive association" = c(0, 0.20, 0.40, 0.60, 0.80),
  "Emax"                 = c(0, 0.50, 0.80, 0.85, 0.86),
  "non-monotonic"        = c(0, 0.50, 0.90, 0.40, -0.10),
  "inverse association"  = c(0, -0.20, -0.40, -0.60, -0.80),
  "no association"       = c(0, 0.01, 0.00, 0.01, 0.00)
)
for (nm in names(shapes)) {
  y <- shapes[[nm]]
  plot(0:4, y, type = "n", ylim = c(-0.95, 0.95), axes = FALSE, ann = FALSE)
  abline(h = 0, col = grey, lty = 2, lwd = 2)
  lines(0:4, y, col = navy, lwd = 4)
  points(0:4, y, pch = 21, bg = red, col = "white", cex = 1.5, lwd = 2)
  box(col = band, lwd = 3)
  title(main = nm, cex.main = 1.45, col.main = navy, font.main = 2)
  mtext("dose ->", side = 1, line = 0.8, cex = 0.75, col = grey)
}
dev.off()

cat("wrote tutorial/fig-levels.png, fig-curve.png, fig-shapes.png\n")
