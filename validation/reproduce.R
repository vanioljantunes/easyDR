# easyDR: reproduce the published alcohol and colorectal cancer analysis.
# Eight prospective cohorts, analysed in Orsini et al. (2012) and bundled as
# alcohol_cancer. Paste the whole file into R; no data files are needed.
# install.packages("remotes"); remotes::install_github("vanioljantunes/easyDR")

library(easyDR)

# One call. No formula argument and no logrr column: sm = "RR" tells easyDR to
# read the counts and the reported risk ratios, derive the contrasts, and
# reconstruct the within-study covariances.
fit <- easydr(alcohol_cancer, sm = "RR", shape = "linear",
              dose.unit = "g/day", at = c(10, 25, 50))
fit

# The coefficient on the log scale, the goodness of fit, and what each cohort
# contributed. This is the number to quote in a methods section.
summary(fit)

# The rows behind the forest plot, with weights. The figure is drawn and the
# table prints itself; the stored values keep full precision.
easydr_forest(fit)$table

# Small-study effects. Egger's regression, printed as a test result with its
# power warning. Deeks' test belongs to diagnostic test accuracy.
easydr_funnel(fit)$egger

# The two checks easyDR ran before and after the model.
fit$curve_label
fit$viability

# The curve, with every argument the README figure used.
easydr_curve(fit,
             main  = "Alcohol and colorectal cancer",
             xlab  = "Alcohol intake (g/day)",
             ylab  = "Risk ratio",
             round = 2)
