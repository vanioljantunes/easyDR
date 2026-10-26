suppressMessages(library(dosresmeta))
data(alcohol_crc); data(ari)
alcohol_cancer <- data.frame(
  studlab = as.character(alcohol_crc$id),
  dose    = alcohol_crc$dose,
  TE      = round(exp(alcohol_crc$logrr), 4),
  lower   = round(exp(alcohol_crc$logrr - qnorm(0.975) * alcohol_crc$se), 4),
  upper   = round(exp(alcohol_crc$logrr + qnorm(0.975) * alcohol_crc$se), 4),
  event.e = as.numeric(alcohol_crc$cases),
  n.e     = alcohol_crc$peryears,
  design  = "incidence rate",
  ref     = alcohol_crc$dose == 0,
  stringsAsFactors = FALSE)
aripiprazole <- data.frame(
  studlab = as.character(ari$id),
  dose    = as.numeric(ari$dose),
  mean.e  = ari$y,
  sd.e    = ari$sd,
  n.e     = as.numeric(ari$n),
  ref     = ari$dose == 0,
  stringsAsFactors = FALSE)
# alcohol_cvd mixes designs: four case-control studies and two of cumulative
# incidence. The design is read from `type` rather than assumed, because it is
# what selects each study's covariance approximation. Unlike alcohol_cancer the
# ratios are not rounded, so the bundled table reproduces the reference values
# in specs/003-sensitivity-moderators/research.md exactly and the help examples
# print the anchored numbers.
data(alcohol_cvd)
alcohol_cvd <- data.frame(
  studlab = as.character(alcohol_cvd$author),
  dose    = alcohol_cvd$dose,
  TE      = exp(alcohol_cvd$logrr),
  lower   = exp(alcohol_cvd$logrr - qnorm(0.975) * alcohol_cvd$se),
  upper   = exp(alcohol_cvd$logrr + qnorm(0.975) * alcohol_cvd$se),
  event.e = as.numeric(alcohol_cvd$cases),
  n.e     = as.numeric(alcohol_cvd$n),
  design  = ifelse(alcohol_cvd$type == "cc", "case-control",
                   "cumulative incidence"),
  ref     = alcohol_cvd$dose == 0,
  stringsAsFactors = FALSE)
save(alcohol_cvd, file = "data/alcohol_cvd.rda", compress = "bzip2")
save(alcohol_cancer, file = "data/alcohol_cancer.rda", compress = "bzip2")
save(aripiprazole, file = "data/aripiprazole.rda", compress = "bzip2")
cat("alcohol_cancer:", nrow(alcohol_cancer), "rows,", length(unique(alcohol_cancer$studlab)), "studies\n")
cat("aripiprazole:", nrow(aripiprazole), "rows,", length(unique(aripiprazole$studlab)), "studies\n")
cat("alcohol_cvd:", nrow(alcohol_cvd), "rows,",
    length(unique(alcohol_cvd$studlab)), "studies,",
    length(unique(alcohol_cvd$design)), "designs\n")
