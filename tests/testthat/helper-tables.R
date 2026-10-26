# Test tables built from the engine's own published datasets, so every
# expectation traces to a source a reader can check.

crc_table <- function(intake = c("effect_counts", "counts", "effect")) {
  intake <- match.arg(intake)
  data(alcohol_crc, package = "dosresmeta", envir = environment())
  d <- get("alcohol_crc", envir = environment())
  tab <- data.frame(
    studlab = as.character(d$id),
    dose    = d$dose,
    design  = "incidence rate",
    stringsAsFactors = FALSE
  )
  if (intake %in% c("effect_counts", "effect")) {
    tab$TE    <- exp(d$logrr)
    tab$lower <- exp(d$logrr - stats::qnorm(0.975) * d$se)
    tab$upper <- exp(d$logrr + stats::qnorm(0.975) * d$se)
  }
  if (intake %in% c("effect_counts", "counts")) {
    tab$event.e <- d$cases
    tab$n.e     <- d$peryears
  }
  tab
}

# alcohol_cvd is four case-control studies and two cumulative-incidence ones.
# An earlier version of this helper labelled all six as cumulative incidence,
# which fed the wrong `type` to the engine for four of them. The design comes
# from the data.
cvd_table <- function() {
  data(alcohol_cvd, package = "dosresmeta", envir = environment())
  d <- get("alcohol_cvd", envir = environment())
  data.frame(
    studlab = as.character(d$author),
    dose    = d$dose,
    TE      = exp(d$logrr),
    lower   = exp(d$logrr - stats::qnorm(0.975) * d$se),
    upper   = exp(d$logrr + stats::qnorm(0.975) * d$se),
    event.e = d$cases,
    n.e     = d$n,
    design  = ifelse(d$type == "cc", "case-control", "cumulative incidence"),
    stringsAsFactors = FALSE
  )
}

ari_table <- function() {
  data(ari, package = "dosresmeta", envir = environment())
  d <- get("ari", envir = environment())
  data.frame(
    studlab = as.character(d$id),
    dose    = d$dose,
    mean.e  = d$y,
    sd.e    = d$sd,
    n.e     = d$n,
    stringsAsFactors = FALSE
  )
}
