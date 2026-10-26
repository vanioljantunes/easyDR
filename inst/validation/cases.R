# The battery of calls that must keep answering identically.
#
# This file defines the cases once and is sourced by both the snapshot builder
# and the test that compares against the snapshot, so the two cannot drift
# apart on which calls they cover.
#
# A refactor is allowed to change how easydr computes an answer. It is not
# allowed to change the answer, or the words easydr prints with it. Every case
# below therefore records both: the printed output as text, and the numbers
# underneath it at full precision, since a print rounds to three or four
# digits and would hide a difference in the seventh.

dr_validation_cases <- function() {
  cvd <- easyDR::alcohol_cvd
  crc <- easyDR::alcohol_cancer
  ari <- easyDR::aripiprazole

  # A counts-only and an effect-only view of the same table, so every intake
  # branch of the preparation code is walked.
  crc_counts <- crc[, c("studlab", "dose", "event.e", "n.e", "design", "ref")]
  crc_effect <- crc[, c("studlab", "dose", "TE", "lower", "upper", "design",
                        "ref")]

  # A study-level covariate for the continuous side, which the bundled trial
  # table does not carry.
  ari_mod <- ari
  u <- unique(ari_mod$studlab)
  ari_mod$arm <- ifelse(ari_mod$studlab %in% u[seq_len(2)], "fixed",
                        "flexible")

  list(
    list(label = "bin linear, effect and counts",
         f = function() easyDR:::easydrbin(crc, sm = "IRR",
                                          dose.unit = "g/day",
                                          at = c(10, 25, 50))),
    list(label = "bin linear, default doses",
         f = function() easyDR:::easydrbin(crc, sm = "IRR",
                                          dose.unit = "g/day")),
    list(label = "bin linear, counts intake only",
         f = function() easyDR:::easydrbin(crc_counts, sm = "IRR",
                                          dose.unit = "g/day")),
    list(label = "bin linear, effect intake only",
         f = function() easyDR:::easydrbin(crc_effect, sm = "IRR",
                                          dose.unit = "g/day")),
    list(label = "bin linear, one stage ML",
         f = function() easyDR:::easydrbin(crc, sm = "IRR", proc = "1stage",
                                          method.tau = "ML")),
    list(label = "bin rcs",
         f = function() easyDR:::easydrbin(crc, sm = "IRR",
                                          dose.unit = "g/day",
                                          shape = "rcs")),
    list(label = "bin rcs, stated knots",
         f = function() easyDR:::easydrbin(crc, sm = "IRR", shape = "rcs",
                                          knots = c(5, 15, 40))),
    list(label = "bin fp",
         f = function() easyDR:::easydrbin(crc, sm = "IRR", shape = "fp",
                                          powers = c(-1, 0))),
    list(label = "bin quadratic",
         f = function() easyDR:::easydrbin(crc, sm = "IRR",
                                          shape = "quadratic",
                                          at = c(10, 25))),
    list(label = "bin log",
         f = function() easyDR:::easydrbin(crc, sm = "IRR", shape = "log")),
    list(label = "bin curve grid with extrapolation",
         f = function() easyDR:::easydrbin(crc, sm = "IRR", shape = "rcs",
                                          grid = c(0, 20, 40, 200))),
    list(label = "bin RR on the mixed-design table",
         f = function() easyDR:::easydrbin(cvd, sm = "RR",
                                          dose.unit = "g/day",
                                          at = c(10, 25, 50))),
    list(label = "bin OR from counts",
         f = function() {
           cc <- crc_counts
           cc$design <- "case-control"
           easyDR:::easydrbin(cc, sm = "OR")
         }),
    list(label = "cont MD",
         f = function() easyDR:::easydrcont(ari, sm = "MD",
                                           dose.unit = "mg/day",
                                           at = c(5, 10, 20))),
    list(label = "cont SMD default doses",
         f = function() easyDR:::easydrcont(ari, sm = "SMD",
                                           dose.unit = "mg/day")),
    list(label = "cont SMD hedges",
         f = function() easyDR:::easydrcont(ari, sm = "SMD",
                                           method.smd = "hedges")),
    list(label = "cont quadratic",
         f = function() easyDR:::easydrcont(ari, sm = "SMD",
                                           shape = "quadratic",
                                           at = c(5, 10, 20))),
    list(label = "cont rcs",
         f = function() easyDR:::easydrcont(ari, sm = "SMD", shape = "rcs")),
    list(label = "loo on a line",
         f = function() easyDR:::easydr_loo(
           easyDR:::easydrbin(cvd, sm = "RR", dose.unit = "g/day"),
           quiet = TRUE)),
    list(label = "loo on a curve at a named dose",
         f = function() easyDR:::easydr_loo(
           easyDR:::easydrbin(cvd, sm = "RR", shape = "quadratic",
                             at = c(10, 25)),
           at = 25, quiet = TRUE)),
    list(label = "loo with failing refits",
         f = function() {
           three <- cvd[cvd$studlab %in% c("Bianchi", "Bobak",
                                           "Vliegenthart-wine"), ]
           easyDR:::easydr_loo(easyDR:::easydrbin(three, sm = "RR",
                                            shape = "quadratic"),
                          quiet = TRUE)
         }),
    list(label = "cumulative, default order",
         f = function() easyDR:::dr_cumul(
           easyDR:::easydrbin(cvd, sm = "RR", dose.unit = "g/day"),
           quiet = TRUE)),
    list(label = "cumulative, stated order",
         f = function() easyDR:::dr_cumul(
           easyDR:::easydrbin(cvd, sm = "RR", dose.unit = "g/day"),
           order = rev(unique(cvd$studlab)), quiet = TRUE)),
    list(label = "loo on continuous",
         f = function() easyDR:::easydr_loo(
           easyDR:::easydrcont(ari, sm = "SMD", dose.unit = "mg/day"),
           quiet = TRUE)),
    list(label = "moderator, categorical",
         f = function() easyDR:::easydr_metareg(
           easyDR:::easydrbin(cvd, sm = "RR", dose.unit = "g/day"),
           ~ design)),
    list(label = "moderator, continuous and centred",
         f = function() {
           d <- crc
           d$year <- c(1990, 1995, 1998, 2001, 2004, 2007, 2010, 2013)[
             match(d$studlab, unique(d$studlab))]
           easyDR:::easydr_metareg(easyDR:::easydrbin(d, sm = "IRR"), ~ year)
         }),
    list(label = "moderator, two covariates",
         f = function() {
           d <- crc
           uu <- unique(d$studlab)
           d$region <- ifelse(d$studlab %in% uu[seq_len(4)], "europe",
                              "america")
           d$year <- c(1990, 1995, 1998, 2001, 2004, 2007, 2010, 2013)[
             match(d$studlab, uu)]
           easyDR:::easydr_metareg(easyDR:::easydrbin(d, sm = "IRR"),
                                ~ region + year)
         }),
    list(label = "moderator on a curve",
         f = function() easyDR:::easydr_metareg(
           easyDR:::easydrbin(cvd, sm = "RR", shape = "quadratic"), ~ design)),
    list(label = "moderator on continuous",
         f = function() easyDR:::easydr_metareg(
           easyDR:::easydrcont(ari_mod, sm = "SMD"), ~ arm)),
    list(label = "summary of a fit",
         f = function() summary(easyDR:::easydrbin(crc, sm = "IRR",
                                                  dose.unit = "g/day"))),
    list(label = "summary of a curve",
         f = function() summary(easyDR:::easydrbin(crc, sm = "IRR",
                                                  shape = "rcs")))
  )
}

# Every message easydr refuses with. A refactor that merged two validation
# branches would most likely show up here first, as a changed or vanished
# message, and these messages are part of what the package promises.
dr_validation_errors <- function() {
  crc <- easyDR::alcohol_cancer
  cvd <- easyDR::alcohol_cvd
  counts <- crc[, c("studlab", "dose", "event.e", "n.e", "design", "ref")]

  list(
    list(label = "V1 non-numeric dose", f = function() {
      d <- crc
      d$dose <- as.character(d$dose)
      easyDR:::easydrbin(d, sm = "IRR")
    }),
    list(label = "V2 no reference row", f = function() {
      d <- crc
      d$ref <- FALSE
      easyDR:::easydrbin(d, sm = "IRR")
    }),
    list(label = "V2 two reference rows", f = function() {
      d <- crc
      d$ref <- TRUE
      easyDR:::easydrbin(d, sm = "IRR")
    }),
    list(label = "V13 OR from cumulative-incidence counts", f = function() {
      easyDR:::easydrbin(counts, sm = "OR")
    }),
    list(label = "V13 RR from case-control counts", f = function() {
      d <- counts
      d$design <- "case-control"
      easyDR:::easydrbin(d, sm = "RR")
    }),
    list(label = "V13 IRR without person-time", f = function() {
      d <- counts
      d$design <- "case-control"
      easyDR:::easydrbin(d, sm = "IRR")
    }),
    list(label = "V16 knots outside the dose range", f = function() {
      easyDR:::easydrbin(crc, sm = "IRR", shape = "rcs",
                        knots = c(200, 300, 400))
    }),
    list(label = "V14 shape the data cannot inform", f = function() {
      two <- crc[crc$studlab %in% c("atm", "hpm"), ]
      easyDR:::easydrbin(two, sm = "IRR", shape = "rcs",
                        knots = c(5, 10, 25, 50))
    }),
    list(label = "unknown shape", f = function() {
      easyDR:::easydrbin(crc, sm = "IRR", shape = "wiggly")
    }),
    list(label = "loo needs two studies", f = function() {
      one <- cvd[cvd$studlab == "Bianchi", ]
      easyDR:::easydr_loo(easyDR:::easydrbin(one, sm = "RR"), quiet = TRUE)
    }),
    list(label = "V22 cumulative order mismatch", f = function() {
      easyDR:::dr_cumul(easyDR:::easydrbin(cvd, sm = "RR"),
                       order = c("Bianchi", "Hokusai"))
    }),
    list(label = "V22 cumulative order repeats", f = function() {
      easyDR:::dr_cumul(easyDR:::easydrbin(cvd, sm = "RR"),
                       order = c(unique(cvd$studlab), "Bianchi"))
    }),
    list(label = "V19 moderator varies inside a study", f = function() {
      d <- cvd
      d$site <- seq_len(nrow(d))
      easyDR:::easydr_metareg(easyDR:::easydrbin(cvd, sm = "RR"), ~ site,
                           data = d)
    }),
    list(label = "V20 moderator constant across studies", f = function() {
      d <- cvd
      d$one <- "same"
      easyDR:::easydr_metareg(easyDR:::easydrbin(cvd, sm = "RR"), ~ one, data = d)
    }),
    list(label = "V21 too few studies for the moderator model",
         f = function() {
           four <- cvd[cvd$studlab %in% c("Bianchi", "Bobak",
                                          "Vliegenthart-wine",
                                          "Vliegenthart-beer"), ]
           easyDR:::easydr_metareg(
             easyDR:::easydrbin(four, sm = "RR", shape = "quadratic"),
             ~ design)
         }),
    list(label = "moderator column absent", f = function() {
      easyDR:::easydr_metareg(easyDR:::easydrbin(cvd, sm = "RR"), ~ nope)
    }),
    list(label = "moderator formula two-sided", f = function() {
      easyDR:::easydr_metareg(easyDR:::easydrbin(cvd, sm = "RR"), y ~ design)
    }),
    list(label = "refit may not override the model", f = function() {
      fit <- easyDR:::easydrbin(cvd, sm = "RR")
      inh <- easyDR:::dr_inherited(fit)
      easyDR:::dr_refit(inh, inh$studies, list(sm = "OR"))
    })
  )
}

# Printed output carries an elapsed time for the analyses that loop, and that
# number is different on every run. It is the one thing in the output that is
# allowed to vary, so it is replaced by a placeholder before comparison. Doing
# it here, in one function both sides call, keeps the builder and the test
# honest with each other.
dr_validation_print <- function(x) {
  out <- utils::capture.output(print(x))
  out <- sub("[0-9.]+ seconds", "<elapsed> seconds", out)
  sub("^[0-9.]+ seconds", "<elapsed> seconds", out)
}

# Pull every number out of a result, whatever kind of object it is. Printed
# output rounds; this does not.
dr_validation_numbers <- function(x) {
  keep <- function(v) if (is.null(v)) NULL else v
  num <- function(d) {
    # Not every slot is a data frame: an empty tie group comes back as a bare
    # vector, and a slot can be absent entirely.
    if (is.null(d)) return(NULL)
    if (!is.data.frame(d)) return(if (is.numeric(d)) d else NULL)
    if (nrow(d) == 0L || ncol(d) == 0L) return(NULL)
    d <- d[, vapply(d, is.numeric, logical(1)), drop = FALSE]
    if (ncol(d) == 0L) return(NULL)
    as.matrix(d)
  }

  if (inherits(x, "easydr_result")) {
    return(list(
      slope = keep(unlist(x$slope[c("estimate", "se", "lb", "ub", "p")])),
      het = keep(unlist(x$heterogeneity[c("Q", "df", "p", "i2")])),
      psi = keep(as.numeric(x$Psi)),
      coef = stats::coef(x$fit),
      vcov = as.numeric(stats::vcov(x$fit)),
      predicted = num(x$predicted),
      curve = num(x$curve),
      nonlin = keep(unlist(x$nonlinearity[c("chi2", "df", "p")])),
      gof = keep(unlist(x$gof[c("D", "df", "p", "r2")])),
      ref = x$reference$curve,
      counts = unlist(x$counts)
    ))
  }
  if (inherits(x, "easydr_sensitivity")) {
    return(list(table = num(x$table), full = unlist(x$full[c("raw", "est")]),
                dose = x$dose, tried = x$tried, fitted = x$fitted,
                failed = x$failed))
  }
  if (inherits(x, "easydr_moderator")) {
    return(list(terms = num(x$terms),
                test = keep(unlist(x$test[c("chi2", "df", "p")])),
                het = keep(unlist(x$heterogeneity[c("Q", "df", "p", "i2")])),
                before = keep(unlist(
                  x$heterogeneity_before[c("Q", "df", "p", "i2")])),
                centres = unlist(x$centres),
                coef = stats::coef(x$fit)))
  }
  if (inherits(x, "summary.easydr_result")) {
    return(list(per_study = num(x$per_study)))
  }
  list(unclassed = NULL)
}
