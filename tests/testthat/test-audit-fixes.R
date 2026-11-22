# Regression tests for the 2026-11-22 calculation audit, one block per
# finding (F1 to F18). Each reuses the table the audit reproduced the finding
# on, and each checks the fixed path against the engine called directly where
# the engine can be called directly.

# alcohol_cvd with every reference dose except Bianchi's moved to 5, so the
# full fit's reference (0) and a subset's own lowest reference (5) differ.
cvd_shifted_refs <- function() {
  d <- alcohol_cvd
  for (s in setdiff(unique(d$studlab), "Bianchi")) {
    d$dose[d$studlab == s & d$ref] <- 5
  }
  d
}

# --- F1: moderator columns in a one-stage fit --------------------------------

test_that("F1: a one-stage moderator test covers the moderator columns only", {
  f1 <- easydr(alcohol_cvd, sm = "RR", shape = "linear", proc = "1stage")
  f2 <- easydr(alcohol_cvd, sm = "RR", shape = "linear")
  m1 <- suppressWarnings(easydr_metareg(f1, ~ design))
  m2 <- easydr_metareg(f2, ~ design)

  # The engine names them differently; easydr must read both the same way.
  expect_identical(names(coef(m1$fit)),
                   c("dose", "dose:designcumulative incidence"))
  expect_identical(m1$terms$moderator, c(FALSE, TRUE))
  expect_identical(m2$terms$moderator, c(FALSE, TRUE))
  expect_identical(m1$test$df, 1L)
  expect_identical(m2$test$df, 1L)
  expect_identical(unname(m1$test$terms), 2L)

  # The one-stage test is dosresmeta's Wald test on that one column.
  w <- dosresmeta::waldtest(b = coef(m1$fit), Sigma = vcov(m1$fit), Terms = 2)
  expect_equal(m1$test$chi2, as.numeric(w$chitest["chi2"]))

  # And the labels say what the one-stage coefficients are.
  expect_match(m1$terms$label[1], "slope at design = case-control")
  expect_match(m1$terms$label[2], "change in slope, cumulative incidence")
})

test_that("F1: a one-stage multi-term moderator fit is read term by term", {
  f <- suppressWarnings(easydr(alcohol_cvd, sm = "RR", shape = "quadratic",
                               proc = "1stage"))
  m <- suppressWarnings(easydr_metareg(f, ~ design))
  expect_identical(m$terms$moderator, c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(m$test$df, 2L)
  expect_match(m$terms$label[1], "^the dose term at design")
  expect_match(m$terms$label[4], "^change in the dose_sq term")
})

test_that("F1: the one-stage subgroup test has the moderator's degrees of freedom", {
  f <- suppressWarnings(easydr(alcohol_cvd, sm = "RR", shape = "linear",
                               proc = "1stage", subgroup = "design"))
  expect_true(f$subgroup$test$estimable)
  expect_identical(f$subgroup$test$df, 1L)
})

# --- F2: profiles run under one-stage ML -------------------------------------

test_that("F2: fractional polynomial powers do not depend on proc or method.tau", {
  a <- easydr(alcohol_cvd, sm = "RR", shape = "fp")
  b <- suppressWarnings(easydr(alcohol_cvd, sm = "RR", shape = "fp",
                               proc = "1stage", method.tau = "ML"))
  expect_identical(a$shape$powers, b$shape$powers)
  expect_identical(a$provenance$proc, "2stage")
  expect_identical(a$provenance$method, "reml")

  # The chosen pair is the one-stage ML winner over the whole grid, each
  # candidate fitted by dosresmeta directly on easydr's own basis.
  frame <- a$frame
  set   <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  pairs <- expand.grid(p1 = set, p2 = set)
  pairs <- pairs[pairs$p1 <= pairs$p2, ]
  ll <- apply(pairs, 1, function(p) {
    shp <- easyDR:::dr_shape_fp(frame$dose, as.numeric(p))
    dat <- data.frame(y = frame$y, id = frame$id, v = frame$v,
                      as.data.frame(shp$basis(frame$dose)))
    fit <- try(suppressWarnings(dosresmeta::dosresmeta(
      y ~ fp1 + fp2, id = id, v = v, type = frame$type, cases = frame$cases,
      n = frame$n, data = dat, covariance = "gl", method = "ml",
      proc = "1stage", control = list(maxiter = 2000))), silent = TRUE)
    if (inherits(fit, "try-error") || isFALSE(fit$converged)) return(-Inf)
    as.numeric(logLik(fit))
  })
  best <- as.numeric(pairs[which.max(ll), ])
  expect_equal(a$shape$powers, best)
})

test_that("F2: the Emax profile does not depend on proc or method.tau", {
  a <- suppressWarnings(easydr(alcohol_cvd, sm = "RR", shape = "emax"))
  b <- suppressWarnings(easydr(alcohol_cvd, sm = "RR", shape = "emax",
                               proc = "1stage", method.tau = "ML"))
  expect_equal(a$shape$ed50, b$shape$ed50)
  expect_equal(a$shape$steepness, b$shape$steepness)
  # The winner is then fitted under the caller's settings.
  expect_identical(a$provenance$proc, "2stage")
  expect_identical(a$provenance$method, "reml")
  expect_identical(b$provenance$proc, "1stage")
})

# --- F3: the marked reference of a continuous outcome ------------------------

test_that("F3: a non-lowest reference row is honoured for continuous outcomes", {
  d <- aripiprazole
  s1 <- which(d$studlab == unique(d$studlab)[1])
  d$ref[s1] <- FALSE
  d$ref[s1[3]] <- TRUE
  moved <- suppressWarnings(easydr(d, sm = "MD", shape = "linear"))
  plain <- easydr(aripiprazole, sm = "MD", shape = "linear")

  # The engine takes the first row of each study as the reference, so the
  # reference is the row easydr put first.
  expect_equal(unname(moved$reference$per_study[1]), d$dose[s1[3]])
  expect_true(moved$frame$is_ref[moved$frame$studlab == d$studlab[s1[1]]][1])

  # The fit is the one dosresmeta gives with that row first in the study.
  ord <- order(moved$frame$id, !moved$frame$is_ref, moved$frame$dose)
  dd  <- moved$frame[ord, ]
  ref <- dosresmeta::dosresmeta(y ~ dose, id = id, sd = cases, n = n,
                                data = dd, covariance = "md", method = "reml")
  expect_equal(moved$slope$estimate, as.numeric(coef(ref)), tolerance = 1e-10)
  expect_equal(moved$slope$se, sqrt(as.numeric(vcov(ref))), tolerance = 1e-10)
  expect_equal(moved$fit$Slist[[1]], ref$Slist[[1]])

  # Which is not the fit against the lowest dose.
  expect_false(isTRUE(all.equal(moved$slope$se, plain$slope$se)))
  expect_false(isTRUE(all.equal(moved$fit$Slist[[1]], plain$fit$Slist[[1]])))
})

# --- F4: a risk ratio from person-time counts --------------------------------

test_that("F4: RR from incidence-rate counts is refused, IRR is accepted", {
  tab <- crc_table("counts")
  expect_error(easydrbin(tab, sm = "RR"), "person-time denominator gives a rate ratio")
  expect_s3_class(easydrbin(tab, sm = "IRR"), "easydr_result")
  # With the reported effect, RR only names the scale and is allowed.
  expect_s3_class(easydrbin(crc_table("effect_counts"), sm = "RR"),
                  "easydr_result")
})

# --- F5: one rule for the spline knots ---------------------------------------

test_that("F5: every spelling of three knots gives the same percentiles", {
  d <- alcohol_cvd$dose
  want <- as.numeric(stats::quantile(d, c(0.10, 0.50, 0.90), names = FALSE))
  expect_equal(easyDR:::dr_shape("rcs", d)$knots, want)
  expect_equal(easyDR:::dr_shape("rcs", d, knots = 3)$knots, want)
  expect_equal(easyDR:::dr_shape("rcs", d, knots = c(.1, .5, .9))$knots, want)
  # The reference zero is part of the distribution, as in dosresmeta.
  expect_equal(easyDR:::dr_shape("rcs", aripiprazole$dose)$knots,
               c(0, 10, 30))
  # Harrell's positions for four and five knots.
  expect_equal(easyDR:::dr_shape("rcs", d, knots = 4)$knots,
               as.numeric(stats::quantile(d, c(.05, .35, .65, .95),
                                          names = FALSE)))
  expect_equal(easyDR:::dr_shape("rcs", d, knots = 5)$knots,
               as.numeric(stats::quantile(d, c(.05, .275, .5, .725, .95),
                                          names = FALSE)))
  # And the fitted curve is the same whichever way it was asked for.
  a <- easydr(alcohol_cvd, sm = "RR", shape = "rcs")
  b <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", knots = 3)
  c3 <- easydr(alcohol_cvd, sm = "RR", shape = "rcs", knots = c(.1, .5, .9))
  expect_equal(a$predicted$est, b$predicted$est)
  expect_equal(a$predicted$est, c3$predicted$est)
})

# --- F6: the bootstrap keeps what the caller fixed ---------------------------

test_that("F6: a bootstrap draw holds a user-fixed ED50 and searches the rest", {
  f <- suppressWarnings(easydr(alcohol_cvd, sm = "RR", shape = "emax",
                               ed50 = 20))
  expect_equal(f$shape$ed50, 20)
  expect_equal(f$args$ed50, 20)
  expect_equal(f$shape$profiled_params, 1L)
  # A draw the engine cannot fit is NULL by design; the ones that fit must
  # keep ED50 where the caller put it.
  set.seed(1)
  draws <- Filter(Negate(is.null),
                  lapply(1:6, function(i) easyDR:::dr_boot_once(f, TRUE)))
  expect_gte(length(draws), 3L)
  for (d in draws) {
    expect_equal(d$shape$ed50, 20)
    expect_true(d$shape$steepness %in% c(0.5, 0.75, 1, 1.5, 2, 3))
  }
})

# --- F7: refits keep the full fit's reference dose ---------------------------

test_that("F7: leave-one-out rows are contrasts against the full fit's reference", {
  d <- cvd_shifted_refs()
  f <- suppressWarnings(easydr(d, sm = "RR", shape = "quadratic"))
  expect_equal(f$reference$curve, 0)
  l <- suppressWarnings(easydr_loo(f, plot = FALSE, quiet = TRUE))
  for (nm in names(l$refits)) {
    r <- l$refits[[nm]]
    if (is.null(r)) next
    expect_equal(r$reference$curve, 0, info = nm)
    row <- l$table[l$table$label == nm, ]
    # The row is RR(1 vs 0) from the refit's own coefficients.
    want <- as.numeric(easyDR:::dr_contrast(r$shape, 1, 0) %*% coef(r$fit))
    expect_equal(row$estimate_raw, want, info = nm)
  }
  # Omitting Bianchi leaves no study with reference 0; the refit says so.
  r <- l$refits[["Bianchi"]]
  expect_true(any(grepl("inherited reference 0", r$provenance$warnings)))
})

test_that("F7: subgroup levels are contrasts against the overall reference", {
  d <- cvd_shifted_refs()
  f <- suppressWarnings(easydr(d, sm = "RR", shape = "linear",
                               subgroup = "design"))
  for (lv in names(f$subgroup$fits)) {
    r <- f$subgroup$fits[[lv]]
    if (is.null(r)) next
    expect_equal(r$reference$curve, f$reference$curve, info = lv)
  }
})

test_that("F7: a refit may not be handed a different reference", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "linear")
  inh <- easyDR:::dr_inherited(f)
  expect_error(easyDR:::dr_refit(inh, unique(alcohol_cvd$studlab)[-1],
                                 extra = list(reference = 5)),
               "may not override")
})

# --- F8: one definition of a target dose -------------------------------------

test_that("F8: target doses are read inside the range for every shape", {
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax", ed50 = 5,
              dose.unit = "mg/day")
  got <- ed_gamma(f, gamma = c(0.5, 0.8))
  # E(d) is proportional to d / (d + 5). Inside 0 to 30 the largest effect
  # is at 30, so ED(g) solves d / (d + 5) = g * 30 / 35.
  top  <- 30 / 35
  want <- 5 * (c(0.5, 0.8) * top) / (1 - c(0.5, 0.8) * top)
  expect_equal(unname(got$ed), want, tolerance = 1e-3)
  expect_equal(got$ed_max, 30)
  # The Emax parameter is kept apart and labelled as the parameter.
  expect_equal(got$asymptotic$ed50, 5)
  expect_equal(unname(got$asymptotic$ed), 5 * c(0.5, 0.8) / (1 - c(0.5, 0.8)))
  td <- suppressWarnings(easydr_target_dose(f, B = 50, seed = 1))
  expect_equal(td$table$Dose[td$table$Quantity == "ED50"], unname(got$ed[1]))
  expect_output(print(td), "Emax parameter ED50 (asymptotic", fixed = TRUE)
  expect_output(print(f), "Emax parameter ED50 (asymptotic): 5.00",
                fixed = TRUE)
  # A spline carries no asymptotic parameter.
  s <- easydr(alcohol_cancer, sm = "RR", shape = "rcs")
  expect_null(ed_gamma(s)$asymptotic)
})

# --- F9: counts on the effect path are checked -------------------------------

test_that("F9: a zero cell beside a reported effect is corrected, not crashed", {
  d <- alcohol_cvd
  i <- which(d$studlab == "Bianchi" & !d$ref)[1]
  d$event.e[i] <- 0
  expect_warning(f <- easydr(d, sm = "RR", shape = "linear"), "zero cases")
  expect_true(all(is.finite(f$fit$Slist[[1]])))
  expect_true("Bianchi" %in% f$provenance$continuity)
  # The reported effects are untouched; only the covariance counts move.
  o <- d[d$studlab == "Bianchi", ]
  o <- o[order(o$dose), ]
  b <- f$frame[f$frame$studlab == "Bianchi", ]
  expect_equal(b$y, log(o$TE))
  expect_equal(b$cases, o$event.e + 0.5)
  expect_equal(b$n, o$n.e + 1)
})

test_that("F9: impossible counts beside a reported effect stop with a reason", {
  d <- alcohol_cvd
  i <- which(d$studlab == "Bianchi" & !d$ref)[1]
  d$event.e[i] <- d$n.e[i] + 10
  expect_error(easydr(d, sm = "RR", shape = "linear"), "impossible counts")
  d <- alcohol_cvd
  d$event.e[i] <- -1
  expect_error(easydr(d, sm = "RR", shape = "linear"), "impossible counts")
})

test_that("F9: a study with missing counts is said to have independent contrasts", {
  d <- alcohol_cvd
  i <- which(d$studlab == "Bianchi" & !d$ref)[1]
  d$event.e[i] <- NA
  expect_warning(f <- easydr(d, sm = "RR", shape = "linear"),
                 "missing counts, so its contrasts are treated as independent")
  expect_match(f$provenance$covariance_label, "independent contrasts for the 1 study")
  expect_true(any(grepl("missing counts", f$provenance$warnings)))
})

# --- F10: no silent REML fallback for AIC and BIC ----------------------------

test_that("F10: AIC and BIC are NA when the ML refit is missing", {
  f <- easydr(alcohol_cvd, sm = "RR", shape = "linear")
  expect_identical(easyDR:::dr_ic(f, "AIC", NULL), NA_real_)
  expect_identical(easyDR:::dr_ic(f, "BIC", NULL), NA_real_)
  ml <- easyDR:::dr_refit_ml(f)
  expect_equal(easyDR:::dr_ic(f, "AIC", ml), stats::AIC(ml$fit))
  # A result the refit cannot be built from gives NA, not the REML value.
  g <- f
  g$data <- NULL
  expect_null(easyDR:::dr_refit_ml(g))
  expect_identical(easyDR:::dr_ic(g, "AIC", easyDR:::dr_refit_ml(g)), NA_real_)
})

# --- F11: the zero-cell correction is one set of counts ----------------------

test_that("F11: counts path corrects effects and covariance from the same counts", {
  d <- alcohol_cvd
  d$TE <- NULL; d$lower <- NULL; d$upper <- NULL
  d$design <- "cohort"
  i <- which(d$studlab == "Bianchi" & !d$ref)[2]
  d$event.e[i] <- 0
  expect_warning(f <- easydr(d, sm = "RR", shape = "linear"),
                 "add = 1/2, to = all")
  b  <- f$frame[f$frame$studlab == "Bianchi", ]
  o  <- d[d$studlab == "Bianchi", ]
  o  <- o[order(o$dose), ]
  ir <- which(o$ref)
  k  <- nrow(o)
  es <- metafor::escalc("RR", ai = o$event.e, bi = o$n.e - o$event.e,
                        ci = rep(o$event.e[ir], k),
                        di = rep(o$n.e[ir] - o$event.e[ir], k),
                        add = 1 / 2, to = "all")
  expect_equal(b$y[-ir], as.numeric(es$yi)[-ir])
  expect_equal(b$v[-ir], as.numeric(es$vi)[-ir])
  expect_equal(b$cases, o$event.e + 0.5)
  expect_equal(b$n, o$n.e + 1)
})

# --- F12: BIC charges profiled parameters with the engine's n ----------------

test_that("F12: the BIC add-on uses the number of non-reference contrasts", {
  f  <- suppressWarnings(easydr(alcohol_cvd, sm = "RR", shape = "hyperbolic_emax"))
  ml <- easyDR:::dr_refit_ml(f)
  nobs <- attr(stats::logLik(ml$fit), "nobs")
  expect_equal(nobs, sum(!f$frame$is_ref))
  expect_equal(easyDR:::dr_ic(f, "BIC", ml),
               stats::BIC(ml$fit) + log(nobs) * 1)
  expect_false(isTRUE(all.equal(easyDR:::dr_ic(f, "BIC", ml),
                                stats::BIC(ml$fit) + log(f$counts$categories))))
})

# --- F13: a seeded call leaves the caller's stream alone ---------------------

test_that("F13: a seeded bootstrap restores the global RNG state", {
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax", ed50 = 5)
  set.seed(99); r_before <- runif(1)
  set.seed(99)
  invisible(suppressWarnings(easydr_target_dose(f, B = 50, seed = 7)))
  expect_equal(runif(1), r_before)
  set.seed(99)
  invisible(easyDR:::dr_boot_curve(f, 10, B = 50, reprofile = FALSE, seed = 7))
  expect_equal(runif(1), r_before)
  td <- suppressWarnings(easydr_target_dose(f, B = 50, seed = 7))
  expect_equal(td$seed, 7)
  expect_output(print(td), "seed 7")
})

# --- F14: only the flat curve is caught --------------------------------------

test_that("F14: a flat curve raises its own condition class", {
  f <- easydr(aripiprazole, sm = "MD", shape = "hyperbolic_emax", ed50 = 5)
  f$fit$coefficients[] <- 0
  expect_error(ed_gamma(f), class = "easydr_flat_curve")
  expect_error(ed_gamma(f), "flat fit has no target dose")
})

# --- F15: the summary names the first coefficient of a curve -----------------

test_that("F15: summary labels a slope as a slope and a curve's first term as such", {
  lin <- paste(capture.output(print(summary(
    easydr(alcohol_cancer, sm = "RR", shape = "linear")))), collapse = "\n")
  expect_match(lin, "slope (log scale)", fixed = TRUE)
  rcs <- paste(capture.output(print(summary(
    easydr(alcohol_cancer, sm = "RR", shape = "rcs")))), collapse = "\n")
  expect_match(rcs, "first coefficient (dose) (log scale)", fixed = TRUE)
  expect_false(grepl("slope (log scale)", rcs, fixed = TRUE))
  md <- paste(capture.output(print(summary(
    easydr(aripiprazole, sm = "MD", shape = "linear")))), collapse = "\n")
  expect_match(md, " slope ")
  expect_false(grepl("log scale", md, fixed = TRUE))
})

# --- F16: inferring the reference --------------------------------------------

test_that("F16: a second effect of 1 with limits does not take the reference", {
  d <- alcohol_cvd
  d$ref <- NULL
  j <- which(d$studlab == "Bianchi" & d$dose > 0)[1]
  d$TE[j] <- 1
  d$lower[j] <- 0.8; d$upper[j] <- 1.25
  f <- easydr(d, sm = "RR", shape = "linear")
  expect_equal(unname(f$reference$per_study["Bianchi"]), 0)
  # Two candidates without limits cannot be told apart: the table must say.
  d$lower[j] <- NA; d$upper[j] <- NA
  expect_error(easydr(d, sm = "RR", shape = "linear"),
               "reference category cannot be inferred")
})

# --- F17: subgroup order is locale independent -------------------------------

test_that("F17: subgroup levels are ordered in the C locale", {
  d <- alcohol_cvd
  lv <- c(Bianchi = "b", Bobak = "A", `Malarcher-wine` = "a",
          `Malarcher-beer` = "b", `Vliegenthart-wine` = "A",
          `Vliegenthart-beer` = "a")
  d$grp <- unname(lv[d$studlab])
  f <- suppressWarnings(easydr(d, sm = "RR", shape = "linear", subgroup = "grp"))
  expect_identical(f$subgroup$levels, c("A", "a", "b"))
})

# --- F18: covariance names match the outcome ---------------------------------

test_that("F18: a covariance that belongs to the other outcome is refused by name", {
  expect_error(easydrbin(alcohol_cvd, sm = "RR", shape = "linear",
                         covariance = "md"),
               "covariance must be one of")
  expect_error(easydr(aripiprazole, sm = "MD", shape = "linear",
                      covariance = "gl"),
               "not read")
  expect_s3_class(easydrbin(alcohol_cvd, sm = "RR", shape = "linear",
                            covariance = "h"), "easydr_result")
})
