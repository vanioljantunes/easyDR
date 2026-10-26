# Table intake: column roles, validation rules V1 to V13, and the analysis frame
# handed to the engine. Nothing here is exported and the engine is never called
# from this file.

# Standard column names by role. Auto-detection looks for these and nothing else.
easydr_names <- list(
  studlab = "studlab",
  dose    = "dose",
  ref     = "ref",
  design  = "design",
  event.e = "event.e",
  n.e     = "n.e",
  TE      = "TE",
  lower   = "lower",
  upper   = "upper",
  mean.e  = "mean.e",
  sd.e    = "sd.e"
)

# Engine design codes, keyed by the words a reviewer actually writes.
easydr_designs <- c(
  "cc" = "cc", "case-control" = "cc", "case control" = "cc", "casecontrol" = "cc",
  "ir" = "ir", "incidence rate" = "ir", "incidence-rate" = "ir", "rate" = "ir",
  "person-time" = "ir",
  "ci" = "ci", "cumulative incidence" = "ci", "cumulative-incidence" = "ci",
  "cohort" = "ci"
)

#' Resolve one column role to a column name
#' @keywords internal
#' @noRd
dr_role <- function(data, supplied, role, required = FALSE) {
  if (!is.null(supplied)) {
    if (!is.character(supplied) || length(supplied) != 1L) {
      stop("easydr: `", role, "` must be a single column name.", call. = FALSE)
    }
    if (!supplied %in% names(data)) {
      stop("easydr: column ", dQuote(supplied), " given for `", role,
           "` is not in the data. Columns present: ",
           paste(names(data), collapse = ", "), ".", call. = FALSE)
    }
    return(supplied)
  }
  standard <- easydr_names[[role]]
  if (!is.null(standard) && standard %in% names(data)) return(standard)
  if (required) {
    # Rule V12: say what was looked for and what is present.
    stop("easydr: cannot find a column for `", role, "`. Looked for ",
         dQuote(standard), ". Columns present: ",
         paste(names(data), collapse = ", "), ". Name it explicitly, for ",
         "example ", role, " = ", dQuote("my_column"), ".", call. = FALSE)
  }
  NULL
}

#' Resolve every column role for one entry point
#' @keywords internal
#' @noRd
dr_roles <- function(data, args, kind) {
  if (!is.data.frame(data)) {
    stop("easydr: `data` must be a data frame with one row per dose category.",
         call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("easydr: `data` has no rows.", call. = FALSE)
  }
  roles <- list(
    studlab = dr_role(data, args$studlab, "studlab", required = TRUE),
    dose    = dr_role(data, args$dose,    "dose",    required = TRUE),
    ref     = dr_role(data, args$ref,     "ref"),
    design  = dr_role(data, args$design,  "design")
  )
  if (identical(kind, "bin")) {
    roles$event.e <- dr_role(data, args$event.e, "event.e")
    roles$n.e     <- dr_role(data, args$n.e,     "n.e")
    roles$TE      <- dr_role(data, args$TE,      "TE")
    roles$lower   <- dr_role(data, args$lower,   "lower")
    roles$upper   <- dr_role(data, args$upper,   "upper")
    has_counts <- !is.null(roles$event.e) && !is.null(roles$n.e)
    has_effect <- !is.null(roles$TE) && !is.null(roles$lower) &&
      !is.null(roles$upper)
    if (!has_counts && !has_effect) {
      stop("easydr: no usable outcome columns. Supply counts (", dQuote("event.e"),
           " and ", dQuote("n.e"), ") or a reported effect (", dQuote("TE"), ", ",
           dQuote("lower"), " and ", dQuote("upper"), "). Columns present: ",
           paste(names(data), collapse = ", "), ".", call. = FALSE)
    }
    # A published table usually reports an adjusted effect AND the counts. The
    # effect is the estimate to pool; the counts are what the covariance
    # approximation needs. Recomputing a crude effect from the counts and
    # discarding the adjusted one would change the question being answered.
    roles$intake <- if (has_effect && has_counts) {
      "effect_counts"
    } else if (has_effect) {
      "effect"
    } else {
      "counts"
    }
  } else {
    roles$mean.e <- dr_role(data, args$mean.e, "mean.e", required = TRUE)
    roles$sd.e   <- dr_role(data, args$sd.e,   "sd.e",   required = TRUE)
    roles$n.e    <- dr_role(data, args$n.e,    "n.e",    required = TRUE)
    roles$intake <- "continuous"
  }
  roles
}

#' Map a design word to the engine code
#' @keywords internal
#' @noRd
dr_design_code <- function(x, study) {
  key  <- tolower(trimws(as.character(x)))
  code <- easydr_designs[key]
  if (length(code) == 0L || any(is.na(code))) {
    stop("easydr: study ", study, " has design ", dQuote(x),
         ", which is not one of: case-control, incidence rate, ",
         "cumulative incidence.", call. = FALSE)
  }
  unname(code[1])
}

#' Rule V13: the summary measure must suit the study design
#'
#' The rule depends on who computes the effect. Where easydr computes it from
#' counts, the measure has to be one those counts can yield: a case-control
#' study gives an odds ratio and not a risk ratio, and a rate ratio needs
#' person-time. Where the authors already reported the effect, `sm` names the
#' scale their ratios are pooled on, the design only selects the covariance
#' approximation, and designs may differ between studies. Pooling a mixture is
#' standard in this literature, so the rule allows it and the print names it.
#' A person-time denominator stays a fact about the data rather than a
#' reporting choice, so the rate-ratio check applies either way.
#'
#' @keywords internal
#' @noRd
dr_check_sm_design <- function(sm, design, study, intake = "counts") {
  if (identical(sm, "IRR") && !identical(design, "ir")) {
    stop("easydr: sm = ", dQuote("IRR"), " but study ", study, " is not marked ",
         "as incidence-rate data. A rate ratio needs a person-time ",
         "denominator. Set design = ", dQuote("incidence rate"),
         " for that study.", call. = FALSE)
  }
  if (!identical(intake, "counts")) {
    return(invisible(TRUE))
  }
  if (identical(sm, "OR") && !identical(design, "cc")) {
    stop("easydr: sm = ", dQuote("OR"), " but study ", study, " has a ", design,
         " design and easydr is computing the effect from its counts. An odds ",
         "ratio belongs to a case-control study. Set the design column, or ",
         "choose sm = ", dQuote("RR"), ".", call. = FALSE)
  }
  if (identical(sm, "RR") && identical(design, "cc")) {
    stop("easydr: sm = ", dQuote("RR"), " but study ", study, " has a ",
         "case-control design and easydr is computing the effect from its ",
         "counts, which can only give an odds ratio. Use sm = ", dQuote("OR"),
         ", or supply the reported effect with its limits.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Infer which row of a study is its reference category
#' @keywords internal
#' @noRd
dr_infer_ref <- function(data, roles, rows, d) {
  if (!is.null(roles$TE)) {
    te <- data[[roles$TE]][rows]
    at_one <- which(!is.na(te) & abs(te - 1) < 1e-8)
    if (length(at_one) == 1L) return(at_one)
    at_zero <- which(!is.na(te) & te == 0)
    if (length(at_zero) == 1L) return(at_zero)
  }
  which.min(d)
}

#' Resolve one study design code
#' @keywords internal
#' @noRd
dr_study_design <- function(data, roles, rows, study) {
  if (!is.null(roles$design)) {
    return(dr_design_code(data[[roles$design]][rows][1], study))
  }
  if (identical(roles$intake, "continuous")) return(NA_character_)
  # Without a design column, a table of events over totals is cumulative
  # incidence. The default is recorded in provenance and printed.
  "ci"
}

#' Validate the table and build the analysis frame
#'
#' Applies rules V1 to V13 from the data model, then returns the frame the engine
#' wants: one row per dose category with the reference row kept, carrying
#' `y = 0` and `v = NA`. The engine rejects a frame whose reference rows are
#' missing: dropping them gives "'x' must have positive length", because the
#' covariance approximation needs every category of a study.
#'
#' @keywords internal
#' @noRd
dr_prepare <- function(data, roles, sm) {
  prov <- list(intake = roles$intake, warnings = character(),
               inferred_ref = character(), continuity = character(),
               designs = character())

  studlab <- as.character(data[[roles$studlab]])
  dose    <- data[[roles$dose]]

  # V1: the dose must be a number the reviewer assigned.
  if (!is.numeric(dose)) {
    stop("easydr: the dose column ", dQuote(roles$dose), " is not numeric. ",
         "easydr needs the numeric dose assigned to each category, not a ",
         "category label such as low or high. Assign a dose, then call again.",
         call. = FALSE)
  }
  if (anyNA(dose)) {
    bad <- which(is.na(dose))
    stop("easydr: missing dose in rows ", paste(bad, collapse = ", "),
         " (studies ", paste(unique(studlab[bad]), collapse = ", "), ").",
         call. = FALSE)
  }

  # V3: a dose-response needs more than one dose.
  if (length(unique(dose)) < 2L) {
    stop("easydr: every row has dose ", unique(dose)[1],
         ". No dose-response is estimable from a single dose value.",
         call. = FALSE)
  }

  # V4: duplicate study and dose pairs. Two different studies sharing one label
  # are indistinguishable from one study and cannot be detected here.
  key <- paste(studlab, dose, sep = " @ ")
  if (anyDuplicated(key) > 0L) {
    dup <- unique(key[duplicated(key)])
    stop("easydr: duplicated study and dose pairs: ",
         paste(dup, collapse = "; "),
         ". Each row must be one dose category of one study.", call. = FALSE)
  }

  studies <- unique(studlab)
  pieces  <- vector("list", length(studies))

  for (i in seq_along(studies)) {
    s    <- studies[i]
    rows <- which(studlab == s)
    d    <- dose[rows]

    # V2: exactly one reference row, or infer one and say so.
    if (!is.null(roles$ref)) {
      flag <- as.logical(data[[roles$ref]][rows])
      flag[is.na(flag)] <- FALSE
      if (sum(flag) != 1L) {
        stop("easydr: study ", s, " has ", sum(flag), " reference rows marked ",
             "in ", dQuote(roles$ref), " (rows ", paste(rows, collapse = ", "),
             "). Exactly one row per study must be the reference category.",
             call. = FALSE)
      }
      iref <- which(flag)
    } else {
      iref <- dr_infer_ref(data, roles, rows, d)
      prov$inferred_ref <- c(prov$inferred_ref, s)
    }

    design <- dr_study_design(data, roles, rows, s)
    if (!is.null(sm) && !identical(roles$intake, "continuous")) {
      dr_check_sm_design(sm, design, s, roles$intake)
    }
    prov$designs <- c(prov$designs, stats::setNames(design, s))

    piece <- dr_study_frame(data, roles, rows, iref, s, sm, design, prov)
    prov  <- piece$prov
    pieces[[i]] <- piece$frame
  }

  frame <- do.call(rbind, pieces)
  frame$id <- match(frame$studlab, studies)
  frame <- frame[order(frame$id, frame$dose), , drop = FALSE]
  rownames(frame) <- NULL

  # V11: three dose levels are needed before a study can inform a nonlinearity
  # test, which a one-term linear fit cannot run in any case.
  levels_per_study  <- tapply(frame$dose, frame$studlab, function(x) length(unique(x)))
  prov$dose_levels  <- levels_per_study
  prov$n_with_three <- sum(levels_per_study >= 3L)

  # V10: a single study supports no heterogeneity and no prediction interval.
  prov$single_study <- length(studies) == 1L
  if (prov$single_study) {
    prov$warnings <- c(prov$warnings, paste0(
      "only one study supplied, so heterogeneity and the prediction interval ",
      "are not estimable"))
  }

  # V9: dose units are the reviewer's responsibility. Flag the obvious slip.
  maxima   <- tapply(frame$dose, frame$studlab, max)
  positive <- maxima[maxima > 0]
  if (length(positive) > 1L && max(positive) / min(positive) > 100) {
    hi  <- names(positive)[which.max(positive)]
    lo  <- names(positive)[which.min(positive)]
    msg <- paste0("per-study dose ranges differ by more than a factor of 100: ",
                  hi, " reaches ", signif(max(positive), 4), " while ", lo,
                  " reaches ", signif(min(positive), 4), ". Check that every ",
                  "study uses the same dose unit. easydr does not convert units")
    warning("easydr: ", msg, call. = FALSE)
    prov$warnings <- c(prov$warnings, msg)
  }

  list(frame = frame, provenance = prov, studies = studies)
}

#' Build the analysis frame for one study
#'
#' The reference row is kept with `y = 0` and `v = NA`, which is the shape the
#' engine's own datasets use.
#'
#' @keywords internal
#' @noRd
dr_study_frame <- function(data, roles, rows, iref, study, sm, design, prov) {
  d <- data[[roles$dose]][rows]
  n <- length(rows)
  y <- rep(NA_real_, n)
  v <- rep(NA_real_, n)
  cases <- rep(NA_real_, n)
  total <- rep(NA_real_, n)

  if (identical(roles$intake, "counts")) {
    cases <- as.numeric(data[[roles$event.e]][rows])
    total <- as.numeric(data[[roles$n.e]][rows])

    # V5: counts must be coherent.
    bad <- which(is.na(cases) | is.na(total) | cases < 0 | total <= 0 |
                   (identical(design, "ci") & cases > total))
    if (length(bad) > 0L) {
      stop("easydr: study ", study, " has impossible counts in rows ",
           paste(rows[bad], collapse = ", "),
           ". Cases must be present, non-negative, and no larger than the ",
           "denominator.", call. = FALSE)
    }

    # V6: a zero cell needs a correction, and the reviewer is told which. The
    # correction has to reach the counts handed to the engine as well, because
    # the covariance approximation divides by them and fails on a zero.
    es <- dr_escalc_counts(cases, total, iref, sm, design)
    y  <- es$y
    v  <- es$v
    if (any(cases == 0)) {
      msg <- paste0("study ", study, " has a dose category with zero cases, so ",
                    "the escalc correction add = 1/2, to = only0 was applied ",
                    "to that study, including the counts used for the ",
                    "covariance")
      warning("easydr: ", msg, call. = FALSE)
      prov$warnings   <- c(prov$warnings, msg)
      prov$continuity <- c(prov$continuity, study)
      cases <- cases + 0.5
    }
  } else if (roles$intake %in% c("effect", "effect_counts")) {
    te <- as.numeric(data[[roles$TE]][rows])
    lo <- as.numeric(data[[roles$lower]][rows])
    hi <- as.numeric(data[[roles$upper]][rows])

    # V8: limits present for every non-reference row, or for none.
    missing_limits <- setdiff(which(is.na(lo) | is.na(hi)), iref)
    if (length(missing_limits) > 0L && length(missing_limits) < n - 1L) {
      stop("easydr: study ", study, " reports confidence limits for some dose ",
           "categories but not others (rows ",
           paste(rows[missing_limits], collapse = ", "),
           "). Supply limits for every non-reference category, or none.",
           call. = FALSE)
    }
    if (length(missing_limits) > 0L) {
      stop("easydr: study ", study, " has no usable confidence limits.",
           call. = FALSE)
    }

    # V7: the reported effect must lie inside its own limits.
    inconsistent <- which(!is.na(te) & !is.na(lo) & !is.na(hi) &
                            (te < lo | te > hi))
    if (length(inconsistent) > 0L) {
      stop("easydr: study ", study, " reports an effect outside its own ",
           "confidence limits in rows ",
           paste(rows[inconsistent], collapse = ", "), ".", call. = FALSE)
    }
    if (any(!is.na(te) & te <= 0)) {
      stop("easydr: study ", study, " reports a non-positive ratio measure in ",
           dQuote(roles$TE), ". A risk, odds or rate ratio must exceed zero.",
           call. = FALSE)
    }

    y <- log(te)
    v <- ((log(hi) - log(lo)) / (2 * stats::qnorm(0.975)))^2
    if (!is.null(roles$event.e) && !is.null(roles$n.e)) {
      cases <- as.numeric(data[[roles$event.e]][rows])
      total <- as.numeric(data[[roles$n.e]][rows])
    }
  } else {
    m  <- as.numeric(data[[roles$mean.e]][rows])
    sdv <- as.numeric(data[[roles$sd.e]][rows])
    nn <- as.numeric(data[[roles$n.e]][rows])
    bad <- which(is.na(m) | is.na(sdv) | is.na(nn) | sdv < 0 | nn <= 0)
    if (length(bad) > 0L) {
      stop("easydr: study ", study, " has unusable mean, standard deviation or ",
           "group size in rows ", paste(rows[bad], collapse = ", "), ".",
           call. = FALSE)
    }
    # The engine takes the group mean itself and differences each category
    # against the study reference internally, deriving the covariance from the
    # standard deviations and group sizes. Differencing here first gives an
    # identical slope, so easydr does not compute an effect at all.
    y <- m
    v <- rep(NA_real_, n)
    cases <- sdv
    total <- nn
  }

  if (!identical(roles$intake, "continuous")) {
    # A log contrast against the study reference is zero by construction, and the
    # engine expects that row present with no variance.
    y[iref] <- 0
    v[iref] <- NA_real_
  }

  frame <- data.frame(
    studlab = rep(study, n),
    dose    = d,
    y       = y,
    v       = v,
    cases   = cases,
    n       = total,
    type    = rep(design, n),
    is_ref  = seq_len(n) == iref,
    stringsAsFactors = FALSE
  )
  list(frame = frame, prov = prov)
}

#' Per-contrast effect sizes for counts, through metafor
#'
#' The contrast is each category against that study's reference category. The
#' formulas are escalc's and are never written out here (research decision D4).
#'
#' @keywords internal
#' @noRd
dr_escalc_counts <- function(cases, total, iref, sm, design) {
  measure <- switch(sm, "OR" = "OR", "IRR" = "IRR", "RR" = "RR", "RR")
  n <- length(cases)
  y <- rep(NA_real_, n)
  v <- rep(NA_real_, n)
  add <- if (any(cases == 0)) 1 / 2 else 0
  for (i in seq_len(n)) {
    if (i == iref) next
    es <- if (identical(measure, "IRR")) {
      metafor::escalc(measure = "IRR", x1i = cases[i], x2i = cases[iref],
                      t1i = total[i], t2i = total[iref],
                      add = add, to = "only0")
    } else if (identical(measure, "OR")) {
      metafor::escalc(measure = "OR", ai = cases[i], bi = total[i] - cases[i],
                      ci = cases[iref], di = total[iref] - cases[iref],
                      add = add, to = "only0")
    } else {
      metafor::escalc(measure = "RR", ai = cases[i], bi = total[i] - cases[i],
                      ci = cases[iref], di = total[iref] - cases[iref],
                      add = add, to = "only0")
    }
    y[i] <- as.numeric(es$yi)
    v[i] <- as.numeric(es$vi)
  }
  list(y = y, v = v)
}
