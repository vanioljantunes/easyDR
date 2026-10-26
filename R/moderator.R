# Meta-regression on study-level covariates.
#
# The thing to get right here is not the arithmetic, it is what the numbers
# mean. The engine regresses each study's dose-response coefficients on the
# covariate, so a moderator coefficient is a change in the dose-response and
# not an effect on the outcome. For a one-term shape the engine labels the dose
# slope "(Intercept)" and the moderator effect by the covariate name, which a
# reader would naturally take for a baseline risk and a risk difference. Both
# are slopes. So easydr builds every label itself, and the printed output says
# what the coefficients are before it shows any.

#' Reduce a table to one row per study, checking the covariates can mean what
#' they appear to mean
#'
#' @keywords internal
#' @noRd
dr_mod_study_table <- function(data, studlab, vars) {
  missing <- setdiff(vars, names(data))
  if (length(missing) > 0L) {
    stop("easydr: the moderator ",
         if (length(missing) > 1L) "columns " else "column ",
         paste(dQuote(missing), collapse = ", "), " ",
         if (length(missing) > 1L) "are" else "is",
         " not in the data. Columns present: ",
         paste(names(data), collapse = ", "), ".", call. = FALSE)
  }
  studies <- unique(as.character(data[[studlab]]))
  rows <- lapply(studies, function(s) {
    sub <- data[as.character(data[[studlab]]) == s, , drop = FALSE]
    one <- lapply(vars, function(v) {
      vals <- unique(sub[[v]])
      # V19: a study-level covariate cannot vary inside a study. If it does, it
      # is a category-level covariate and this is a different analysis.
      if (length(vals) > 1L) {
        stop("easydr: the moderator ", dQuote(v), " takes ", length(vals),
             " different values inside study ", s, " (",
             paste(utils::head(as.character(vals), 4), collapse = ", "),
             "). A moderator is a characteristic of a study, so it must take ",
             "one value per study. A covariate that varies by dose category ",
             "is a different analysis and easydr does not fit it.",
             call. = FALSE)
      }
      vals[1]
    })
    names(one) <- vars
    as.data.frame(c(list(.study = s), one), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  for (v in vars) {
    # V20: a covariate with one value across every study explains nothing and
    # is not estimable. Reporting a coefficient of zero would be worse.
    if (length(unique(out[[v]])) < 2L) {
      stop("easydr: the moderator ", dQuote(v), " takes the same value (",
           as.character(out[[v]][1]), ") in all ", nrow(out),
           " studies, so it cannot explain any difference between them and ",
           "is not estimable.", call. = FALSE)
    }
  }
  out
}

#' Rule V23: warn when a continuous moderator's range rests on one study
#'
#' A coefficient driven by a single study is not wrong, but a reader cannot see
#' it from the estimate or its interval.
#'
#' @keywords internal
#' @noRd
dr_mod_narrow <- function(tab, vars) {
  narrow <- character()
  for (v in vars) {
    x <- tab[[v]]
    if (!is.numeric(x) || length(x) < 3L) next
    rng <- range(x)
    for (end in c("min", "max")) {
      at_end <- if (identical(end, "min")) x == rng[1] else x == rng[2]
      if (sum(at_end) != 1L) next
      rest <- x[!at_end]
      gap  <- if (identical(end, "min")) {
        min(rest) - rng[1]
      } else {
        rng[2] - max(rest)
      }
      if (gap > diff(range(rest))) {
        narrow <- c(narrow, v)
        warning("easydr: the range of the moderator ", dQuote(v),
                " rests on one study: ", tab$.study[at_end][1], " sits at ",
                x[at_end][1], " and the gap to the next study is wider than ",
                "the whole spread of the rest. Its coefficient is driven by ",
                "that one study.", call. = FALSE)
        break
      }
    }
  }
  unique(narrow)
}

#' Centre the continuous moderators at their mean across studies
#'
#' Two reasons, and the second one is the one that bites. A coefficient for an
#' uncentred covariate reports the dose-response at a covariate value of zero,
#' which for a calendar year is an extrapolation two thousand years out and is
#' not a quantity anyone wants. And the raw values ill-condition the second
#' stage: on the bundled data a quadratic shape with a raw year fails outright
#' with a singular system, while the centred fit succeeds.
#'
#' The moderator effect and its test are unchanged by this, because shifting a
#' covariate cannot change its slope. Only the reference-level coefficient
#' moves, from a value at zero to a value at the mean, and the label says which
#' mean it is.
#'
#' @keywords internal
#' @noRd
dr_mod_centre <- function(tab, vars) {
  centres <- list()
  for (v in vars) {
    if (!is.numeric(tab[[v]])) next
    centres[[v]] <- mean(tab[[v]])
    tab[[v]] <- tab[[v]] - centres[[v]]
  }
  list(table = tab, centres = centres)
}

#' The reference level of each categorical moderator
#'
#' The engine reports a coefficient per non-reference level and never names the
#' level they are against, so easydr reads it off the factor and prints it. A
#' coefficient without its reference level cannot be interpreted at all.
#'
#' @keywords internal
#' @noRd
dr_mod_references <- function(tab, vars) {
  out <- list()
  for (v in vars) {
    x <- tab[[v]]
    if (is.numeric(x)) next
    out[[v]] <- levels(factor(x))[1]
  }
  out
}

#' Build the labels the engine will not
#'
#' Every label is constructed from the shape terms and the moderator columns,
#' so a one-term fit reads as a slope and a change in slope rather than as an
#' intercept and a bare covariate name.
#'
#' @keywords internal
#' @noRd
dr_mod_labels <- function(coef_names, shape_terms, mod_cols, vars, one_term,
                          centres = list(), refs = list()) {
  # Where the reference coefficient sits has to be stated, and it is a
  # different thing for each kind of covariate: a level for a categorical one,
  # the mean for a continuous one that was centred.
  parts <- vapply(vars, function(v) {
    if (!is.null(refs[[v]])) return(paste0(v, " = ", refs[[v]]))
    if (!is.null(centres[[v]])) {
      return(paste0("the mean of ", v, " (",
                    format(centres[[v]], digits = 6), ")"))
    }
    paste0(v, " = 0")
  }, character(1))
  at <- if (length(parts) == 0L) {
    "the reference level"
  } else {
    paste(parts, collapse = " and ")
  }
  esc <- function(z) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", z)
  vapply(coef_names, function(nm) {
    hit <- mod_cols[vapply(mod_cols, function(mc)
      identical(mc, nm) || endsWith(nm, paste0(".", mc)), logical(1))]
    if (length(hit) == 0L) return(nm)
    mc   <- hit[which.max(nchar(hit))]
    term <- if (identical(mc, nm)) {
      NA_character_
    } else {
      sub(paste0("\\.", esc(mc), "$"), "", nm)
    }
    piece <- if (one_term || is.na(term) || !term %in% shape_terms) {
      "slope"
    } else {
      paste0("the ", term, " term")
    }
    if (identical(mc, "(Intercept)")) {
      return(paste0(piece, " at ", at))
    }
    var <- vars[vapply(vars, function(v) startsWith(mc, v), logical(1))]
    var <- if (length(var) == 0L) mc else var[which.max(nchar(var))]
    lvl <- sub(paste0("^", esc(var)), "", mc)
    if (nzchar(lvl)) {
      paste0("change in ", piece, ", ", lvl)
    } else {
      paste0("change in ", piece, " per unit of ", var)
    }
  }, character(1), USE.NAMES = FALSE)
}

# The meta-regression, and the bubble plot that draws it, are kept internal.
# The source this package teaches from does not present them, and a reader who
# met a moderator model here would have nowhere to read about it. The code and
# its tests stay: exposing it again is one roxygen tag.

#' Fit a study-level moderator through the engine
#' @keywords internal
#' @noRd
easydr_metareg <- function(object, mods, data = NULL, level = NULL) {
  inh <- dr_inherited(object)
  if (!inherits(mods, "formula")) {
    stop("easydr: `mods` must be a one-sided formula of study-level ",
         "covariates, such as ", dQuote("~ design"), ".", call. = FALSE)
  }
  if (length(mods) != 2L) {
    stop("easydr: `mods` must be one-sided: write ", dQuote("~ design"),
         " and not ", dQuote("y ~ design"), ". A moderator explains the ",
         "dose-response; it is not an outcome.", call. = FALSE)
  }
  if (is.null(data)) data <- inh$data
  if (is.null(level)) level <- inh$level
  vars <- all.vars(mods)

  raw_tab   <- dr_mod_study_table(data, inh$studlab, vars)
  narrow    <- dr_mod_narrow(raw_tab, vars)
  refs      <- dr_mod_references(raw_tab, vars)
  centred   <- dr_mod_centre(raw_tab, vars)
  study_tab <- centred$table
  centres   <- centred$centres

  mm        <- stats::model.matrix(mods, data = study_tab)
  n_studies <- nrow(study_tab)
  n_coef    <- inh$shape$terms * ncol(mm)
  # V21: the second stage needs more studies than coefficients, or there is
  # nothing left to estimate the between-study variance with.
  if (n_coef >= n_studies) {
    stop("easydr: this moderator model needs ", n_coef, " coefficients (",
         inh$shape$terms, " shape term",
         if (inh$shape$terms == 1L) "" else "s", " times ", ncol(mm),
         " moderator column", if (ncol(mm) == 1L) "" else "s",
         ") and there are only ", n_studies, " studies. It needs more studies ",
         "than coefficients. Fit fewer moderators, or a simpler shape.",
         call. = FALSE)
  }

  keep  <- c("dose", "studlab", "event.e", "n.e", "TE", "lower", "upper",
             "ref", "design", "mean.e", "sd.e")
  roles <- dr_roles(data, as.list(inh$args[intersect(names(inh$args), keep)]),
                    kind = if (identical(inh$entry, "easydrcont")) {
                      "cont"
                    } else {
                      "bin"
                    })
  prep <- dr_prepare(data, roles, inh$sm)
  shp  <- dr_shape(inh$shape, prep$frame$dose, studies = prep$frame$studlab)

  # The covariates are study-level, so they ride on every row of their study.
  # The engine takes them from the data frame it is handed.
  idx      <- match(as.character(prep$frame$studlab), study_tab$.study)
  mod_data <- study_tab[idx, vars, drop = FALSE]
  rownames(mod_data) <- NULL

  fit <- dr_fit(prep$frame, sm = inh$sm, intake = roles$intake,
                proc = inh$proc, method = inh$method,
                covariance = inh$covariance,
                method.smd = inh$args$method.smd, shape = shp,
                mod = mods, mod_data = mod_data)

  b     <- coef(fit)
  se    <- sqrt(diag(as.matrix(vcov(fit))))
  zcrit <- stats::qnorm(1 - (1 - level) / 2)
  labels <- dr_mod_labels(names(b), colnames(shp$basis(prep$frame$dose)),
                          colnames(mm), vars, one_term = shp$terms == 1L,
                          centres = centres, refs = refs)
  is_mod <- !vapply(names(b), function(nm)
    identical(nm, "(Intercept)") || endsWith(nm, ".(Intercept)"), logical(1))

  terms_tab <- data.frame(
    coefficient = names(b), label = labels, estimate = as.numeric(b),
    se = as.numeric(se), lb = as.numeric(b) - zcrit * as.numeric(se),
    ub = as.numeric(b) + zcrit * as.numeric(se), moderator = unname(is_mod),
    stringsAsFactors = FALSE
  )
  rownames(terms_tab) <- NULL

  structure(
    list(
      fit                  = fit,
      terms                = terms_tab,
      test                 = dr_moderator_test(fit, which(is_mod)),
      heterogeneity        = dr_heterogeneity(fit, FALSE),
      heterogeneity_before = object$heterogeneity,
      reference_levels     = refs,
      narrow               = narrow,
      moderators           = vars,
      types                = vapply(vars, function(v)
        if (is.numeric(study_tab[[v]])) "continuous" else "categorical",
        character(1)),
      study_table          = raw_tab,
      centres              = centres,
      inherited            = inh,
      sm                   = object$sm,
      sm_label             = object$sm_label,
      ratio                = object$ratio,
      dose_unit            = object$dose_unit,
      level                = level,
      call                 = match.call()
    ),
    class = "easydr_moderator"
  )
}

#' Print a moderator analysis
#'
#' Registered so the object prints, but given no help page: the analysis it
#' belongs to is internal.
#'
#' @keywords internal
#' @noRd
#' @exportS3Method print easydr_moderator
print.easydr_moderator <- function(x, digits = 6, ...) {
  fmt <- function(v) dr_fmt(v, digits)

  cat("Moderator analysis (", x$sm_label, ")\n", sep = "")
  for (v in x$moderators) {
    ref <- x$reference_levels[[v]]
    cat("Moderator: ", v, " (", x$types[[v]],
        if (!is.null(ref)) paste0(", reference level ", ref) else "",
        ")\n", sep = "")
  }
  cat("Inherited: ", dr_inherited_line(x$inherited), "\n\n", sep = "")
  cat("A moderator coefficient is an effect on the dose-response, not on the\n")
  cat("outcome.\n\n")

  body <- data.frame(
    term = x$terms$label,
    estimate = fmt(x$terms$estimate),
    ci = paste0(fmt(x$terms$lb), " to ", fmt(x$terms$ub)),
    stringsAsFactors = FALSE
  )
  names(body) <- c("term", "estimate", paste0(format(100 * x$level), "% CI"))
  print(body, row.names = FALSE, right = FALSE)

  cat("\n")
  if (isTRUE(x$test$estimable)) {
    cat("Moderator test: X2 = ", format(round(x$test$chi2, 2), nsmall = 2),
        " (df = ", x$test$df, "), p = ",
        format.pval(x$test$p, digits = 3, eps = 1e-4), "\n", sep = "")
    cat("  Wald test on the ", length(x$test$terms), " moderator column",
        if (length(x$test$terms) == 1L) "" else "s", "\n", sep = "")
  } else {
    cat("Moderator test: not estimable\n  ", x$test$reason, "\n", sep = "")
  }

  cat("\n")
  h <- x$heterogeneity
  if (isTRUE(h$estimable)) {
    cat("Heterogeneity: Q = ", format(round(h$Q, 2), nsmall = 2), " (df = ",
        h$df, "), p = ", format.pval(h$p, digits = 3, eps = 1e-4), "\n",
        sep = "")
  } else {
    cat("Heterogeneity: not estimable\n  ", h$reason, "\n", sep = "")
  }
  b <- x$heterogeneity_before
  if (isTRUE(b$estimable)) {
    cat("  without the moderator: Q = ", format(round(b$Q, 2), nsmall = 2),
        " (df = ", b$df, "), p = ",
        format.pval(b$p, digits = 3, eps = 1e-4), "\n", sep = "")
  }

  if (length(x$centres) > 0L) {
    cat("
Defaults applied:
")
    for (v in names(x$centres)) {
      cat("  - ", v, " was centred at its mean across studies (",
          format(x$centres[[v]], digits = 6), "), so the reference ",
          "coefficient is the dose-response there rather than at ", v,
          " = 0. The moderator effect and its test are unchanged by this.
",
          sep = "")
    }
  }

  if (length(x$narrow) > 0L) {
    cat("\nThe range of ", paste(x$narrow, collapse = ", "),
        " rests on one study, so that coefficient is driven by it.\n", sep = "")
  }

  cat("\nThis is an observational association across studies, not a ",
      "randomised\ncomparison. A moderator that explains heterogeneity has ",
      "not been shown to\ncause it.\n", sep = "")
  invisible(x)
}
