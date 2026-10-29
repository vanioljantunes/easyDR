# What a sensitivity analysis inherits, and how a refit is assembled.
#
# All three analyses of this feature answer a question about a fit that already
# exists. None of them is allowed to change it. The risk is not that a refit
# fails loudly, it is that a refit succeeds with a model nobody asked for: a
# spline whose knots were re-placed on fewer doses, or a column auto-detected
# differently because a subset left it empty. Both would print cleanly and both
# would be wrong. So the settings are read off the fitted object and handed
# back, never derived again.

#' The one path from a table to a fitted result
#'
#' Both entry points did the same five things in the same order: check the
#' level, translate `method.tau`, resolve the column roles, prepare the
#' analysis frame, resolve the shape, fit, and build the result. Only the
#' column roles and two arguments differed. Keeping two copies of that
#' sequence meant any change to the pipeline had to be made twice, which is
#' how the two entry points would eventually have come to disagree about
#' something neither page documented.
#'
#' The public functions keep their own signatures, defaults and `match.arg`
#' calls: those are the grammar a reader knows from `meta`, and they differ
#' between the two on purpose. The `call` is passed in from the public
#' function so the printed call is the one the user wrote.
#'
#' @keywords internal
#' @noRd
dr_entry <- function(kind, data, call, roles_args, sm, method.tau, proc,
                     covariance, method.smd, shape, knots, powers,
                     ed50 = NULL, steepness = NULL, grid, at,
                     dose.unit, level, common, random) {
  dr_check_level(level)
  method <- dr_method_tau(method.tau, common, random)

  roles <- dr_roles(data, roles_args, kind = kind)
  prep  <- dr_prepare(data, roles, sm)
  shp   <- dr_shape(shape, prep$frame$dose, knots = knots, powers = powers,
                    ed50 = ed50, steepness = steepness,
                    studies = prep$frame$studlab)
  # A shape with a parameter that is not linear in the model, which is the
  # Emax family, is resolved by fitting over a grid and keeping the best.
  # Every other shape passes straight through.
  shp   <- dr_profile_shape(shp, prep$frame, sm = sm, intake = roles$intake,
                            proc = proc, method = method,
                            covariance = covariance, method.smd = method.smd)
  fit   <- dr_fit(prep$frame, sm = sm, intake = roles$intake, proc = proc,
                  method = method, covariance = covariance,
                  method.smd = method.smd, shape = shp)

  dr_result(fit, prep, sm = sm, dose_unit = dose.unit, at = at, level = level,
            proc = proc, method = method, call = call, shape = shp,
            grid = grid, data = data,
            entry = if (identical(kind, "cont")) "easydrcont" else "easydrbin",
            args = dr_entry_args(roles, sm, method.tau, proc, covariance, shp,
                                 dose.unit, level, common, random,
                                 method.smd = method.smd))
}

#' Collect the arguments a refit has to reuse
#'
#' Called by the entry points so the result carries everything needed to call
#' the same entry point again. The columns are recorded by the name that was
#' resolved, not left to be auto-detected a second time: a subset can make a
#' column look different from the one the full fit used.
#'
#' @keywords internal
#' @noRd
dr_entry_args <- function(roles, sm, method.tau, proc, covariance, shape,
                          dose.unit, level, common, random,
                          method.smd = NULL) {
  cols <- roles[setdiff(names(roles), "intake")]
  cols <- cols[!vapply(cols, is.null, logical(1))]
  model <- list(sm = sm, method.tau = method.tau, proc = proc, shape = shape,
                dose.unit = dose.unit, level = level, common = common,
                random = random)
  if (!is.null(covariance)) model$covariance <- covariance
  if (!is.null(method.smd)) model$method.smd <- method.smd
  c(cols, model)
}

#' The settings a sensitivity analysis inherits from a fitted result
#'
#' One place reads them, so the three analyses cannot drift apart in what they
#' carry. The shape is the resolved object rather than its name, which is what
#' keeps a refit on the fitted basis.
#'
#' @keywords internal
#' @noRd
dr_inherited <- function(object) {
  if (!inherits(object, "easydr_result")) {
    stop("easydr: this analysis takes a fitted result from easydrbin() or ",
         "easydrcont().", call. = FALSE)
  }
  if (is.null(object$data) || is.null(object$entry) || is.null(object$args)) {
    stop("easydr: this result was produced before refitting was supported, so ",
         "the table it was fitted from is not available. Refit it with ",
         "easydrbin() or easydrcont() and try again.", call. = FALSE)
  }
  p <- object$provenance
  list(
    entry      = object$entry,
    shape      = object$shape,
    sm         = object$sm,
    proc       = p$proc,
    method     = p$method,
    covariance = p$covariance,
    reference  = object$reference$curve,
    level      = p$level,
    dose_unit  = object$dose_unit,
    studlab    = object$args$studlab,
    studies    = unique(as.character(object$data[[object$args$studlab]])),
    args       = object$args,
    data       = object$data,
    call       = object$call
  )
}

#' One line naming everything a refit kept
#'
#' Printed by all three analyses, so a reader can see the inherited settings
#' instead of being asked to trust that they were inherited.
#'
#' @keywords internal
#' @noRd
dr_inherited_line <- function(inh) {
  unit <- if (is.null(inh$dose_unit)) "dose unit" else inh$dose_unit
  paste0(inh$shape$label, "; ",
         if (identical(inh$proc, "1stage")) "one-stage " else "two-stage ",
         tolower(inh$method), "; reference ", inh$reference, " ", unit)
}

#' Refit the inherited model on a subset of the studies
#'
#' The only thing that differs from the original call is which studies are
#' present. The shape travels as an object, so its knots, powers, shift and
#' scale are the fitted ones; the columns travel by name; the pooling stage and
#' the covariance approximation are the ones that were used. A subset the shape
#' no longer suits raises its own error, which the caller records against that
#' row rather than swallowing.
#'
#' @param extra arguments to override for this refit only, used for the dose a
#'   curve is compared at. Nothing about the model may be overridden this way.
#'
#' @keywords internal
#' @noRd
dr_refit <- function(inh, keep, extra = NULL) {
  sub <- inh$data[as.character(inh$data[[inh$studlab]]) %in% keep, ,
                  drop = FALSE]
  args <- inh$args
  if (length(extra) > 0L) {
    forbidden <- intersect(names(extra),
                           c("sm", "proc", "method.tau", "shape", "covariance",
                             "method.smd", "level"))
    if (length(forbidden) > 0L) {
      stop("easydr: a refit may not override ",
           paste(forbidden, collapse = ", "),
           ". Those are inherited from the fit.", call. = FALSE)
    }
    args[names(extra)] <- extra
  }
  do.call(inh$entry, c(list(data = sub), args))
}
