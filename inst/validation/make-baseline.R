# Build the baseline snapshot.
#
# Run this ONCE, on code whose answers are known good, before a refactor
# starts. After that, tests/testthat/test-invariance.R compares every run
# against what this wrote.
#
#   Rscript inst/validation/make-baseline.R
#
# Re-running it after a refactor would overwrite the evidence with whatever the
# refactored code now produces, which is the one thing that must not happen.
# The script therefore refuses to overwrite an existing baseline unless
# EASYDR_BASELINE_OVERWRITE is set, and records the git commit it was built
# from so a disagreement can be traced to a known state of the code.

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/validation/cases.R")

out_path <- "inst/validation/baseline.rds"
if (file.exists(out_path) &&
    !nzchar(Sys.getenv("EASYDR_BASELINE_OVERWRITE"))) {
  stop("easydr: ", out_path, " already exists. A baseline records what the ",
       "code answered BEFORE a change, so overwriting it destroys the only ",
       "thing the comparison rests on. Set EASYDR_BASELINE_OVERWRITE=1 if ",
       "you genuinely intend to replace it.", call. = FALSE)
}

sha <- tryCatch(
  trimws(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)),
  error = function(e) NA_character_
)

cases <- dr_validation_cases()
snap  <- vector("list", length(cases))
for (i in seq_along(cases)) {
  k <- cases[[i]]
  cat(sprintf("[%2d/%d] %s\n", i, length(cases), k$label))
  obj <- withCallingHandlers(
    k$f(),
    warning = function(w) invokeRestart("muffleWarning")
  )
  snap[[i]] <- list(
    label   = k$label,
    printed = dr_validation_print(obj),
    numbers = dr_validation_numbers(obj),
    classes = class(obj)
  )
}
names(snap) <- vapply(cases, function(k) k$label, character(1))

errs <- dr_validation_errors()
emsg <- vector("list", length(errs))
for (i in seq_along(errs)) {
  k <- errs[[i]]
  cat(sprintf("[err %2d/%d] %s\n", i, length(errs), k$label))
  emsg[[i]] <- list(
    label   = k$label,
    message = tryCatch({
      k$f()
      NA_character_
    }, error = conditionMessage)
  )
}
names(emsg) <- vapply(errs, function(k) k$label, character(1))

baseline <- list(
  built_on   = format(Sys.Date()),
  commit     = sha,
  r_version  = paste(R.version$major, R.version$minor, sep = "."),
  dosresmeta = format(utils::packageVersion("dosresmeta")),
  easyDR     = format(utils::packageVersion("easyDR")),
  results    = snap,
  errors     = emsg
)
saveRDS(baseline, out_path, version = 2)

cat("\nwrote", out_path, "\n")
cat("commit", sha, "| cases", length(snap), "| error messages",
    length(emsg), "\n")
