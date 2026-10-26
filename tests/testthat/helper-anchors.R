# Reference values for the sensitivity and moderator analyses.
#
# Every value here is a captured run of the engine on the bundled alcohol_cvd
# table, R 4.6.0 with dosresmeta 2.2.0, captured 2026-10-24. None of them is a
# published figure, and none should be cited as one. They prove easydr agrees
# with dosresmeta; they do not prove either agrees with a paper. Where a
# publication reports the same analysis, its figure should replace the captured
# value here.
#
# The linear paths of features 001 and 002 are anchored to published results
# instead: Orsini et al. (2012) for binary outcomes, Crippa and Orsini (2016)
# for continuous ones. Those live in test-published-bin.R and
# test-published-cont.R.

# Studies in the order they appear in the data, which is the default order the
# cumulative analysis uses.
cvd_studies <- c("Bianchi", "Bobak", "Malarcher-wine", "Malarcher-beer",
                 "Vliegenthart-wine", "Vliegenthart-beer")

# The full linear fit on the log scale.
cvd_full_slope <- -0.0043654117582
cvd_full_Q     <- 14.1463500517
cvd_full_df    <- 5L
cvd_full_Psi   <- 1.005709157e-4

# One slope per leave-one-out refit, named by the study left out. Bianchi is
# the influential study: removing it nearly triples the slope, which is why
# this table is the one the sensitivity examples use.
cvd_loo_slopes <- c(
  "Bianchi"           = -0.0117619010539,
  "Bobak"             = -0.000827522339765,
  "Malarcher-wine"    = -0.00428650150779,
  "Malarcher-beer"    = -0.00509218818642,
  "Vliegenthart-wine" = -0.00264041801244,
  "Vliegenthart-beer" = -0.00321896617509
)

# The study design as a moderator, which is a real study-level covariate on
# this table: four case-control studies and two of cumulative incidence. The
# reference level is case-control, being first alphabetically and so first in
# levels(). The engine names the second coefficient
# "designcumulative incidence", which is why easydr builds its own labels.
cvd_mod_reference    <- "case-control"
cvd_mod_slope_ref    <- -0.0006838
cvd_mod_slope_change <- -0.0099074
cvd_mod_wald_chi2    <- 0.591579
cvd_mod_wald_df      <- 1L
cvd_mod_wald_p       <- 0.441809
cvd_mod_Q            <- 13.01563
cvd_mod_df           <- 4L
cvd_mod_p            <- 0.0111997
