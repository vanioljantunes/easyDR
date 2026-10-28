# easyDR: the cubic polynomial and the one-power fractional polynomial

Progress: 15/15 (100%)
Started: 2026-10-28 02:55 BRT   Updated: 2026-10-28 03:40 BRT

## Goal

`shape = "cubic"` fits dose, dose squared and dose cubed, closing the
polynomial family that stopped at degree 2. `powers` accepts a single value as
well as a pair, so a one-power fractional polynomial can be fitted; Royston and
Altman treat FP1 as the first candidate, and the current code refuses it with
"`powers` must be a pair".

The operator found both gaps: the package named five shapes, and the only
"cubic" in it was the restricted cubic spline, which is a different model.

After this there are six named shapes: linear, log, quadratic, cubic, rcs, fp.

Not doing, by decision: a linear spline for thresholds, and an unrestricted
B-spline. The B-spline duplicates `rcs`.

## Steps

- [x] 1. R/shapes.R: dr_shape_cubic(), dose first so the Wald test stays a test of departure from linearity
- [x] 2. R/shapes.R: add "cubic" to match.arg and the switch
- [x] 3. R/shapes.R: let dr_shape_fp() take one power or two, terms follows (w=2)
- [x] 4. R/easydr.R: the missing-shape error lists six shapes
- [x] 5. roxygen for shape and powers in easydr.R, then document()
- [x] 6. tests: cubic fits, nests the quadratic, FP1 fits, FP with 3 powers still refused (w=2)
- [x] 7. inst/validation: re-record the one refusal whose message lists the shapes
- [x] 8. tutorial/figures.R: the fitted grid becomes six panels
- [x] 9. tutorial/api.html: the shapes table gains cubic, "five" becomes six
- [x] 10. README section 2: the shapes table, the grid caption and the "all five" sentence
- [x] 11. devtools::test() and R CMD check to 0/0/0
- [x] 12. read the rebuilt grid figure back
- [x] 13. commit and push
- [x] 14. (added) fix the slope label on one-term nonlinear shapes, found while adding FP1
- [x] 15. (added) re-record the two baseline answers that changed, guarded

## Notes

Every file that lists the shapes by name, from the grep: `R/easyDR-package.R`,
`R/easydr.R`, `R/shapes.R`, nine test files, `tutorial/api.html`,
`tutorial/figures.R`, `tutorial/_out/shape-error.txt`,
`inst/validation/cases.R`, `README.md`, and two `man/` pages that
`document()` regenerates.

A cubic on a dose reaching 60 puts design values near 2e5, which is the
conditioning problem `dr_shape_fp()` already solves by scaling. Watch for it:
if the fit is ill-conditioned, scale the cubic the same way rather than
leaving a fit that looks hundreds of AIC better than it is.
