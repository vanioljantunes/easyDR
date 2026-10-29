# Multi-shape curves, and alignment with Crippa and Orsini (2016)

Progress: 17/47 (36%)
Started: 2026-10-29 16:10 BRT   Updated: 2026-10-29 18:30 BRT

## Goal

Two features, requested back to back, that meet in the same section of
Chapter 4.

**A. Several shapes in one call, and one figure.** `easydr()` accepts a vector
of shapes and returns an `easydr_multi` holding one fit per shape, in the
order given. `easydr_curve()` gains `mult.curves` and `curve.panels` and draws
them as panels on shared axes with one summary table underneath. Panel
layouts: 2 side by side, 3 side by side, 4 as two rows of two, 5 as three then
two, 6 as two rows of three. The other functions refuse a multi object and say
how to pick one fit.

**B. Crippa and Orsini (2016).** The Emax family, target dose estimation with
bootstrap intervals, the continuous covariance formulas, goodness-of-fit
reporting, and the matching checklist items in Chapter 5.

Decided with the author:

- A multi object; the other functions refuse it rather than silently using the
  first fit.
- Panels share axes, and one table sits under the whole figure.
- Emax is fitted by profiling, and every interval for the Emax family comes
  from a parametric bootstrap that re-profiles on each draw, so the
  uncertainty in the profiled parameter is carried rather than ignored.

Done means: the suite is green, the book renders, the saved figures are
inspected by eye, README and both chapters are updated, committed and pushed.

## Steps

### A. Several shapes, one figure

- [x] 1. `easydr()` takes a shape vector and returns `easydr_multi` (w=2)
- [x] 2. `print.easydr_multi`: one row per shape, terms, departure from
      linearity, deviance, AIC
- [x] 3. `dr_require_single()` guard on forest, funnel, loo and table, naming
      how to reach one fit
- [x] 4. Split `easydr_curve()` so the panel is built separately from being
      drawn (w=2)
- [x] 5. `mult.curves` and `curve.panels`, shared axes, the six layouts (w=3)
- [x] 6. One summary table under the panels, a column per shape
- [x] 7. File export on the panel figure, size from the row and column count
- [x] 8. Tests, and the constitution test updated (w=2)
- [x] 9. Chapter 4: the six-shapes section uses the feature (w=2)
- [x] 10. README: the feature, with the layout table

### B. The Emax family

- [ ] 11. `dr_shape_hyperbolic_emax()`: linear in d/(d+ED50) for a fixed
      ED50, profiled over a grid on deviance (w=3)
- [ ] 12. `dr_shape_emax()`: the same over a two-dimensional grid of ED50 and
      steepness (w=3)
- [ ] 13. Parametric bootstrap that re-profiles each draw, and replaces the
      reported intervals for the Emax family (w=3)
- [ ] 14. The two shapes reachable from `easydr()`, `easydr_curve()` and the
      shape table
- [ ] 15. Correct the comment in `R/curve_shape.R` that says the package does
      not fit Emax, which this feature makes untrue

### C. Target dose

- [ ] 16. `ed_gamma()`: the dose reaching a stated fraction of the maximum
      effect inside the observed range (w=2)
- [ ] 17. `easydr_target_dose()`: ED50, ED80 and EDmax with parametric
      bootstrap intervals (w=3)
- [ ] 18. Exported, documented, and refused for a shape with no maximum
- [ ] 19. `easydr_table()` gains optional target-dose rows (w=2)

### D. Goodness of fit

- [ ] 20. Deviance against the saturated model, R-squared and adjusted
      R-squared in `print()` and `summary()` (w=2)
- [ ] 21. The same statistics available to `easydr_table()`

### E. Tests

- [ ] 22. Emax family: fit, profile, bootstrap, refusals (w=3)
- [ ] 23. Target dose: known closed-form cases, bootstrap behaviour,
      refusals (w=3)
- [ ] 24. Constitution and grammar tests updated for the new exports (w=2)

### F. The book

- [ ] 25. Chapter 4: covariance reconstruction for MD and SMD, formulas and
      code (w=2)
- [ ] 26. Chapter 4: target dose, the definition and the bootstrap (w=2)
- [ ] 27. Chapter 4: the Emax family in the model comparison panel
- [ ] 28. Chapter 5: Table 5.1 gains the target-dose and goodness-of-fit
      items (w=2)
- [ ] 29. Chapter 5, section 5.3: bounded and ordinal scales, with the PANSS
      rescaling example
- [ ] 30. Chapter 5, section 5.7: knot-location sensitivity across the 10th,
      25th, 50th, 75th and 90th percentiles
- [ ] 31. README, full suite, render, figures inspected, commit, push (w=2)

## Blocked

Step 28 and the rest of section F wait on one thing: the third request ended
at the heading "### 3. Chapter 5 & Checklist Updates" with nothing under it.
The chapter 5 edits from the second request are known and can go ahead; the
ones this heading was about are not.

A fourth group of work arrived with that request and is not yet in the step
list above: `easydr_compare_shapes()`, `easydr_knot_sensitivity()`, the
"A Structured Framework for Model Choice" section of chapter 4, and the
AIC/BIC treatment of non-nested shapes. The comparison table that
`print.easydr_multi` already prints carries AIC, BIC, deviance and
p-nonlinearity, so `easydr_compare_shapes()` is a thin wrapper over the multi
object rather than new machinery.

## Notes

Weights in brackets; unmarked steps count 1. Total 47.

`dosresmeta` 2.2.0 exports `covar.smd` and `gof`, which are the Crippa and
Orsini continuous covariance and goodness-of-fit routines, so steps 20, 21 and
25 build on the engine rather than reimplementing it. It has no Emax, so the
Emax family is easyDR's own.

Another session commits to this repo. Fetch before each commit.
