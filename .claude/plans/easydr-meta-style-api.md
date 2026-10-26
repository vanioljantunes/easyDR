# easyDR meta-style API

Progress: 19/38 weighted (50%)
Started: 2026-10-26 BRT   Updated: 2026-10-26 BRT

## Goal

One exported family in the easyDR repo, named and argued like the meta package:
easydr(), easydr_select(), easydr_forest(), easydr_curve(), easydr_loo(),
easydr_metareg(), easydr_bubble(), easydr_funnel(). The caller passes sm = "RR"
and a long table of dose levels; easydr computes the contrasts, calls
dosresmeta once, and folds the viability check and the curve label into the
returned object. dr_viability() and curve_shape() stop being exported.
Shape is always explicit: no auto-selection inside easydr().
Done when R CMD check is 0/0/0, the Orsini 2012 anchors reproduce to the digits
already published in the README, and every function appears in the README with
a rendered picture.

## Decisions taken

- Home: the easyDR repo (not metarank4). Code relicensed MIT -> GPL-2 because
  dosresmeta is GPL-2 and moves to Imports. Book stays CC BY-NC-SA 4.0.
- Shape always explicit. The AIC loop stays its own call, easydr_select().
- Publication bias: funnel of study-specific trend estimates against their
  standard errors plus Egger regression. Deeks is diagnostic-accuracy only.

## Steps

### Phase 1 - repo and licence
- [x] 1. DESCRIPTION: GPL-2, dosresmeta/metafor/mixmeta/Hmisc to Imports
- [x] 2. LICENSE replaced with GPL-2 pointer; LICENSE-text.md untouched
- [ ] 3. README licence section rewritten to say GPL-2 for code
- [x] 4. .Rbuildignore covers tutorial, book, docs, .claude

### Phase 2 - port the validated engine (w=2 each)
- [x] 5. port prepare.R engine.R shapes.R result.R interval.R (w=2)
- [x] 6. port inherit.R select.R sensitivity.R moderator.R (w=2)
- [x] 7. port data/ and data-raw/ and R/data.R (w=2)
- [x] 8. port tests/ and inst/validation/ with calls rewritten to easydr() (w=2)
- [x] 9. devtools::document + load_all clean

### Phase 3 - easydr()
- [x] 10. easydr(): sm drives the role set, shape required, one call
- [x] 11. fold viability in as internal, carried on the result
- [x] 12. fold curve label in as internal, carried on the result
- [x] 13. print method shows viability line and curve label
- [x] 14. easydr_select() exported, AIC loop, warns it is one-stage ML only
- [x] 15. tests: easydr() reproduces the Orsini anchors exactly
- [x] 16. tests: shape omitted is an error, never a silent default

### Phase 4 - the plots
- [ ] 17. easydr_forest(): per-study trend, pooled diamond, meta layout
- [ ] 18. easydr_curve(): main, xlab, ylab, xlim, ylim, round, prediction
- [ ] 19. easydr_loo(): wraps loo_dr, returns table plus refits
- [ ] 20. easydr_metareg(): wraps moderator_dr, meta-style printout
- [ ] 21. easydr_bubble(): moderator on x, study trend on y, size = precision
- [ ] 22. easydr_funnel(): trend against SE plus Egger regression test
- [ ] 23. tests for all six, including the Egger numbers
- [ ] 24. snapshot-free plot tests (return value, not pixels)

### Phase 5 - docs and pictures
- [ ] 25. roxygen pages for the eight exported functions
- [ ] 26. tutorial/plot.R extended: forest, curve, bubble, funnel figures
- [ ] 27. tutorial/api.html: easydr() and easydr_select()
- [ ] 28. tutorial/plots.html: the five plot functions
- [ ] 29. render both to man/figures
- [ ] 30. README rewritten around the new API with the new pictures
- [ ] 31. README validation section rerun against easydr(), numbers verified

### Phase 6 - green and shipped
- [ ] 32. devtools::test 0 failures
- [ ] 33. R CMD check 0 errors 0 warnings 0 notes
- [ ] 34. commit and push to origin main

## Notes

Port verdict: of the 32 recorded baseline cases, 31 printed byte-identically and
one differed by a single line, the adoption hint that now reads
easydr(data, ..., shape = "linear"). No number moved. The 18 recorded refusals
all still refuse with the same message.

easydr() vs the raw dosresmeta call in the old README: slope 0.0064373 against
0.0064375, because easydr pools the reported rate ratio and uses the counts only
for the covariance, whereas the old chunk pooled dosresmeta's own logrr column.
Rate ratios agree to five decimals. The README validation numbers have to be
regenerated from easydr(), not copied across.

Weights: steps 5-8 are w=2 (ports of 3500 lines plus 18 test files).
Total weight 38; percentage is weighted.
