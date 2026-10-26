# easyDR round 2: published shape names, no select, easyDTA-grade plots

Progress: 27/28 weighted (96%)
Started: 2026-10-26 BRT   Updated: 2026-10-26 BRT

## Goal

Four corrections the operator asked for, verified rather than asserted:
the curve labels become the ones the literature actually publishes;
easydr_select() is deleted because it confused more than it helped; the
tutorial teaches easydr_curve() well enough to replicate from the picture;
and the four plots are rebuilt to the easyDTA standard, with metafor as the
reference for the funnel and the bubble. Every picture carries its code
beside the graphic, and every README picture is followed by the chunk.
Done when R CMD check is 0/0/0 and no page claims anything unverified.

## What the research found

The labels did NOT follow Crippa and Orsini, which the help page claimed.
The published trichotomy is plateau after an initial increase, still
increasing at the maximum dose, and bell-shaped. Emax is a pharmacology
MODEL name that Crippa and Orsini compare splines against, not a shape
label. The claim is removed, not reworded.

Operator chose: monotonic increasing, monotonic decreasing, plateau,
bell-shaped, U-shaped, flat. Bell and U were previously collapsed.

## What went wrong on the way

Three silent no-op edits: a python replace whose anchor had already changed
printed "ok" without asserting, so the funnel kept drawing an old footer twice
over. Every later edit script asserts the replacement landed.

Defects the rendering found, none of which a test would have caught:
the leave-one-out table drew on top of its own plot, because a null unit in
the stack collapses the viewport height; the forest's effect panel carried an
axis title the text panels did not, so the panel areas differed and the pooled
diamond sat half a row high; the I-squared line overflowed its panel; the
categorical bubble clipped the very interval that ci = TRUE had added.

The superscript two went in as literal UTF-8 and failed R CMD check. It is now
"I-squared" and "X2", which is what the package's print methods already say.

## Steps

### Phase 1 - the labels
- [x] 1. curve_shape(): six labels, bell-shaped split from U-shaped
- [x] 2. remove the unsupported attribution; cite Leucht 2020 and Crippa 2016
- [x] 3. tests: one per label, plus the bell/U split
- [x] 4. fig-shapes.png redrawn with six panels

### Phase 2 - delete easydr_select
- [x] 5. delete R/select.R, its tests and its man page
- [x] 6. easydr() missing-shape error lists the shapes, names no other function
- [x] 7. drop the select case from cases.R and from the baseline
- [x] 8. package page, NAMESPACE and the export-surface test follow

### Phase 3 - plots to the easyDTA standard (w=2 each)
- [x] 9. ggplot2 and gridExtra to Imports (w=2)
- [x] 10. easydr_forest(): panelled, zebra rows, diamond, I-squared footer (w=2)
- [x] 11. easydr_loo(): gains the plot, curves per omitted study plus table (w=2)
- [x] 12. easydr_funnel(): contour-enhanced 90/95/99, legend, stats table (w=2)
- [x] 13. easydr_bubble(): metafor regplot style, band on the line (w=2)
- [x] 14. easydr_curve(): arguments documented to replicate from the picture
- [x] 15. tests updated for every changed return value

### Phase 4 - the pages
- [x] 16. every tutorial picture shows the code beside the graphic
- [x] 17. api.html gains the easydr_curve() call and its argument table
- [x] 18. README: a copy-paste chunk after every picture
- [x] 19. regenerate figures and console blocks from one run
- [x] 20. render the pages, read every one back
- [x] 21. devtools::test 0 failures
- [x] 22. R CMD check 0/0/0
- [ ] 23. commit and push

## Notes

Weights: steps 9 to 13 are w=2. Total weight 28.
