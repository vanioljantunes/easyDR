# easyDR subgroups

Progress: 32/34 (94%)
Started: 2026-10-27 01:10 BRT   Updated: 2026-10-28 02:20 BRT

## Goal

`easydr(data, subgroup = "design")` takes a column name and fits the model
again inside each level of it, inheriting every model argument. `easydr_forest()`
and `easydr_curve()` gain `subgroup = TRUE/FALSE`: TRUE draws the levels, FALSE
draws the overall fit alone. The curve's corner result box moves out of the
panel and becomes a table under the figure, below the footer, the way the
funnel already puts Egger's numbers under its panel. With subgroups that table
has one row per quantity, one column per level, and a last column for the
overall. Curves differ by colour and line type, marked points by shape. The
forest groups its rows the way `meta` does, with `subgroup.name = FALSE`, so a
block header reads `case-control`, not `design = case-control`.

Decided with the operator:

- Subgroup estimates come from a **stratified refit** per level, reusing the
  inheritance `easydr_loo()` already uses. Each level keeps its own tau.
- **No test for subgroup differences.** The table reports the estimates side
  by side and nothing more.
- A subgroup must be **constant within a study**, the rule moderators already
  follow. A column that changes inside a study is refused.
- A level that **cannot be fitted** becomes a row carrying the reason, the way
  a failed leave-one-out refit does. It does not kill the call.
- `subgroup = TRUE` is the **default** in both plots when the fit carries one.

Demo data: `alcohol_cvd`, column `design`, levels case-control (4 studies) and
cumulative incidence (2 studies). Both levels fit on their own: slopes
-0.00044816 and -0.010657.

The README is rebuilt around six sections the operator set, in this order:
install and what it is with every main function and its main arguments; data
preparation for binary and continuous outcomes with a worked `easydr()`,
`easydr_curve()` and `easydr_forest()`; subgroups; leave-one-out; publication
bias; validation. The five shapes appear as one grid figure of five fitted
panels, drawn with `result = FALSE` so each panel is a curve and a title
rather than a curve and an unreadable table. It is written once, after the
subgroup feature exists, because section 3 has nothing to show until then.

## Steps

- [x] 1. R/subgroup.R: validate the column, refuse one that varies within a study
- [x] 2. R/subgroup.R: dr_subgroups() stratified refits, reason rows on failure (w=2)
- [x] 3. easydr(): the `subgroup` argument, wired in, stored on the result
- [x] 4. print.easydr_result(): name the subgroup column and its levels
- [x] 5. plot-utils: dr_summary_table(), quantities by level with an overall column (w=2)
- [x] 6. plot-curve: corner box out of the panel, table under the footer (w=2)
- [x] 7. plot-curve: one curve per level, colour and line type, point shapes (w=3)
- [x] 8. plot-curve: legend naming the levels
- [x] 9. plot-forest: subgroup blocks, meta style, subgroup.name = FALSE (w=3)
- [x] 10. plot-forest: a subtotal diamond per level, overall underneath
- [x] 11. tests: validation, refits, both plots, and the no-subgroup path unchanged (w=2)
- [x] 12. roxygen for every changed argument, then document()
- [x] 13. tutorial/figures.R: a subgroup forest and a subgroup curve
- [x] 14. tutorial/figures.R: one grid figure, five fitted shapes, result = FALSE (w=2)
- [x] 15. README 1: install, what it is, every main function and its main arguments
- [x] 16. README 2: data prep binary and continuous, easydr, curve, forest, the shape grid (w=2)
- [x] 17. README 3: subgroups, forest and curve
- [x] 18. README 4: leave-one-out
- [x] 19. README 5: publication bias, the funnel
- [x] 20. README 6: validation, validation/reproduce.R and the expected values table
- [!] 21. tutorial pages still follow the old four-page structure (asked the operator; plots.png now unreferenced)
- [x] 22. devtools::test() and R CMD check to 0/0/0
- [ ] 23. commit and push

## Notes

The corner box moves for every curve, with or without a subgroup: the operator
asked for that first and separately.

Out of scope, not asked for: subgroups in `easydr_funnel()` and
`easydr_loo()`.
