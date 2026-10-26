# easyDR round 3: clearer forest columns, I2, metareg hidden, humanized text

Progress: 25/25 weighted (100%)
Started: 2026-10-26 BRT   Updated: 2026-10-26 BRT

## Goal

Five corrections. The forest says what its columns count. I-squared is written
I2. The curve carries its own result in the corner. The meta-regression and its
bubble plot leave the public surface, because the source this package teaches
from does not present them. And every word a reader sees is rewritten against
the humanizer patterns.

## Steps

### A - the Levels column
- [x] 1. header becomes Dose levels, and Dose range carries its unit
- [x] 2. a legend under the forest says what each column counts
- [x] 3. test pins the headers

### B - I2
- [x] 4. forest writes I2, not I-squared
- [x] 5. test pins it

### C - the curve carries its result
- [x] 6. annotation box top left: effect, CI, p value, I2
- [x] 7. the box is the fit's own numbers, not recomputed; test pins that
- [x] 8. a curve has no single effect, so the box names the dose it used

### D - hide the meta-regression and the bubble
- [x] 9. easydr_metareg and easydr_bubble unexported, code and tests kept
- [x] 10. package page, NAMESPACE and the export-surface test follow
- [x] 11. plots.html drops its section 3; overview drops two rows
- [x] 12. figures.R stops building fig-bubble; README chunk follows

### E - humanizer pass
- [x] 13. README prose (w=2)
- [x] 14. overview.html and data.html prose (w=2)
- [x] 15. api.html and plots.html prose (w=2)
- [x] 16. roxygen description and section blocks (w=2)

### Ship
- [x] 17. figures and pages rebuilt from one run
- [x] 18. devtools::test 0 failures
- [x] 19. R CMD check 0/0/0
- [x] 20. read every rendered page back
- [x] 21. commit and push

## What went wrong on the way

Truncating edit: replacing a roxygen block by slicing from the file start to
the function name deleted the five dr_mod_* helpers above easydr_metareg().
The suite caught it at once, 10 failures and 13 errors, and git restored the
file. The slice now runs from the blank line before the block, not from zero.

ggplot2 4.x dropped label.size on annotate("label"), so the curve warned on
every draw. A test asserting the call is silent found it.

dr_p() returns "< 0.001", so every caller writing "p = " in front produced
"p = < 0.001". dr_p_lab() now writes the operator itself.

## Humanizer patterns to hunt in my own text

1  not-X-but-Y            used constantly: "a trend, not an effect"
2  one-line closers       "which is the point of it", "that is the usual mistake"
4  staged run-ups         "A dose-response review goes wrong twice."
6  forced triads          three-item lists with two real members
7  repeated openings      runs of sentences opening "The" or "A curve"
9  stacked qualifiers     "almost never", "quietly"
15 -ing riders            ", making it", ", allowing"
18 replacing is/are       "serves as", "provides"
26 wrong reader           re-explaining what the page just said

Weights: steps 13 to 16 are w=2. Total weight 25.
