# easyDR book, Chapter 4: statistics theory and code
Progress: 21/22 (95%)
Started: 2026-10-29 09:40 BRT   Updated: 2026-10-29 09:59 BRT

## Goal
`book/04-statistics.qmd` renders clean under Quarto, about 6,000 words, 14 H2
sections (one H2 = one slide), theory interleaved with runnable easyDR code.
Flow decided with the author: fundamental principles, the data the method
needs, what a dose curve is (within-study then between-study), the six shapes
with splines and knots in depth, an `Exploring heterogeneity` section holding
subgroups and leave-one-out, a `Publication bias` section holding the funnel
and Egger. Done means: Quarto renders, style_check passes, six new verified
references in references.bib, committed and pushed.

## Decisions (author, 2026-10-29)
- Replace the old 15-section skeleton. One-stage, robust error meta-regression
  and Bayesian survive as callouts, not H2 sections.
- Re-show all six shapes on one panel even though Chapter 1 showed four.
- Examples: `aripiprazole` (continuous) and `alcohol_cvd` (binary) as in
  Chapter 2; `alcohol_cancer` for the shape gallery and the funnel.
- Equations shown, then read aloud in prose, then the R output that the
  formula produces.
- Six references added, each verified before citing.

## Steps
- [x] 1. Verify environment: R 4.6.0, easyDR installed, five functions run (done in session)
- [x] 2. Probe every output the chapter quotes (continuous linear, rcs, subgroup, loo, egger)
- [x] 3. Verify the six new references (Durrleman 1989, Desquilbet 2010, Discacciati 2017, Harrell 2015, Sappani 2025, Jiang 2025)
- [x] 4. Add the verified entries to book/references.bib
- [x] 5. Draft H2 1 to 4: principles, the data, why estimates correlate, reconstructing the covariance (w=2)
- [x] 6. Draft H2 5 to 7: what a dose curve is, the first stage, the second stage (w=2)
- [x] 7. Draft H2 8: the six shapes fitted side by side
- [x] 8. Draft H2 9 to 10: the spline basis, and where the knots go (w=2)
- [x] 9. Draft H2 11: departure from linearity and goodness of fit
- [x] 10. Draft H2 12: Exploring heterogeneity (subgroups, leave one out)
- [x] 11. Draft H2 13: Publication bias (funnel, Egger)
- [x] 12. Draft H2 14: writing it up, with the reporting box
- [x] 13. Add the furniture: Test your knowledge, Summary bullets
- [x] 14. Render the chapter with quarto.exe and fix every chunk error (w=2)
- [x] 15. Read the rendered figures, confirm each H2 has one
- [x] 16. Run style_check.py, fix flags, re-run
- [x] 17. Word count check against the 6,000 target
- [ ] 18. Commit and push to origin on the easyDR repo

## Notes
Old skeleton named a physical-activity review (@garcia2023) that has no
packaged data, so it could never have rendered a live chunk. Dropped per the
chapter contract rule that an example is only usable if the source prints what
the method needs.
