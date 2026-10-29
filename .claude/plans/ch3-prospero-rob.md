# easyDR book, Chapter 3: PROSPERO and risk of bias
Progress: 8/10 (80%)
Started: 2026-10-29 08:14 BRT   Updated: 2026-10-29 09:02 BRT

## Goal
`book/03-prospero.qmd` written in full, about 3,000 words, six content H2s that
match the class deck one to one, checked by `style_check.py`, rendering clean,
and the deck outline patched so book and class agree. Done means: style check
passes, `quarto render` succeeds, both repos pushed.

## Decisions taken with the author (2026-10-29)
1. Six content H2s; the deck gains two slides rather than the chapter losing
   nine sections.
2. `ROBINS-II` on deck slide 5 is wrong and becomes `RoB 2`, in the book and in
   the deck.
3. The statistical analysis plan is carried by prose and quoted protocol
   blocks. No decision table, no flowchart.
4. No screenshots and no plots. Text furniture only: callouts, quoted protocol
   paragraphs, and `kable` tables.
5. Two worked examples, the same two as Chapter 2: `aripiprazole` (continuous,
   five randomized trials, RoB 2) and `alcohol_cvd` (binary, four cohorts,
   ROBINS-E plus JBI for the case-control studies).

## Steps
- [x] 1. Read the chapter contract, chapters 1 and 2, and the deck outline
- [x] 2. Ground the statistics section: easyDR README and Crippa and Orsini 2016 via docling
- [x] 3. Grill the author on structure, tool naming, stats-plan form, worked examples
- [x] 4. Draft 03-prospero.qmd, six H2s plus summary and quiz (w=3)
- [x] 5. style_check.py on the draft, analysis mode
- [x] 6. style_check.py --setup on the guidance sections
- [x] 7. quarto render, confirm the chapter builds and the two kables print
- [x] 8. Patch the deck outline: RoB 2, and the two new slides
- [ ] 9. Commit and push easyDR
- [ ] 10. Commit and push claudeOS for the deck change

## Notes
The deck lives at `claudeOS/metahub/slides/metarank4/dose-response/dose-response/`,
titled "4.1.3 Prospero and Risk of Bias", and is a different folder from
`4.1 Dose-response`, which is the Chapter 2 deck.
