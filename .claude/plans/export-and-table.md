# easydr_table + file export on every plot

Progress: 20/20 (100%)
Started: 2026-10-29 14:05 BRT   Updated: 2026-10-29 15:20 BRT

## Goal

Two things, both reaching the README and the book:

1. A new exported `easydr_table()` that takes an `easydr()` fit, reports the
   effect at any doses the user names, and exports to xlsx or PDF.
2. Every plot function (`easydr_curve`, `easydr_forest`, `easydr_funnel`,
   `easydr_loo`, internal `easydr_bubble`) gains `filetype` (png default, pdf,
   tiff), `path` (working directory default), `file.name`, `width` and
   `height`. Width and height are worked out from the figure itself when the
   user does not give them, and the automatic size has to fit each plot type
   with no clipping and no dead space. Verified by eye on a batch of saved
   files, not by code reading.

Done means: the exported surface stays coherent, the constitution test is
updated, every saved example is inspected, README and the book chapters that
mention plots are updated, committed and pushed.

## Steps

- [x] 1. Shared saver in plot-utils.R: file-argument validation, natural-size
      measurement on a null device, device open and close
- [x] 2. `dr_draw_stack()` routes the five call sites through that one place
- [x] 3. `easydr_curve()` file args + roxygen
- [x] 4. `easydr_forest()` file args + roxygen (body_w and body_h are already
      known, so its natural size is exact)
- [x] 5. `easydr_funnel()` file args + roxygen
- [x] 6. `easydr_loo()` file args + roxygen
- [x] 7. `easydr_bubble()` file args (internal, same surface)
- [x] 8. `easydr_table()` in a new R/table.R: doses, main, digits, filetype
      xlsx or pdf, path, file.name (w=2)
- [x] 9. writexl to Suggests, roxygen run, NAMESPACE regenerated
- [x] 10. Constitution test updated for the new export; tests for the new args
      in test-plots.R (w=2)
- [x] 11. Save a batch of examples across every plot type and shape, read the
      PNGs, fix the sizing until each one fits (w=3)
- [x] 12. README section on saving figures and on easydr_table
- [x] 13. Book: the statistics chapter and chapter 5 Table 2 use
      easydr_table; mention the save arguments where plots appear
- [x] 14. Full test suite, commit, push

## Notes

Weights: steps 8, 10 count 2 and step 11 counts 3, so the total is 20.

The package print method emits a 186 to 251 character "Defaults applied" line
that overflows the book's content column. A `strwrap()` fix works but
`test-invariance.R` records the printed lines verbatim, so it would mean
re-recording seven captured calls. Left alone and reported.

Chapter 5 carries no visible code by the author's instruction, so the saving
arguments are described there in prose and shown as code in Chapter 4 and the
README.

Another session has been committing to this repo. Fetch before each commit.
