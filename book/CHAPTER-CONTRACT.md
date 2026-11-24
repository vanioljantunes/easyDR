# How a chapter of this book is built

Agreed with the author on 2026-10-28, while planning Chapter 2. It holds for
the remaining chapters of this book and for the easyPR and easyGWA books,
until the author changes it here.

## Voice

The `doing-meta-style` skill, installed at
`metahub/.claude/skills/doing-meta-style/`. Invoked before drafting, not
after. Checked before delivery:

```
python <skill>/scripts/style_check.py book/0N-chapter.qmd
```

Analysis prose carries the inclusive **we**. Setup and guidance sections carry
**you**, and are checked with `--setup`.

## Structure

Changed by the author on 2026-11-24: **one H2 is one item of the deck's
agenda, and one H3 is one slide inside it.** In 4.1.4 and 4.1.5 the agenda
items are the block eyebrows. The H2 titles show in the left column of the
book, under the chapter list, so the reader sees the class agenda while
reading. Finer points live as H4 subsections or callouts inside a slide's H3.

A chapter starts from its deck (`deck.html`, what the class shows), with
`outline.md` as the brief. The deck is no longer regenerated from the
chapter: `scripts/chapter-to-outline.mjs` still reads H2 as a slide and would
now produce one slide per block, so do not run it on these chapters.

## The step-by-step

Steps are **interleaved with the examples**, not stated once and demonstrated
later. Each step is named, then immediately shown on every worked example
before the next step begins. The reader sees the same step done on continuous
data and on binary data side by side.

## Worked examples

Two per chapter where the data allows it: **one continuous, one binary**. Both
come from the meta-analyses the matching class used, so the chapter and the
class teach off the same papers. An example the class never showed is not a
substitute; if its data is missing, the data gets extracted and packaged.

Where the class made an extraction screenshot, the chapter uses the data in
that screenshot, so the reader who watched the class meets the same rows.

An example is only usable if the source prints what the method needs. Chapter
2 lost Furukawa 2019 this way: the article pools 77 studies, and its
supplementary appendix gives each arm a dose and a denominator but never a
responder count, so no binary extraction is possible from any published
source. The appendix is kept at `ma-examples/Furukawa-2019-appendix.pdf` as
evidence. Check this before promising an example, not after.

## Figures

**Live R output, produced by runnable chunks in the chapter.** Tables, plots
and printed objects are generated at render time so they cannot go stale and
the reader reproduces them.

The one exception is a tool that lives outside R. A spreadsheet extraction
template is shown as a screenshot, because a picture is the honest way to show
a spreadsheet.

Figures cropped from a published paper are allowed only when the point is what
that paper printed.

## Columns

The book teaches the column names of the class workbook
(`book/data/easydr-data.xlsx`, the same file the 4.1.4 class uses). easyDR
reads them as they are.

| Outcome | Columns |
|---|---|
| Binary | `studlab`, `type`, `dose`, `event`, `n`, `TE`, `lower`, `upper` |
| Continuous | `studlab`, `dose`, `mean`, `sd`, `n` |

`type` is `cc` (case-control), `ci` (cumulative incidence) or `ir`
(incidence rate). There is no reference column: the reference row is the one
with `TE` = 1 and empty limits, and easyDR finds it. Any extra column (such as
`nordic`) holds a subgroup. The packaged data sets (`alcohol_cvd`,
`aripiprazole`) keep their older names (`event.e`, `n.e`, `design`, `ref`),
which easyDR also accepts.

## Chapter boundaries

A chapter does not spend the next chapter's material. Chapter 2 fits a model
only far enough to report that it **converged and produced a finite estimate
with a usable interval**. It does not name a shape, test non-linearity, or
read a slope aloud. That is Chapter 4.

## Furniture

A chapter is as long as its material. There is no word target. There was one
until 2026-10-29, "about 3,000 words", and it was dropped because it was a bad
rule: a chapter that walks a form field by field is longer than one that fits a
model, and padding or cutting to reach a number serves nobody.

Every chapter ends with:

- a **Test your knowledge** box, four questions, answers in the appendix
- **Summary** bullets, one per key idea, no new information
- a **reporting box** where the chapter produces something that goes into a
  paper: the methods sentence, quoted, as the reader would write it

Callout boxes follow the style's types: a package and data box when a data set
first appears, an important box for a limitation or a misuse.
