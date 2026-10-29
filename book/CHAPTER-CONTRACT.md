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

**One H2 is one slide, and the chapter's H2 sections are the deck's content
slides, one to one.** The deck is regenerated from the chapter, so the two
cannot drift. Finer points live as H3 subsections or callouts inside an H2,
never as an H2 of their own.

A chapter therefore starts from its deck's `outline.md`, not from a blank
page.

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

The book teaches **easyDR's own column names**, which are the names the
packaged data sets already carry.

| Outcome | Columns |
|---|---|
| Binary | `studlab`, `dose`, `TE`, `lower`, `upper`, `event.e`, `n.e`, `design`, `ref` |
| Continuous | `studlab`, `dose`, `mean.e`, `sd.e`, `n.e`, `ref` |

`ref` is `TRUE` on each study's reference level. `design` is one of
`case-control`, `incidence rate`, `cumulative incidence`.

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
