# easyDR

Dose-response meta-analysis, made teachable: an R package and the book that
documents it.

- Book: <https://vanioljantunes.github.io/easyDR/>
- Issues and questions: <https://github.com/vanioljantunes/easyDR/issues>

## What this is

Most meta-analysis teaching stops at a single pooled estimate, which answers
"treated versus untreated" and nothing else. A dose-response meta-analysis asks
a different question: how does the effect change as the dose changes? The answer
is a curve, not a number.

The book teaches that method end to end, in seven chapters that follow the order
of the work: what it is, whether your studies can support one, registering it,
the theory, the code, the checklist, and how to write it up.

The package carries the two steps that students repeatedly get wrong, and leaves
the modelling to [dosresmeta](https://github.com/alecri/dosresmeta):

| Function | Does |
|---|---|
| `dr_viability()` | Checks whether a dataset can support a dose-response model at all, before any fitting |
| `curve_shape()` | Names the shape of a fitted curve in the vocabulary the literature uses |

## Install

```r
# install.packages("remotes")
remotes::install_github("vanioljantunes/easyDR")
```

## Build the book locally

```bash
quarto render book
```

The rendered site lands in `book/_book`. Rendering executes the R chunks, so the
package needs to be installed first:

```bash
R CMD INSTALL .
```

## Licence

MIT for the code, CC BY-NC-SA 4.0 for the book text. See `LICENSE` and
`LICENSE-text.md`.

## Related

Part of the MetaRank 4 catalog from the [Meta-Analysis
Academy](https://metahub.metaanalysisacademy.com): one design, taught end to
end, per repo. The book is free. The video course, templates and support are the
paid half.
