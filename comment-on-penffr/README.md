# Comment on Tamo Tchomgui et al. (2026, *ADAC*) — materials

This folder contains a Comment on

> J. S. Tamo Tchomgui, J. Jacques, G. Fraysse, V. Barriac, S. Chretien (2026).
> *A penalized spline estimator for functional linear regression with functional response.*
> Advances in Data Analysis and Classification. DOI 10.1007/s11634-026-00681-w.
> Code: https://github.com/Orange-OpenSource/penffr

and the reproducible re-analysis supporting it.

## The three points

1. **Not novel.** "FFR/PenFFR" is penalized-spline (P-spline) function-on-function regression in a mixed-model representation — the method `pffr` (`refund`) and the FAMM framework have implemented since 2015. Only the conformal/optimal-transport prediction bands are arguably new, and they are orthogonal to the estimator.
2. **Under-cited.** The foundational references (Eilers & Marx 1996; Ruppert, Wand & Carroll 2003; Wood 2017 / Wood, Pya & Säfken 2016; Scheipl, Staicu & Greven 2015; Scheipl, Gertheiss & Greven 2016; Greven & Scheipl 2017; Malfait & Ramsay 2003; Currie, Durbán & Eilers 2006) are not cited.
3. **Unfair comparison.** `pffr` is mis-specified in the authors' code (scalar `lat`/`lon` entered as 50-dim time-varying functional terms; `bs="cr"` and `k`s that do not match the paper's Table 3; full-range `pffr` integral vs historical PenFFR). Correcting it removes much of `pffr`'s reported error; and the authors' own Hawaii table already ranks `pffr` first.

## Key reproduced numbers

| Experiment | Result |
|---|---|
| Canadian Weather, concurrent `pffr`, authors' code | ISE **89.5 (52.1)** — reproduces paper's 89.31 (52.03) |
| Canadian Weather, concurrent `pffr`, `lat`/`lon` fixed | ISE 64.8 (46.3) |
| Hawaii Ocean, concurrent `pffr` | ISE×10² **0.52** — exact match to paper; **best method in the paper's Table 5** |
| Simulation | penalized `pffr` is flat in basis dim `k`; error tracks noise (vs the paper's constant `pffr` MRPE) |

(See `comment.md` for the full argument and `results/` for tables/figures.)

## Reproducing

Requirements: R (≥4.1), `mgcv`. The scripts prefer an installed `refund`; if `refund`
is not installed they source `pffr` directly from a local `refund` checkout (set
`REFUND_R` to its `R/` directory). Datasets are fetched from the CRAN GitHub mirror
on first run (`fda::CanadianWeather`; `FRegSigCom`'s `ocean`).

```
Rscript analysis/00_setup.R          # loads pffr, fetches data, defines helpers
Rscript analysis/01_canadianweather.R
Rscript analysis/02_hawaii_ocean.R
Rscript analysis/03_simulation.R
```

## Contents

- `comment.md` — the Comment (submission-ready source); `comment.pdf` — typeset version.
- `cover-letter-editor.md` — cover note to the ADAC Editor-in-Chief.
- `analysis/` — reproduction scripts.
- `results/` — output tables and figures.
- `reference/` — extracted paper text, transcribed reported numbers, and the authors' own
  Canadian-Weather analysis script (`authors_RD_Canada.Rmd`).
- `vendor-penffr/` — the authors' released package sources (GPL-2.0-or-later, Orange SA),
  retained for verification of the method↔`pffr` correspondence.
