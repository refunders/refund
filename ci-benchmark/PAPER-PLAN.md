# Paper Plan: Sandwich Confidence Intervals for `pffr`

Working plan for turning the existing `pffr` sandwich-CI simulation report into a
**methodological statistics paper** (target: JCGS / CSDA / Biometrics-tier).

## Context

In Feb–Mar 2026 we developed and merged cluster-robust ("sandwich") confidence
intervals for `pffr()` (penalized function-on-function regression) into `refund`
(shipped in 0.1-40). Two simulation studies are run and a JSS-format report
skeleton exists with all figures rendered but **all prose sections are empty TODO
stubs**. Lifting this to a methods paper requires: sharper novelty framing, a
broader set of competitor methods, ideally a theoretical justification, and a
real-data application.

This paper is developed in the `ci-paper` git worktree. The heavy result
directories (~350 MB) are symlinked to the main checkout and git-ignored; only
the report sources and this plan are tracked.

---

## PART A — What we already have

### A1. Implementation (merged to `master`, refund 0.1-40, 2026-03-17)
Three sandwich covariance estimators for `pffr`, clustering **by curve**:
- **`"hc"`** — observation-level Huber-White (heteroskedasticity only).
- **`"cluster"`** (CR1, now the **default**) — Liang–Zeger cluster-robust,
  clusters = curves; handles within-curve heteroskedasticity + autocorrelation.
- **`"cl2"`** — Bell–McCaffrey/CR2 leverage-adjusted cluster sandwich for small G.

Key code (`R/pffr-core.R` unless noted): `apply_sandwich_correction()` (dispatcher),
`gam_sandwich_cluster()` (CR1), `gam_sandwich_cluster_cl2()` (CL2),
`build_cluster_id()`, `compute_gaulss_scores()`, `assemble_cluster_sandwich()`,
`build_cl2_working_standard/gaulss()`, `sym_inv_sqrt()`. API: `sandwich=` in
`pffr()` and `coef.pffr()`; `coef.pffr()` also gained
`ci = c("none","pointwise","simultaneous")`, `level`, `n_sim`, `sim_seed`
(`R/pffr-methods.R`). Works for Gaussian, `gaulss`, and non-Gaussian GLMs.
`pffr_gls()`/`pffrGLS()` deprecated → error. Tested in `tests/testthat/test-pffr.R`.
Provenance: `c89ab52` (HC) → `4bfc493` (CR1) → `82e1931` (gaulss/GLM) →
`3548a86` (CL2) → mega-refactor `c6541ec` (merged 2026-03-17).

### A2. Simulation studies (complete; results archived, symlinked)
- **Study 1 — non-Gaussian families.** 2 families (Poisson, Binomial) × 3
  correlations (IID, AR1(0.9), Fourier⁺(0.3)) × 2 sizes (n=200,400) = 12 cells ×
  150 reps. Methods: default/HC/cluster/CL2. Terms: ff, linear, intercept.
  → `study1-nongaussian/` (+ `study1-nongaussian-cl2/`).
- **Study 2 — grid refinement.** Gaussian; 3 corr × 3 n (20,40,80) × 2 SNR
  (3,25) × 3 response grids (ny=40,80,120, nx=60) = 54 cells × 50 reps.
  → `study2-grid-refinement/` (+ `main_summary.csv`, `cov_quality.rds`,
  `study2-cl2-timing/`). Older Gaussian round in `results/` (analysis.Rmd).

### A3. The report
- **`pffr-ci-report.qmd`** → rendered **`pffr-ci-report.pdf`** (JSS Quarto,
  helpers `_pffr-ci-helpers.R`, refs `pffr-ci-refs.bib`). Figures + analysis code
  all present and rendered; **every prose section is a TODO stub**.
- Superseded: `analysis.Rmd` / `analysis.html` (earlier combined analysis).

### A4. Established findings
1. Default Bayesian CIs collapse under within-curve correlation and **worsen on
   finer response grids** (ff coverage 62%→43% under AR1) — headline novel result.
2. Cluster sandwich is immune (stable ~86–87% across grids); HC insufficient.
3. **CL2 is G-dependent**: negligible at G=200–400 (Study 1), +5–10pp over cluster
   at G=20–80 (Study 2).
4. Residual ~84–89% cluster-coverage gap — penalization/smoothing bias, not
   finite-sample sandwich bias (CL2 doesn't close it).
5. `sandwich::vcovCL(type="HC2")` is broken for GAMs (unpenalized hat matrix vs.
   penalized bread) — motivates our custom penalized implementation.

---

## PART B — What remains to be done

### B1. Writing
Fill every TODO stub in `pffr-ci-report.qmd`, generating prose **from computed
summaries** (`results: asis` + `cat(sprintf(...))`; never speculative). Add:
model/notation, precise estimator definitions (CR1/CL2 in the penalized-GAM
setting, with the `B2 = Vp − Ve` Bayesian term), grid-collapse mechanism,
coverage-gap discussion, conclusions.

### B2. New experiments / comparisons (to lift to a methods paper)
Add as competitors; re-summarize coverage / width / SD(z):
- **Curve/cluster bootstrap CIs** via `pffr_coefboot()`.
- **Simultaneous bands**: joint coverage (`ci="simultaneous"`) alongside the
  current pointwise measure — for default, cluster, CL2, boot.
- **t_{G−1} critical values** for cluster/CL2 (vs. z=1.96) to isolate the small-G
  reference-distribution effect.
- **Literature scan** for further competitors to cite/compare: Goldsmith–Greven–
  Crainiceanu (2013) corrected PC-based bands; Krivobokova–Kneib–Claeskens (2010)
  and Degras (2011) simultaneous bands; Choi–Reimherr (2018) joint regions;
  Greven–Scheipl (2017) FAMM; wild bootstrap for FDA.

### B3. Real-data application
Pick a dataset with clear within-curve dependence so the sandwich visibly
changes conclusions vs. defaults. Candidates: **DTI** (refund), gasoline/octane,
a weather/FoF example, or our own data. Decision deferred to first work session
(depends on B2 lit scan).

### B4. Optional theory
Sketch that the by-curve cluster sandwich is consistent for the sampling variance
of penalized FoF coefficients under within-curve dependence, and why CR2 helps at
small G. Scope decision deferred.

### B5. Reproducibility/housekeeping
Report sources tracked; heavy results git-ignored (symlinked). Pin a seed/README
so figures regenerate; re-render after prose is added.

---

## PART C — Proposed paper structure

1. **Introduction** — under-studied inference in penalized FoF regression; default
   mgcv Bayesian CIs assume correct error structure; within-curve dependence is the
   FDA norm. Contributions: (i) cluster/CL2 sandwich for penalized FoF in `refund`;
   (ii) the grid-refinement SE-collapse phenomenon; (iii) practical guidance.
2. **Model & default inference** — the `pffr` model, penalized-GAM representation,
   Bayesian posterior covariance `Vp=(X'WX+S)^{-1}`, Marra–Wood coverage argument.
3. **Robust covariance estimators** — HC, CR1 cluster (by curve), CR2/CL2; the
   penalized-bread subtlety (why `sandwich::vcovCL` fails); pointwise vs.
   simultaneous; t vs. z reference. (Optional consistency sketch — B4.)
4. **Simulation study** (ADEMP) — Study 1 (non-Gaussian) + Study 2 (grid
   refinement), now including bootstrap + simultaneous + t_{G−1} competitors.
   Headline figure: grid-dependent default collapse. Coverage-gap analysis.
5. **Application** (B3) — real data; CI/conclusion changes default→cluster→CL2.
6. **Discussion** — when each method is needed; the residual gap; limitations
   (penalization bias, cluster-count dependence, computational cost).
7. **Software** — the `sandwich=`/`ci=` API in `refund` (short subsection).
8. **Appendix** — estimator algebra, computational environment, extra figures.

---

## PART D — Worktree setup (DONE)

- Worktree `ci-paper` created off `master` at `~/fda/refund-worktrees/ci-paper`.
- Report sources + Quarto JSS `_extensions/` copied in; heavy result dirs
  symlinked to the main checkout and git-ignored.
- Verified: `load_study1()`/`load_study2()` resolve through symlinks (21600 /
  64800 rows); `quarto render pffr-ci-report.qmd --to jss-pdf` succeeds.

## Process notes
- Convene the **council-of-bots** on scope/structure and each major prose draft,
  iteratively.
- Generate all reported numbers from code; no speculative findings.
- Only `.R` files get `air format`; never the `.qmd`.
