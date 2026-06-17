# Paper Plan & Handoff: Sandwich Confidence Intervals for `pffr`

Handoff document for a **methodological statistics paper** (target: JCGS / CSDA /
Biometrics-tier) built from the `pffr` sandwich-CI simulation work. Anyone
picking this up — Claude instance, agent, or co-author — should be able to read
this and know where everything is and what to do next. **Keep it current.**

## STATUS — 2026-06-17

- **Implementation:** merged to `master`, shipped in `refund` 0.1-40 (2026-03-17). Done.
- **Simulation studies (Study 1, Study 2):** run; results archived (see Data status). Done.
- **Report prose:** all TODO stubs in `pffr-ci-report.qmd` now drafted
  (Introduction/related-work, Methods/estimator definitions, both study
  Discussions, cross-study Discussion, Conclusions). Numbers in the Discussions
  are **provisional** — drafted from established findings (Part A4) and flagged
  `DRAFT/verify` in-source; must be regenerated from summaries.
- **Literature search:** complete; `pffr-ci-refs.bib` has 49 verified entries.
- **New experiments (Part B2):** plan written
  (`ci-benchmark/NEW-EXPERIMENTS-PLAN.md` on branch `ci-experiments`); coding not
  yet started. This is the current work front.
- **Council review (2026-06-17):** Claude + Codex legs returned (Gemini failed,
  empty). Findings + required fixes in section "Council review" below.
- **Real-data application (B3), theory (B4):** not started.

## Repository / branch / worktree map  ← READ FIRST

Everything lives under `~/fda/`. The work is split across branches/worktrees:

| Path | Branch | Holds |
|------|--------|-------|
| `~/fda/refund` (main checkout) | `claude/practical-wright-mljg9q` | **Paper prose**: `ci-benchmark/pffr-ci-report.qmd` + `pffr-ci-refs.bib` + this plan. The 3 prose commits (52e8ff4, 0e4105f, f2aaa4f) are here; pushed. |
| `~/fda/refund-worktrees/ci-paper` | `ci-paper` | Earlier paper worktree; has the **rendered** `pffr-ci-report.pdf` and the result-dir **symlinks**. |
| `~/fda/refund-worktrees/ci-experiments` | `ci-experiments` (off `pffr-refactor`) | **New-experiments dev**: production harness (`confint-benchmark.R`, `benchmark-utils.R`, `sim-study-*.R`, `run-production.sh`) + `NEW-EXPERIMENTS-PLAN.md`. |
| `~/fda/refund-worktrees/pffr-sandwich` | `pffr-sandwich` | Older sandwich-benchmark scripts (historical). |
| `pffr-refactor` branch | — | Canonical home of the production benchmark harness + merged package code (bootstrap, sandwich, simultaneous). |

**Note (to consolidate later):** the paper prose lives on a different branch than
the benchmark harness. The paper branch has the package code but not the
`sim-study-*.R` drivers; `ci-experiments` (off `pffr-refactor`) has both. New
result dirs produced by the experiments must be made loadable by the report
(`load_study1()`/`load_study2()` in `_pffr-ci-helpers.R`).

**Data status (ACTION NEEDED):** the heavy result directories
(`study1-nongaussian[-cl2]`, `study2-grid-refinement`, `study2-cl2-timing`,
`results`) are **git-ignored and currently absent** in the checkouts on this
machine — the `ci-paper` worktree symlinks to `~/fda/refund/ci-benchmark/<dir>`,
which do not resolve. They must be **located or regenerated** (harness on
`pffr-refactor`/`ci-experiments`) before the provisional report numbers can be
verified and before `quarto render` will run the analysis chunks. The
`study2-cl2-timing/` dir is needed for the `TODO(data)` computational-cost table.

---

## PART A — What we have

### A1. Implementation (merged to `master`, refund 0.1-40, 2026-03-17)
Three sandwich covariance estimators for `pffr`, clustering **by curve**:
- **`"hc"`** — observation-level Huber-White (heteroskedasticity only).
- **`"cluster"`** (CR1, **default**) — Liang–Zeger cluster-robust, clusters = curves.
- **`"cl2"`** — Bell–McCaffrey/CR2 leverage-adjusted cluster sandwich for small G.

Key code (`R/pffr-core.R`): `apply_sandwich_correction()`, `gam_sandwich_cluster()`
(CR1), `gam_sandwich_cluster_cl2()` (CL2), `assemble_cluster_sandwich()`,
`build_cluster_id()`, `compute_gaulss_scores()`, `build_cl2_working_*()`,
`sym_inv_sqrt()`. API: `sandwich=` in `pffr()`/`coef.pffr()`; `coef.pffr()` also
has `ci = c("none","pointwise","simultaneous")`, `level`, `n_sim`, `sim_seed`
(`R/pffr-methods.R`). Bootstrap: `pffr_coefboot()` (`R/pffr-robust.R`). Works for
Gaussian, `gaulss`, and non-Gaussian GLMs.

### A2. Simulation studies (results archived — see Data status)
- **Study 1 — non-Gaussian.** Poisson, Binomial × {IID, AR1(0.9), Fourier⁺(0.3)}
  × n∈{200,400} = 12 cells × 150 reps. Methods: default/HC/cluster/CL2. G=200–400.
  → `study1-nongaussian/` (+ `study1-nongaussian-cl2/`).
- **Study 2 — grid refinement.** Gaussian; 3 corr × n∈{20,40,80} × SNR∈{3,25} ×
  ny∈{40,80,120} (nx=60) = 54 cells × 50 reps. G=20–80.
  → `study2-grid-refinement/` (+ `study2-cl2-timing/`).
- Seeds: Study 1 base 3001, Study 2 base 4001; per-rep `base + 1000*dgp_id + rep_id`.

### A3. The report
`pffr-ci-report.qmd` → `pffr-ci-report.pdf` (JSS Quarto; helpers
`_pffr-ci-helpers.R`; refs `pffr-ci-refs.bib`). Figures render from the result
dirs. **Prose now drafted** (was all TODO). Superseded: `analysis.Rmd`.

### A4. Established findings (the Discussion draws on these; VERIFY when data is back)
1. Default Bayesian CIs collapse under within-curve correlation and **worsen on
   finer grids** (ff coverage ~62%→~43% under AR1) — headline novel result.
2. Cluster sandwich largely resistant (~86–87% across grids); HC insufficient.
3. CL2 G-dependent: negligible at G=200–400 (Study 1), +5–10pp at G=20–80 (Study 2).
4. Residual ~84–89% cluster-coverage gap — leading explanation penalization/
   smoothing bias (strongest evidence: persists at large G). NOT yet identified
   vs. alternatives (see council M5).
5. `sandwich::vcovCL(type="HC2")` broken for GAMs (unpenalized hat/bread) —
   motivates the custom penalized implementation.

---

## PART B — Remaining

### B1. Writing — mostly done; verification + council fixes outstanding
Prose drafted. Outstanding: (a) regenerate all provisional numbers from
summaries and remove `DRAFT/verify` flags; (b) fill the `TODO(data)`
computational-cost table from `study2-cl2-timing/`; (c) apply council fixes below.

### B2. New experiments / competitors — PLAN DONE, coding next
See `ci-benchmark/NEW-EXPERIMENTS-PLAN.md` (branch `ci-experiments`). Three
competitors, added via extension scripts that re-run the same seeds:
- **E1 — t_{G-1} critical values** for cluster/CL2 (cheap; local).
- **E2 — simultaneous bands** joint coverage via `coef(ci="simultaneous")` (local).
- **E3 — curve/cluster bootstrap** via `pffr_coefboot()` (heavy; **LRZ**).
Compute constraints: ≤3 local cores, no local job >2h; E3 on CoolMUC-4 via
`lrz-remote`. Open budget decisions: bootstrap B, resample method, CI type,
whether E1/E2 run on both studies.

### B3. Real-data application — NOT STARTED
Pick data with clear within-curve dependence where sandwich changes conclusions.
Candidates: DTI (refund), gasoline/octane, weather FoF. The council flags this as
important for the methods-paper framing (else narrow to a simulation/software paper).

### B4. Optional theory — NOT STARTED
Consistency sketch for the by-curve cluster sandwich; why CR2 helps at small G.

### B5. Reproducibility/housekeeping
Report sources tracked; results git-ignored/symlinked. Pin seeds/README; re-render
after numbers verified. Consolidate the branch split (B1/B2 on one branch).

---

## PART C — Proposed paper structure
1. Introduction (done) — under-studied inference in penalized FoF; contributions.
2. Model & default inference (done) — pffr model, Vp/Ve, Marra–Wood.
3. Robust covariance estimators (done) — HC, CR1, CR2/CL2; penalized-bread; t vs z.
4. Simulation study (ADEMP) — Study 1 + Study 2 + **B2 competitors** (in progress).
5. Application (B3) — real data. **Missing.**
6. Discussion (drafted) — when each method is needed; residual gap; limitations.
7. Software — the `sandwich=`/`ci=` API (short).
8. Appendix — estimator algebra, environment, extra figures.

---

## Council review — 2026-06-17 (Claude + Codex; Gemini failed)
Raw: `/tmp/council-ci-paper-codex.txt`; Claude leg in this session's transcript.
Required fixes before submission (convergent across both reviewers unless noted):

1. **Vp formula too literal** (both, HIGH). `Vp=(X'WX+S)^{-1}` / `Ve=Vp(X'WX)Vp`
   ignore reparameterization/constraints; `X'WX+S` not invertible in raw basis.
   Rephrase as "the posterior covariance returned by `mgcv`" or add generalized-
   inverse caveat. Also state `B2 = Vp − Ve = Vp S Vp ⪰ 0` (PSD because S⪰0) and
   call it the penalty-induced (Marra–Wood squared-bias) term, not "Bayesian
   correction" (Claude M2).
2. **Provisional numbers stated as findings** (both, HIGH). Keep figures as
   ranges and explicitly provisional until recomputed; the `t_{G-1}` attribution
   for the CL2 gain is unsupported until E1 is run — mark as conjecture.
3. **Over-strong claims** (both, MED). "immune" → "largely resistant"; the
   residual-gap attribution is not identified — enumerate alternatives (small-G
   reference dist, smoothing-parameter uncertainty, basis/discretization) and
   foreground the strongest evidence (gap persists at large G). Replace the
   "CL2 doesn't close it ⇒ bias" eliminative argument (Claude M5).
4. **Novelty/framing** (both). "no implementation available" → historical
   ("not available before this contribution"); add real-data app + competitors
   (B2/B3) or narrow framing to simulation/software paper.
5. **CR1 one-cluster guard** (Codex, code bug). `assemble_cluster_sandwich()`
   (`R/pffr-core.R:1137`) divides by `G-1` with no guard; CL2 checks `G<2` but
   CR1 doesn't → `Inf/NaN` for one cluster. Add a guard / state the ≥2-cluster
   requirement.
6. **Methods precision** (Claude). `gaulss` CL2 is a sign-split pseudo-obs
   factorization, not exact Bell–McCaffrey — don't claim exact CR2 for gaulss;
   `G/(G-1)` is a simplified CR1 factor, NOT "HC1"; mention the leverage cap
   (0.999) / eigenvalue floor; fix the `X_g` symbol collision (raw vs whitened);
   add Nychka (1988) and discuss `mgcv` `unconditional=TRUE` as a baseline.

---

## NEXT STEPS (ordered)
1. Apply council prose fixes 1–4, 6 to `pffr-ci-report.qmd` (+ bib: Nychka 1988).
2. Fix council #5 (CR1 one-cluster guard) in `R/pffr-core.R` + test.
3. Implement E1 (t-crit) and E2 (simultaneous) extension scripts on
   `ci-experiments`; smoke (1–2 reps, ≤3 cores) → R1 code-review gate → pilot.
4. Locate/regenerate the result data; verify & de-provisionalize all report
   numbers; fill the computational-cost table.
5. Size + launch E3 (bootstrap) on LRZ via `lrz-remote` (SLURM array per cell).
6. Decide + build the real-data application (B3).
7. Re-render the report; council R3 results review.

## Process notes
- Convene `council-of-bots` on each major draft and before production runs.
- Generate all reported numbers from code; never speculative.
- `air format` only on `.R`, never `.qmd`.
- Compute: ≤3 local cores, no local job >2h; heavy/long runs on LRZ.
</content>
