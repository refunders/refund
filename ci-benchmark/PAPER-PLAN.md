# Paper Plan & Handoff: Sandwich Confidence Intervals for `pffr`

Handoff document for a **methodological statistics paper** (target: JCGS / CSDA /
Biometrics-tier) built from the `pffr` sandwich-CI simulation work. Anyone
picking this up — Claude instance, agent, or co-author — should be able to read
this and know where everything is and what to do next. **Keep it current.**

## STATUS — 2026-06-17

- **Implementation:** merged to `master`, shipped in `refund` 0.1-40 (2026-03-17). Done.
- **Simulation studies (Study 1, Study 2):** run; results archived (see Data status). Done.
- **Report prose:** all TODO stubs in `pffr-ci-report.qmd` drafted; **council
  peer-review fixes applied** (commit 09791ae8: Vp/B2 precision, de-provisionalized
  numbers, "immune"→"largely resistant", gap-attribution rewrite, gaulss-CL2/HC1/
  leverage-cap caveats, Nychka 1988 + `unconditional=TRUE`, CR1 one-cluster guard
  in `R/pffr-core.R`). Numbers in the Discussions are still **provisional ranges**
  flagged `DRAFT/verify` — now regenerable since the data is local (see below).
- **Rendering:** `quarto render pffr-ci-report.qmd --to jss-pdf` **works** (18-page
  PDF, visually verified). Was broken by the `.gitignore` `*.tex` rule swallowing
  the `_extensions/.../partials/*.tex` JSS/ACM templates so they were never
  committed — fixed in d58bd1e9 (restored partials + `!_extensions/**/*.tex`).
  **Outstanding content gap:** the **abstract and keywords are still placeholders**
  ("—!!!—...required—!!!—") and must be written.
- **Literature search:** complete; `pffr-ci-refs.bib` has 49 verified entries.
- **New experiments (Part B2):** foundation + drivers written on `ci-experiments`
  (foundation b323db6; drivers 87f48652). Decisions: **E3 = B=499, percentile,
  curve resample, both studies (LRZ)**; B3 deferred.
  **R1 council code-review done (2026-06-17, all 3 legs):** E1 (t_{G-1}) and E2
  (simultaneous) drivers validated — Study 2 seed-fidelity correct, G=n/df=G-1 on
  cluster+cl2 only, joint/pointwise coverage correct (confirmed by reading smoke
  output). **3 driver bugs found, fixer in progress:** (1) E3 bootstrap produces
  all-NA rows (extractor fails silently); (2) Study 1 competitors driver
  `run_competitors1_rep` has a truncated/incomplete simulate+fit call;
  (3) `boot-extension.slurm` lacks `set -euo pipefail`, doesn't ensure `logs/`,
  and its Study 2 failure markers are silently dropped by the loader.
- **Council review (2026-06-17):** Claude + Codex legs returned (Gemini failed,
  empty). Findings + required fixes in section "Council review" below.
- **Real-data application (B3):** deferred by decision. **Theory (B4):** not started.

## SESSION HANDOFF — 2026-06-17 (afternoon, RESUME HERE)

**E3 boot bug ROOT-CAUSED & fixed; LRZ jobs submitted; full council done + fixes applied.**

**Done this session:**
- **E3 boot all-NA root cause found** (Codex): the stored model call captured
  `family = fam` (a *local* var in the driver); when `boot::boot()` re-ran the call
  in another frame, `fam` was gone → every refit errored → all-NA. Fix in
  `prepare_modcall_for_bootstrap()` carries `object$family` forward (commit
  `10e0241b` on `ci-experiments`, + regression test). The earlier `refund::pffr`
  qualification (`3af325d6`) was necessary but not sufficient. **Validated on LRZ**:
  smoke now gives non-NA coverage, all replicates succeed.
- **E3 + cl2-timing submitted on LRZ** (cluster serial, partition serial_long/std):
  study1 boot `5267685_[1-12]`, study2 boot `5267686_[1-18]`, cl2-timing `5267687`.
  (Cancelled duplicate older arrays `5266858`/`5266859`.) Fetch when done:
  `rsync -az lrz:refund/ci-benchmark/study{1,2}-boot/ ci-benchmark/study{1,2}-boot/`
  and `study2-cl2-timing/`.
- **Full council review (all 3 legs: Codex, Gemini, Claude)** run on the revised
  paper; findings verified against code and applied (see "Council round 2" below).
  Headline: nominal-level error `z=1.96`→`1.645` (code was right, prose wrong);
  CL2 `z_{gd}` prose factorization corrected; Fourier⁺ PSD-projection note; Study 2
  grid story (`nx=60` vs `90×120` truth grid) clarified; Vp/CR1/HC notation, leverage
  symmetry, `sandwich`/`V_c` citations; duplicated DRAFT blocks removed.
- **Abstract + keywords written** (were placeholders). Keywords must be PLAIN text
  (no `\pkg`/`\proglang` — they leak into `pdfkeywords`/`\Plainkeywords` and break
  LaTeX).

**Branches:** paper prose → `claude/practical-wright-mljg9q`; experiments →
`ci-experiments` (boot fix `10e0241b`).

**IN FLIGHT / BLOCKED:**
- **E4/fastFMM — REFRAMED.** The driver (`sim-study-fastfmm-extension.R`) is
  committed (`c5c7703e`), but `fui()` **cannot fit the existing Study 1/2 DGPs**:
  one curve per subject → degenerate `(1|id)` random effect (lme4 rejects it).
  fastFMM/FUI needs **repeated/longitudinal** functional data. **Decision (user):
  add a repeated-measures function-on-scalar sub-study following Cui et al. (2022,
  JCGS) where FUI is valid — PLANNED in `NEW-EXPERIMENTS-PLAN.md` (E4 section),
  not yet implemented.** fastFMM installs **locally** (PPM binaries); **LRZ Rfast +
  fastFMM install RESOLVED by Codex (2026-06-17)** — removed the `~/.Rprofile`
  stdout banner (was breaking include-path detection), added `<numeric>` to Rfast
  sources, and `PKG_CPPFLAGS =`→`+=` in Rfast `Makevars`. Both load on LRZ now;
  details in `NEW-EXPERIMENTS-PLAN.md` (E4 install status).
- **LRZ E3 boot jobs**: `5267685_[1-12]` (study1), `5267686_[1-18]` (study2)
  queued, run over hours/days. Monitor `squeue`; fetch when done.
- **`study2-cl2-timing/` + `study2-results-extracted/` RECOVERED (2026-06-17)** from
  the "Teleport auto-stash" (`git stash@{0}`, main tree — teleport had `git add`-ed
  everything) via `git restore --source=stash@{0} -- <paths>`; now present in
  `~/fda/refund/ci-benchmark/` and the ci-paper symlinks resolve. The redundant LRZ
  cl2-timing job (`5267687`) was **cancelled**. → The `TODO(data)` cost table is now
  fillable from `study2-cl2-timing/summary_by_method.csv` (+ `summary_overheads.csv`).
  Keep the stash (it holds the full pre-teleport working state, 19,498 files).

**TO RESUME (ordered):**
1. Confirm paper renders clean (PDF) after this session's edits; commit qmd.
2. Fetch LRZ results (E3 boot, cl2-timing) when jobs finish; integrate into report
   (`load_*` helpers) and de-provisionalize numbers; fill `TODO(data)` cost table.
3. Implement the planned E4 **repeated-measures** sub-study (see
   `NEW-EXPERIMENTS-PLAN.md` E4 section): new longitudinal function-on-scalar DGP
   per Cui et al. (2022) where `fui` is valid; pilot locally (fastFMM installs via
   PPM binaries). Resolve LRZ `Rfast` only if LRZ FUI runs are needed.
4. Optional council polish: gaulss `logb` detail (M6), leverage cap/floor wording
   (M7), index `g`/`i` & `n_g`/`D_g` unification (M12/M13), M10 ff-identifiability
   (`∫β ds=0` vs per-t centering) — verify against DGP truth-centering.
5. Pilot E1/E2 locally (≤3 cores) or on LRZ.

---

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

**Data status:** the production result data was on **LRZ** (`~/refund/ci-benchmark/`)
and has been **rsync'd local** into `~/fda/refund/ci-benchmark/`:
`study1-nongaussian/` (1801 rds), `study1-nongaussian-cl2/` (1801),
`study2-grid-refinement/` (2702; per-rep in `main/`, plus `main_results_combined.rds`
and `cov_quality.rds`). The `ci-paper` worktree symlinks now **resolve**, so
`load_study1()`/`load_study2()` and `quarto render` work. STILL MISSING (not on
LRZ): `study2-cl2-timing/` and `study2-results-extracted/` — the
`study2-cl2-timing/` run must be **regenerated** for the `TODO(data)`
computational-cost table. Refetch command: `rsync -az lrz:refund/ci-benchmark/<dir>/
ci-benchmark/<dir>/`.

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
- **E4 — `fastFMM` (Fast Univariate Inference, FUI; Cui et al. 2022, JCGS
  31(1):219–230)** as an external competitor. SCOPED: FUI fits function-on-**scalar**
  / concurrent functional mixed models via massively-univariate GLMMs + smoothing +
  analytic (Gaussian) or subject-bootstrap (non-Gaussian) joint bands. It
  **structurally CANNOT estimate our `ff` function-on-function surface** (stated
  limitation in the paper) → compare head-to-head only on **intercept, linear
  (scalar), and concurrent** terms; the ff inapplicability is itself a finding that
  motivates pffr. FUI assumes within-curve independence, so under AR1/Fourier it
  should be miscalibrated (a point in our favour). Wrinkle: our DGP has one curve =
  one unit (no repeated measures), so the `(1|id)` random effect is degenerate — the
  E4 adapter must verify `fui()` runs (with/without RE) on our data. Cost: Gaussian
  analytic is cheap (local); non-Gaussian needs boot=500 (LRZ). Citation added to
  `pffr-ci-refs.bib` (`CuiLerouxSmirnovaCrainiceanu2022`, `Rpkg_fastFMM`).
  Frame in paper as a **partial/related-work competitor** on the scalar terms.
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

## Author feedback — round 2 (2026-06-17)
1. **Notation:** boldface ALL matrix/vector quantities consistently. DONE (W3 sweep captured)
   (W3 global sweep; W1/W2 did their own sections).
2. **Title:** "Confidence Intervals for Function-on-Function Regression". DONE (aae08bf4).
3. **§2.1 (model):** generalized GLM rewrite (link + random effects, Gaussian-only
   residuals, Greven & Scheipl 2017 notation, explicit conditional-independence-
   on-the-additive-predictor assumption). DONE (aae08bf4).
4. **§2.2 (sandwich):** expanded methodological core incl. detailed gaulss score
   block, CR1/CL2 derivation, leverage cap, penalized-bread. DONE (aae08bf4).
5. **§3 + §4 → one Results section** with ADEMP setup intro. DONE (36b3a50b):
   `## Simulation studies {#sec-results}` with Study 1/2 as `###` subsections,
   anchors preserved, no broken refs; renders 22 pp. Then W3 boldface sweep + council.

## NEXT STEPS (ordered)
1. ~~Apply council prose fixes + Nychka~~ DONE (09791ae8).
2. ~~CR1 one-cluster guard~~ DONE (09791ae8).
3. ~~Fetch production result data~~ DONE (rsync from LRZ; symlinks resolve).
4. E1/E2 driver scripts on `ci-experiments` (in progress, agent) → R1 council
   code-review gate → pilot.
5. E3 bootstrap driver + LRZ SLURM (in progress, agent) → submit on LRZ.
6. **Verify & de-provisionalize report numbers** now that data is local: render
   the analysis chunks, compare computed summaries to the provisional ranges,
   update prose, remove `DRAFT/verify` flags.
7. Regenerate `study2-cl2-timing/` → fill the `TODO(data)` cost table.
8. Real-data application (B3, deferred). Re-render; council R3 results review.

## Process notes
- Convene `council-of-bots` on each major draft and before production runs.
- Generate all reported numbers from code; never speculative.
- `air format` only on `.R`, never `.qmd`.
- Compute: ≤3 local cores, no local job >2h; heavy/long runs on LRZ.
</content>
