# New Experiments Plan (PAPER-PLAN Part B2): competitor CI methods

Adds the competitor methods named in `PAPER-PLAN.md` B2 to the existing two
studies, so the paper can compare **coverage / width / SD(z)** across a complete
method set and isolate the small-G reference-distribution effect. Developed on
branch `ci-experiments` (worktree off `pffr-refactor`, which carries both the
merged package code and the production harness).

## Design principle: extension scripts, not a re-run

We reuse the proven pattern of `sim-study-nongaussian-cl2-extension.R`: re-run the
**same seeds** as the production studies (so the simulated data + base fit are
byte-identical), and extract only the *new* competitor metric into a separate
output directory. This keeps each competitor additive, cheap to re-summarize, and
paired with the existing default/hc/cluster/cl2 results.

- Study 1 seeds: `STUDY1_BASE_SEED (3001) + 1000*dgp_id + rep_id`, `make_study1_settings()`.
- Study 2 seeds: base `4001 + 1000*dgp_id + rep_id`, `make_study2_settings()`.
- Per-(dgp,rep) atomic `.rds` saves + skip-existing resume (as in the cl2 ext).
- Modes: `smoke` (1–2 reps), `pilot` (10 reps), `full` (Study 1: 150, Study 2: 50).

## ADEMP

- **Aims.** Quantify how much of the default→cluster coverage story is changed by
  (i) a small-G reference distribution, (ii) simultaneous vs pointwise targets,
  and (iii) a fully nonparametric resampling alternative.
- **Data-generating mechanisms.** Unchanged — exactly Study 1 (non-Gaussian,
  G=200/400) and Study 2 (Gaussian grid-refinement, G=20/40/80). No new DGPs.
- **Estimands.** Unchanged: ff surface, linear/varying coefficient, intercept,
  E(Y); evaluated on the same `coef.pffr` grids.
- **Methods (new competitors).**
  - **E1 — t_{G-1} critical values** for `cluster` and `cl2`. Same fit, same SEs;
    replace z=qnorm(0.95) with qt(0.95, df=G-1), G = number of curves. Isolates
    the reference-distribution effect (expected negligible at Study 1 G≥200,
    material at Study 2 G≤80). Cheapest; reuses the base fit's covariance.
  - **E2 — simultaneous bands** via `coef(fit, sandwich=m, ci="simultaneous",
    n_sim=2000, sim_seed=...)` for m in {none(default), cluster, cl2}. Adds
    *joint* coverage (band covers truth at ALL grid points simultaneously) and
    mean simultaneous band width, alongside the existing pointwise numbers.
  - **E3 — curve/cluster bootstrap** via `pffr_coefboot(object, method="resample",
    B, conf=0.90, type="percent", parallel="multicore", ncpus=...)`. Percentile
    CIs per term → coverage/width. Fully nonparametric robustness baseline.
- **Performance measures.** Reuse `compute_term_metrics()` outputs (coverage,
  mean_width, z_mean, z_sd, rmse, bias, ...). Add `coverage_joint` and
  `width_sim` for E2. Always report Monte Carlo SEs (coverage:
  sqrt(p(1-p)/R)). Level = 90% throughout (existing `alpha = 0.10`).

## Integration points (verified in `confint-benchmark.R` / package)

- `extract_term_ci_df()` (confint-benchmark.R:755) hardcodes `z_crit <- qnorm(1 -
  alpha/2)` (line 818). **Add an optional `df = Inf` arg**: `crit <- if
  (is.finite(df)) qt(1 - alpha/2, df) else qnorm(1 - alpha/2)`. Backward
  compatible (default Inf → unchanged). Used by E1.
- `coef.pffr(..., ci="simultaneous")` returns `smterms[[k]]$coef$lower/$upper`
  (max-statistic band; pffr-methods.R:827-839, 1290-1334). E2 reads these
  directly instead of recomputing se±crit; joint coverage = `all(truth in
  [lower,upper])` over the term grid.
- `pffr_coefboot()` (R/pffr-robust.R:358) returns bootstrap CIs per coefficient;
  E3 maps its percentile CIs onto the same term grids.

## Compute budget and local-vs-LRZ split (constraints: ≤3 local cores, no local job >2h)

| Exp | Extra cost per (dgp,rep) | Studies | Where |
|-----|--------------------------|---------|-------|
| E1 t-crit | ~0 (one refit, reuse cov; no extra sandwich) | 1 + 2 | **local** (≤3 cores, minutes–<1h) |
| E2 simultaneous | refit + 3× (2000 MVN draws); cheap | 1 + 2 | local **pilot**; full likely local (<2h) or LRZ if Study 1 fits are slow |
| E3 bootstrap | refit + **B** refits per rep → dominant | 1 + 2 | **LRZ** (SLURM array over dgp cells), via `lrz-remote` |

E3 sizing (to decide): Study 1 = 12 cells × 150 reps = 1800 base fits × B; Study 2
= 54 cells × 50 = 2700 × B. With B≈199 and a ~few-second pffr fit, this is
~10^6 fits → SLURM array on CoolMUC-4, one task per dgp cell, `parallel=
"multicore"` within task. Local only for `smoke`/`pilot`.

## Workflow (setup-benchmark gates)

0. Parametrize `extract_term_ci_df()` with `df`; add E2 joint-coverage extractor.
1. `debug_single_fit`-style check: one (dgp,rep) for each of E1/E2/E3.
   **R1 code-review gate** (`/council-of-bots`) on the new extension scripts.
2. `smoke` (1–2 reps) for each, locally, in background.
3. `pilot` (10 reps) E1+E2 locally; E3 pilot on LRZ.
4. Timing → size the E3 SLURM array.
   **R2** review of pilot results.
5. Full runs: E1/E2 local; E3 on LRZ. Incremental saves + resume.
   **R3** results review before folding into the report.

## Open decisions (surfaced to user)

1. **Bootstrap budget B** and resample method (`resample` curve bootstrap vs
   `residual`/`residual.c`). Drives LRZ cost.
2. Run E1/E2 on **both** studies or Study 2 only (small-G is where t and boot
   matter most)?
3. Bootstrap CI type: `percent` (default) vs `bca` (costlier, better).

---

## E4 — fastFMM / FUI comparison (NEW repeated-measures sub-study)

**Finding (2026-06-17).** The first E4 driver
(`sim-study-fastfmm-extension.R`) tried to run `fastFMM::fui()` on the EXISTING
Study 1/2 DGPs. This is **infeasible**: those DGPs have ONE functional curve per
subject (independent curves), so the subject random intercept `(1 | id)` that FUI
requires is unidentifiable (n levels = n observations); lme4 rejects it at the
grouping-factor check *before* variance estimation, and `override_zero_var = TRUE`
does not help. fastFMM/FUI is a functional **mixed-model** method for
**longitudinal / repeated** functional data — its own `fui` example uses
`refund::DTI` with multiple visits per subject (`cca ~ visit + sex + (1|ID)`). So
fastFMM cannot be compared on the function-on-function, independent-curve setting
that is this paper's focus.

**Decision (user, 2026-06-17): add a repeated-measures sub-study where fastFMM is
valid (this section), reusing the fastFMM-paper DGP where possible.**

### Aims
Compare `pffr` (with a functional random intercept) vs `fastFMM::fui` for CI
coverage (pointwise + simultaneous) of the function-on-**scalar** coefficients
both can estimate, on a longitudinal functional DGP where FUI applies. This is the
setting fastFMM was built for; it tests whether `pffr`'s Bayesian/sandwich CIs
match FUI's massively-univariate-LME inference there, and complements the
function-on-function studies (Study 1/2).

### DGP — follow Cui, Leroux, Smirnova & Crainiceanu (2022, JCGS 31(1):219-230)
*First action: check whether the FUI paper's published simulation code (or the
package's `G_generate` / DTI structure) can be reused verbatim to set the random-
effect eigenfunctions and eigenvalues.* Target structure:
- N subjects, J visits each (e.g. N ∈ {50, 100}, J ∈ {3, 5}) → repeated curves per subject.
- `Y_ij(t) = β0(t) + x_ij β1(t) [+ z_ij β2(t)] + b_i(t) + ε_ij(t)`, t on a regular
  grid (L ≈ 50–100).
- `b_i(t) = Σ_{k=1}^{K} ξ_ik ψ_k(t)`, K ≈ 2–4 eigenfunctions with decreasing
  eigenvalues — the within-subject longitudinal correlation FUI models; `ε` white noise.
- Known truth `β0, β1 [, β2]` → coverage computable. Gaussian first; binomial
  variant (`fui(family="binomial", analytic=FALSE)` + bootstrap) as a stretch goal.
- Scope: function-on-**scalar** only (fastFMM cannot fit the `ff` surface). Compare
  the scalar-covariate coefficient function(s) and the functional intercept.

### Methods
- **pffr**: `pffr(Y ~ x + z + s(id, bs = "re"), ...)` (functional random intercept);
  default Bayesian, `cluster`, and `cl2` sandwich CIs (cluster = by subject `id`).
- **fastFMM**: `fui(Y ~ x + z + (1|id), data, analytic = TRUE)`; pointwise CIs from
  `betaHat_var` diagonal; simultaneous via the FUI max-statistic band (recompute
  `qn` at α = 0.10, as the current driver's `compute_qn_at_alpha()` already does).
- Same seeds / same data; paired comparison; 90% nominal; MC SEs over R reps.

### Estimands / measures
`β0(t), β1(t) [, β2(t)]` on the response grid. Pointwise coverage + mean width;
simultaneous joint coverage + band width (per Study 1/2 conventions).

### Compute / where
fastFMM analytic-Gaussian FUI is light (massively-univariate LMEs) → **local**
pilot OK. Binomial/bootstrap FUI and the full pffr×reps grid → **LRZ**. Driver:
new `sim-study-fastfmm-longitudinal.R` (or extend the existing driver with a
repeated-measures DGP path); reuse `extract_fui_metrics()` /
`compute_qn_at_alpha()` from the current driver.

### Install status (2026-06-17)
- fastFMM installs **locally** via PPM binaries (`fastFMM` + `Rfast` OK).
- **LRZ install unresolved**: `Rfast` fails to compile on R 4.3.3 / gcc13 —
  `Random.h` needs `#include <numeric>` for `std::iota`, plus a `LinkingTo`
  include-path issue (RcppArmadillo.h not found though installed). Tried: CRAN
  source (×2), `Ncpus=4` chain, GitHub dev (`RfastOfficial/Rfast`). Not needed for
  the light Gaussian-analytic pilot; resolve before any LRZ FUI runs (try a patched
  tarball with `<numeric>` added, or a spack/conda Rfast).
</content>
