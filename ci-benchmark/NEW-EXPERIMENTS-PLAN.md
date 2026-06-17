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
</content>
