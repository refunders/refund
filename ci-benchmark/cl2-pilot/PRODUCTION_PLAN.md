# CL2 Production Benchmark Plan (ADEMP)

Date: 2026-02-18  
Scope: benchmark-only planning and execution protocol for deciding whether to promote CL2 into core `pffr` API.

## 1) Aims

Primary aim:

1. Quantify whether CL2 improves CI coverage calibration versus current cluster CR1 for low-to-moderate numbers of curves.

Secondary aims:

1. Quantify CL2 interval-width inflation relative to CR1.
2. Assess whether gains concentrate in correlated settings and attenuate under IID.
3. Verify CL2 behavior for non-Gaussian links (Poisson/Binomial) and `gaulss`.

Decision hypothesis:

1. CL2 is worth integrating if it reduces absolute coverage gap to nominal (90%) in correlated cells without unacceptable width inflation or IID over-conservatism.

## 2) Pilot Evidence and Uncertainty Caveat

Current pilots are directional motivation only, not effect-size estimation:

1. Gaussian pilot used `n_rep = 8` per cell.
2. Non-Gaussian pilot used `n_rep = 4` per cell.
3. gaulss quick test used `n_rep = 2` per cell.

Coverage MC precision for those rep counts is coarse:

1. `n_rep=8`: `MC_SE(p)` ranges from `0.106` (`p=0.9`) to `0.177` (`p=0.5`), so `±2*MC_SE` is about `±21pp` to `±35pp`.
2. `n_rep=4`: `MC_SE(p)` ranges from `0.150` to `0.250`, so `±2*MC_SE` is about `±30pp` to `±50pp`.
3. `n_rep=2`: `MC_SE(p)` ranges from `0.212` to `0.354`, so `±2*MC_SE` is about `±42pp` to `±71pp`.

Interpretation:

1. Pilot magnitudes are noisy.
2. Only directional consistency is used here: CL2 was positive in all pilot block aggregates, strongest in correlated low-`n` settings.
3. Production run is required for reliable magnitude estimation.

## 3) ADEMP

### A) Aims

As above.

### D) Data-generating mechanisms

Production benchmark has three blocks and two tiers:

1. Core tier: precision-focused factorial cells.
2. Stress tier: smallest-`n` sentinel cells (`n=20`) to probe CL2 regime where leverage corrections should matter most.

Fixed model/grids across all blocks (unless explicitly stated):

1. `nxgrid = 30`, `nygrid = 40` (coarse grid, matching existing benchmark infra).
2. `ff` basis: `k_ff = c(12, 12)`.
3. `bs.yindex`: `k = 12`, `m = c(2, 1)`.
4. Terms evaluated for decisions: `ff`, `linear` (primary), `E(Y)` (secondary).

Known grid caveat:

1. Coarse-grid intercept bias is known from Study 2.
2. Therefore intercept is diagnostic only, not a primary decision target.
3. Primary decisions are based on `ff` and `linear`.

Correlation comparability design:

1. All blocks share a common comparison slice: `corr_type in {iid, ar1(0.9)}`.
2. Additional correlation levels are block-specific exploratory extensions to keep runtime feasible.

Block G (Gaussian, ff + linear):

1. Core:
   - `n`: `20, 40, 80`
   - `corr_type`: `iid`, `ar1(0.3)`, `ar1(0.9)`, `fourier_pos(0.3)`
   - Cells: `12`
2. Stress:
   - none beyond core (`n=20` already included).

Block NG (Non-Gaussian ff + linear):

1. Core:
   - `family`: `poisson`, `binomial`
   - `n`: `25, 50`
   - `corr_type`: `iid`, `ar1(0.9)`, `fourier_pos(0.3)`
   - Cells: `12`
2. Stress (`n=20` added per review):
   - `family`: `poisson`, `binomial`
   - `n`: `20`
   - `corr_type`: `iid`, `ar1(0.9)` (no fourier for runtime control)
   - Cells: `4`

Block H (`gaulss`, heteroskedastic Gaussian):

1. Core:
   - `n`: `30, 60`
   - `corr_type`: `iid`, `ar1(0.9)`
   - `hetero_type`: `none`, `bump(3.0)`
   - Cells: `8`
2. Stress (`n=20` added per review):
   - `n`: `20`
   - `corr_type`: `iid`, `ar1(0.9)`
   - `hetero_type`: `none`, `bump(3.0)`
   - Cells: `4`

Total planned cells: `12 + 16 + 12 = 40`.

### E) Estimands

Primary:

1. Term-level pointwise CI coverage for `ff` and `linear` at nominal 90%.

Secondary:

1. `E(Y)` coverage.
2. Mean interval width.
3. `z_sd` and `z_mean` calibration diagnostics.

### M) Methods

Block G and Block NG:

1. `default`
2. `hc`
3. `cluster` (current CR1 path)
4. `cl2` (experimental path in `ci-benchmark/cl2-pilot/cl2-utils.R`)

Block H:

1. `gaulss_default`
2. `gaulss_cluster`
3. `gaulss_cl2`

### P) Performance measures

For each `(block, dgp, method, term)`:

1. Coverage.
2. Mean width.
3. RMSE and bias.
4. `z_sd`, `z_mean`.
5. Failure rate.
6. Fit time.
7. CL2 overhead time (post-fit covariance/metric extraction).

MC uncertainty reporting:

1. `MC_SE(coverage) = sqrt(p*(1-p)/n_rep)`.
2. Report `coverage +/- 2*MC_SE`.

## 4) Replication Targets (precision-driven)

Target precision:

1. Core G/NG cells: `MC_SE <= 0.0125` near `p=0.9` -> `n_rep ~= 576`; use `600`.
2. Core H cells: `MC_SE <= 0.015` near `p=0.9` -> `n_rep ~= 400`; use `400`.
3. Stress cells (NG/H `n=20`): `MC_SE <= 0.02` near `p=0.9` -> `n_rep ~= 225`; use `300`.

Planned reps:

1. Block G: `12 * 600 = 7200`.
2. Block NG core: `12 * 600 = 7200`.
3. Block NG stress: `4 * 300 = 1200`.
4. Block H core: `8 * 400 = 3200`.
5. Block H stress: `4 * 300 = 1200`.
6. Total: `20,000` dataset tasks, each analyzed by all applicable methods on the same data.

## 5) Runtime Budget (fit + CL2 overhead)

Observed pilot median fit times (model fitting only):

1. Gaussian: ~`3-6s`.
2. Non-Gaussian: ~`6-13s`.
3. gaulss: ~`40-77s`.

Measured CL2 post-fit overhead (quick timing on representative cells):

1. Gaussian AR1 n=40: cluster metric step `1.09s`, CL2 `1.10s` (ratio `1.01`, delta `+0.01s`).
2. Poisson AR1 n=50: cluster `0.869s`, CL2 `0.876s` (ratio `1.01`, delta `+0.007s`).
3. gaulss AR1 n=60: cluster `1.20s`, CL2 `1.70s` (ratio `1.42`, delta `+0.5s`).

Interpretation:

1. CL2 overhead is modest versus total fit time for G/NG.
2. gaulss has higher relative CL2 overhead in the covariance step, but still small compared with total gaulss fit time.

Approximate walltime with 6 workers:

1. Block G: `~2-4h`.
2. Block NG (core + stress): `~4-7h`.
3. Block H (core + stress): `~12-20h`.
4. Total: `~18-31h`, run budget `~36h` including retries/resume overhead.

## 6) Architecture and reproducibility

One-dataset-all-methods per `(dgp, rep)`:

1. Generate once.
2. Fit all applicable methods.
3. Save per-task RDS immediately.

Seed plan:

1. Fixed base seed per block.
2. Deterministic task seed: `seed = base_seed + 1000*dgp_id + rep_id`.
3. Record seed and error info in every output row.

Resume/dedup:

1. Incremental `dgpXXX_repYYY.rds` saves.
2. On resume, skip existing keys.
3. Before dedup, assert deterministic metric equality excluding timing fields.

## 7) Review gates and execution phases

Phase 0: infrastructure checks

1. Known-answer checks and smoke runs for G/NG/H.

R1: code review gate

1. Review benchmark scripts and CL2 utilities.
2. Resolve high-severity findings.

Phase 1: calibration pilot

1. `n_rep = 30` per core cell.
2. `n_rep = 20` per stress cell.
3. Verify ranking stability, failure rates, and do-no-harm diagnostics.

R2: pre-production review

1. Review pilot summaries and logs.
2. Confirm no pathological cells before scaling.

Phase 2: production

1. Run full core/stress reps.
2. Persist per-block summaries with MC_SE columns.

R3: results review

1. Validate conclusions against pre-specified thresholds only.

## 8) Predefined decision thresholds (no post-hoc moving)

CL2 promotion recommendation if all criteria hold:

1. Correlated-cell aggregate gain:
   - Mean absolute coverage-gap improvement vs CR1 `>= 0.02` for primary terms (`ff`, `linear`).
2. IID guardrail:
   - Mean overshoot above nominal under IID not worse than `+0.02` for primary terms.
3. Per-cell do-no-harm criterion:
   - No primary `(dgp, term)` cell has worsening `coverage_cl2 - coverage_cr1 < -max(0.02, 2*MC_SE_cell)`.
4. Width rule (quantified):
   - Correlated cells with `n >= 30`: median width ratio `CL2/CR1 <= 1.35`.
   - Stress cells with `n <= 25`: width ratio up to `1.80` allowed only if both:
     - coverage gain `>= 0.10` absolute, and
     - post-CL2 coverage `>= 0.85`.
5. Failure-rate guardrail:
   - CL2 failure rate not worse than CR1 by more than `1pp`.

## 9) Deliverables

Per block:

1. `results_combined.rds`.
2. `summary_by_cell.csv` (including MC_SE columns).
3. `summary_cluster_vs_cl2.csv`.
4. `summary_overall.csv`.
5. concise markdown report.

Cross-block:

1. Common-slice comparison table (`iid` + `ar1(0.9)`).
2. Recommendation memo for package integration (`sandwich = "cl2"` yes/no and scope).

## 10) Immediate implementation tasks

1. Build one unified production runner with block/tier configs and a shared schema.
2. Add MC_SE computation directly in summarizers.
3. Add validation checks for:
   - missing `(dgp, rep, method, term)` combinations,
   - duplicate-key deterministic consistency,
   - per-cell do-no-harm criterion,
   - failure-rate guardrail.
4. Add CL2 overhead logging (`cluster_step_time`, `cl2_step_time`) to summaries.
