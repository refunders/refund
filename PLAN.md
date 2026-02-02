# pffr refactor: plan / status

This file tracks what is implemented on this branch, what remains open, and the
next intended refactor steps for the `pffr`/simulation work.

## Status (as of 2026-01-26)

### Implemented

- Package metadata aligned with code:
  - `Depends: R (>= 4.1.0)` (branch uses `\(x)` and `|>` syntax).
  - `DESCRIPTION` `Collate:` updated to include `R/pffr-core.R` and `R/pffr-simulate.R`.
- `pffr()` refactor groundwork:
  - Formula parsing extracted to `R/pffr-formula.R`.
  - Core fitting helpers extracted to `R/pffr-core.R`.
  - `pffr()` uses core helpers for validation, call building, label mapping, and metadata.
- AR(1) support (`rho`) tightened:
  - Defaults to `algorithm = "bam"` when `rho` is used (unless user explicitly sets a non-bam algorithm, which errors).
  - Enforces `method = "fREML"` when `rho` is used.
  - For non-Gaussian families with `rho`, auto-sets `discrete = TRUE` unless user explicitly supplies `discrete = FALSE` (then errors).
- `pffr_simulate()` redesign:
  - Exported API is `pffr_simulate()` (formula interface recommended).
  - Deprecated wrapper `pffrSim()` retained for compatibility.
  - Entry points moved to `R/pffr-simulate.R`; supporting helpers remain in `R/pffr-utilities.R`.
  - Effect key matching for term labels/calls is whitespace-insensitive.
- Dependency hygiene:
  - Removed accidental usage of `cli::cli_abort()` (back to base `stop()` messages).
- Robust inference hardening:
  - Bootstrap CI extraction (`boot::boot.ci()`) wrapped defensively to return `NA` bounds on failures.


### Validation performed

- Targeted tests re-run locally: `tests/testthat/test-pffr-ar.R`, `tests/testthat/test-pffrSim.R`, `tests/testthat/test-pffr.R` (with `NOT_CRAN=true` so `skip_on_cran()` blocks execute).
- `R CMD build .` succeeds.

## Roadmap

### 1) `pffr_simulate()` (formula simulator)

- [x] Formula parser for simulation (`parse_pffr_formula()` in `R/pffr-utilities.R`).
- [x] Preset libraries + random truth generator, including `wiggliness` and `k_truth`.
- [x] Formula simulation for: intercept, linear terms, `c()` constants, `s()` smooths, `te/ti/t2`, `ff()` and `sff()` (treated like `ff()`).
- [x] Effect key normalization:
  - For `effects=...` keys, spacing is ignored for term labels/calls.

Open:
- []  Remove backward compatibility layer
  - `pffrSim()` deprecated wrapper can point to scenario-based simulation via `pffrSim_legacy()`.
- [ ] Use `scenario_to_formula()` (`R/pffr-utilities.R`) for scenario-string to formula conversion.
- [ ] Add formula-interface simulation support for `ffpc()` terms.
- [ ] Add formula-interface simulation support for `pcre()` terms.
- [ ] Add explicit documentation/examples for effect-key conventions beyond whitespace normalization.

### 2) `pffr()` modularization and duplication reduction

- [x] Core helper module introduced (`R/pffr-core.R`) and used by `pffr()`.
- [x] Formula parsing helpers introduced (`R/pffr-formula.R`) and used by `pffr()` and `pffr_gls()`.

Open:
- [x] Replace remaining inline steps in `pffr()` with existing core helpers:
  - dimension detection via `pffr_get_dimensions()`
  - y-index handling via `pffr_setup_yind_dense()` / `pffr_setup_yind_sparse()`
  - response setup via `pffr_setup_response()`
- [x] Extract a single internal “prepare” step (`pffr_prepare()`) that returns:
  - transformed mgcv formula,
  - constructed `pffr_data`,
  - configured mgcv call (from `pffr_build_call()`),
  - postprocessing inputs (term/label mapping inputs, `where_specials`, etc.).
  This should let `pffr()` become “prepare → fit → postprocess”.

Deferred (affects `pffr_gls()`; not implemented in this change set):
- [ ] Unify `pffr_gls()` with `pffr()` by reusing `pffr_prepare()` and related helpers.
- [ ] Decide status of sparse `pffr_gls()` (implement vs explicitly unsupported).

### 3) AR(1) (`rho`) behavior and docs

- [x] Implementation now matches documented mgcv constraints (bam + fREML; discrete for non-Gaussian).

Open:
- [ ] Review user-facing docs in `R/pffr.R`/man pages for any remaining mismatches (esp. around automatic defaults vs user-supplied overrides).

### 4) API consistency (naming/export surface)

Open:
- [ ] Decide on the long-term naming convention for exported helpers:
  - This branch exports new snake_case helpers (with dotted legacy wrappers deprecated).
  - Repo historically uses dotted exports; decide whether to standardize on one convention and update docs/deprecations accordingly.

### 5) Tests and checks

- [x] Added/updated tests for AR(1), simulation, and refactor touchpoints.

Open:
- [ ] Run full `devtools::test()` and `devtools::check()` locally before merge/release and resolve any new NOTES/WARNINGS.
- [ ] Consider adding regression tests around `scenario_to_formula()` if it is kept/used.

## Notes / conventions

### Effect key matching (simulation)

Current matching order in `pffrSim_formula()` for selecting a term’s effect:
1. Exact term label (as produced by `terms.formula()`).
2. Whitespace-normalized term label.
3. If the term has a parsable call: exact `safeDeparse(call)`.
4. If the term has a parsable call: whitespace-normalized `safeDeparse(call)`.
5. Variable name fallback (first argument for `ff()`/`sff()`/`te()` etc).
6. Type-specific default.

### Integration weights (simulation)

The simulator’s `ff()` effect computation uses simple numerical integration
consistent with the `ff()` defaults (rectangular/riemann-style weights).
