# pffr Refactoring Plan

> **Implementation Note**: Use `/r-package-coding` skill for R package standards (naming, roxygen2, checkmate, cli).

---

## Goals

1. Redesign `pffrSim` as a formula-driven data generator with customizable true effects
2. Write comprehensive unit tests with ground truth comparison
3. Refactor pffr formula parsing into modular functions
4. (Future) Better defaults and inference

---

## Step-by-Step Implementation

### Phase 1: pffrSim Redesign

#### 1.1 Core Infrastructure

- [x] Create `parse_pffr_formula()` in `R/pffr-utilities.R`
  - Extract term types and variable names using `terms.formula()` with specials
  - Return structured list: `list(terms = list(type, varname, options), response = "Y")`

- [x] Create preset effect libraries (internal, not exported)
  - `ff_presets`: `cosine`, `product`, `gaussian`, `separable`, `historical`
  - `smooth_presets`: `beta`, `dnorm`, `sine`, `cosine`, `polynomial`, `step`
  - `const_presets`: `constant`, `gaussian_2d`, `linear`
  - `intercept_presets`: `constant`, `beta`, `sine`, `zero`

- [x] Create `resolve_effect(spec, term_type, xind, yind)`
  - Accept: preset string, custom function, or list with parameters
  - Return: evaluated coefficient function/matrix

#### 1.2 Effect Computation Functions

- [x] `compute_ff_effect(X, beta_st, xind)` - Integral ∫X(s)β(s,t)ds
  - Use rectangular integration weights (consistent with `ff()` defaults)

- [x] `compute_linear_effect(x, beta_t, yind)` - Varying coefficient x·β(t)

- [x] `compute_smooth_effect(x, f_xt, yind)` - Smooth f(x,t)

- [x] `compute_const_effect(x, f_x)` - Constant over t

- [x] `compute_intercept(mu_t, n)` - Functional intercept

- [ ] (Later) `compute_sff_effect()`, `compute_ffpc_effect()`, `compute_pcre_effect()`

#### 1.3 Covariate Generation

- [x] `generate_functional_covariate(n, xind, type = "bspline")`
  - Generate smooth random functions via B-spline basis

- [x] `generate_scalar_covariate(n, type)` - Normal, uniform, or factor

#### 1.4 Main pffrSim Function

- [x] Implement new signature:
  ```r
  pffrSim(formula, n = 100, yind = seq(0,1,l=60), xind = seq(0,1,l=40),
          data = NULL, effects = list(), intercept = "beta",
          SNR = 10, family = gaussian(), propmissing = 0, seed = NULL)
  ```

- [x] Add noise generation:
  - Gaussian: `SNR = var(eta) / var(epsilon)`
  - Non-Gaussian: Use `family$linkinv(eta)` then `family$simulate()` if available

- [x] Return structure with `truth` attribute:
  ```r
  attr(result, "truth") <- list(eta, etaTerms, beta, epsilon)
  ```

#### 1.5 Backward Compatibility

- [x] Deprecate `scenario` argument with `.Deprecated()` warning

- [x] Implement scenario-to-formula mappings:
  | Scenario | Formula | Effects |
  |----------|---------|---------|
  | `"int"` | `Y ~ 1` | `intercept = "beta"` |
  | `"ff"` | `Y ~ ff(X1)` | `list(X1 = "cosine")` |
  | `"lin"` | `Y ~ xlin` | `list(xlin = "dnorm")` |
  | `"smoo"` | `Y ~ s(xsmoo)` | `list(xsmoo = "sine")` |
  | `"te"` | `Y ~ c(te(xte1, xte2))` | `list(\`te(xte1,xte2)\` = "gaussian_2d")` |
  | `"const"` | `Y ~ c(xconst)` | `list(xconst = 1.5)` |
  | `"all"` | Full formula | All combined |

- [x] Verify legacy `pffrSim(scenario = "all")` produces identical output (via `pffrSim_legacy`)

#### 1.6 Documentation

- [x] Update roxygen docs in `R/pffr-utilities.R` with new API examples
- [x] Update `NEWS.md` with new interface and deprecation notice
- [x] Run `devtools::document()`

---

### Phase 2: Unit Tests

#### 2.1 Extend `tests/testthat/test-pffr.R`

Add sections with clear headers. Use `skip_on_cran()` for slow tests.

**Term Recovery Tests** (high SNR, check correlation > 0.9):
- [x] `ff()` term recovers β(s,t)
- [x] `ff()` with `limits="s<t"`
- [x] `s()` smooth term varying over t
- [x] `c()` constant-over-t terms
- [x] Linear varying coefficient terms
- [x] `sff()` term
- [x] `ffpc()` term
- [ ] Factor terms (covered in existing tests)
- [ ] `pcre()` random effects (covered in existing tests)

**Algorithm Tests** (`skip_on_cran()`, `skip_if_not_installed()`):
- [x] `algorithm="bam"` for large data
- [x] `algorithm="gamm"` with random effects
- [x] `algorithm="gamm4"`

**Data Structure Tests**:
- [x] Sparse data (`propmissing > 0`)
- [x] Different `xind`/`yind` grid sizes

**Inference Tests**:
- [x] `predict()` dimensions and values
- [x] `predict(type="terms")`
- [x] `coef()` extraction
- [x] `summary()` output

#### 2.2 Test Guidelines

- Use `n = 20-50`, `nygrid = 30-50`, `nxgrid = 20-40`
- Prefer correlation/RMSE over exact coefficient equality:
  ```r
  expect_gt(cor(as.vector(est), as.vector(true)), 0.9)
  ```
- Targeted runtimes: < 3 min locally (i.e.: prefer small data sets with high SNR)

#### 2.3 Golden Snapshots (Optional)

- [ ] Save reference mgcv formula and key data columns for backward compat checks
- [ ] Store in `tests/testthat/fixtures/`

---

### Phase 3: pffr Methods Refactoring & Test Coverage

#### 3.1 Test Coverage Gaps

**3.1.1 Factor Variables with >2 Levels**

Problem: Current tests only use binary factors; no tests for 3+ level factors with default treatment contrasts.

- [ ] Add test: factor with 3 levels as varying coefficient (`xfactor` with levels A/B/C)
- [ ] Add test: factor in interaction with smooth terms
- [ ] Verify `coef.pffr` returns correct number of coefficients per level

NOTE: possibly this requires different handling of centering constraints than other terms. 
If so, defer to future phase.

**3.1.2 Untested Methods**

Problem: Several pffr methods have 0% test coverage.

- [ ] `residuals.pffr`: test dimensions, reformatting, sparse data
- [ ] `plot.pffr`: smoke test (returns without error)
- [ ] `qq.pffr`: smoke test
- [ ] `pffr.check`: smoke test
- [ ] `print.summary.pffr`: capture output, verify key content present

**3.1.3 Edge Cases**

- [ ] `se.fit = TRUE` in predict.pffr - verify SE dimensions and values
- [ ] Very small n (n=5) still fits without error
- [ ] Model with single smooth term

#### 3.2 summary.pffr & getShrtlbls Refactoring

**3.2.1 Replace getShrtlbls with Mapping-Based Approach**

Problem: Current regex-based string parsing in `getShrtlbls()` is fragile and produces non-unique labels.

Files: `R/pffr-utilities.R` (getShrtlbls), `R/pffr.R` (where labelmap is created)

Approach: Store term→label mapping at model fit time instead of reverse-engineering from term strings.

- [ ] In `pffr()`, create `object$pffr$shortlabels` as named vector: `c("s(x,y,by=g).g1" = "g1(t)", ...)`
- [ ] Modify `summary.pffr` to use `object$pffr$shortlabels` directly
- [ ] Deprecate or simplify `getShrtlbls()` to fall back on new mapping
- [ ] Ensure unique labels even for `s(g, bs="re") + s(g, bs="mrf", ...)`

**3.2.2 Test print.summary.pffr Output**

- [ ] Test that formula is printed correctly
- [ ] Test that smooth terms table has expected row names
- [ ] Test that R-squared and deviance explained are shown
- [ ] Test sample size format: "n = 1200 (40 x 30)"

#### 3.3 predict.pffr Improvements

**3.3.1 Implement Prediction with `limits`**

Problem: `ff()` terms with `limits = "s<t"` cannot predict on new data (line 126 TODO in `R/pffr-methods.R`).

- [ ] Implement prediction for ff terms with limits
- [ ] Create test: fit model with `ff(X, limits="s<t")`, predict on new data
- [ ] Verify predicted values match expectations

**3.3.2 Test `se.fit = TRUE`**

- [ ] Test `predict(model, se.fit = TRUE)` returns list with `fit` and `se.fit`
- [ ] Verify SE dimensions match fit dimensions
- [ ] Test with different prediction types (link, response, terms)

#### 3.4 coef.pffr PCRE Fix (Optional)

Problem: `coef.pffr` fails for pcre terms with more than 2 functional principal components (line 606 FIXME in `R/pffr-methods.R`).

**NOTE**: This is a niche edge case. Investigate complexity first - if it requires deep architectural changes, defer or mark as "won't fix".

- [ ] Diagnose the specific failure mode for >2 FPCs
- [ ] If fix is straightforward: implement fix to handle arbitrary number of FPCs
- [ ] If fix is complex: document limitation and defer to future phase
- [ ] Add test: pcre term with 3 FPCs (either verifying fix or documenting expected failure)

#### 3.5 fitted.pffr gaulss Handling

Problem: `fitted.pffr` only returns means for gaulss, silently discarding scale predictions.

- [ ] Add parameter `which = c("mean", "scale", "both")` or similar
- [ ] Document the behavior for location-scale families
- [ ] Add test for gaulss fitted values extraction

#### 3.6 Formula Parsing Refactor

Extract modular functions from `pffr()` for better maintainability and testing.

- [ ] Extract `parse_pffr_formula()` from `pffr()`
- [ ] Extract term-specific transformers: `transform_ff_term()`, `transform_s_term()`, etc.
- [ ] Extract `build_mgcv_data()` and `build_mgcv_formula()`
- [ ] Use feature flag: `options(pffr.use_modular_parser = TRUE)`
- [ ] Maintain identical outputs to existing tests

#### Key Files to Modify

| File | Changes |
|------|---------|
| `R/pffr.R` | Add shortlabels creation at fit time, extract modular parsing functions |
| `R/pffr-utilities.R` | Simplify/deprecate getShrtlbls |
| `R/pffr-methods.R` | predict limits, coef PCRE fix, fitted gaulss |
| `R/pffr-formula.R` | (New) Extracted formula parsing functions |
| `tests/testthat/test-pffr.R` | All new tests |

#### Implementation Order

1. Test coverage (3.1) - add tests without changing code (safety net)
2. Formula parsing refactor (3.6) - biggest structural change, do first to avoid merge conflicts
3. getShrtlbls refactor (3.2) - fits cleanly into new modular structure
4. predict.pffr limits (3.3.1)
5. coef.pffr PCRE fix (3.4) - if complexity is reasonable
6. fitted.pffr gaulss (3.5)
7. se.fit tests (3.3.2)

---

### Phase 4: Better Defaults (Future)

- [ ] Larger default spline bases (current `bs.yindex k=5` is small)
- [ ] Consider `gaulss` or `pffrGLS` by default
- [ ] Smarter automatic algorithm selection

---

### Phase 5: Better Inference (Future)

- [ ] Simultaneous confidence bands
- [ ] Conformal prediction intervals

---

### Phase 6: Statistical Methods (Future)

Additional S3 methods for model diagnostics and comparison.

- [ ] `vcov.pffr`: Check if exists, add tests or implement
- [ ] `confint.pffr`: Confidence intervals for coefficients
- [ ] `anova.pffr`: Model comparison (nested models, likelihood ratio tests)
- [ ] `update.pffr`: Re-fit model with modified formula/data

---

## Verification Checklist

After each phase:
- [ ] `devtools::test()` passes
- [ ] `devtools::check()` no new warnings/errors
- [ ] Example models from documentation run correctly
- [ ] Legacy `pffrSim(scenario = ...)` produces same results

---

## Reference

### Key Files

| File | Lines | Purpose |
|------|-------|---------|
| `R/pffr.R` | 229-931 | Main pffr function |
| `R/pffr-utilities.R` | 152-291 | pffrSim (modify this) |
| `R/pffr-ff.R` | - | ff() term constructor |
| `R/pffr-sff.R` | - | sff() term constructor |
| `R/pffr-ffpc.R` | - | ffpc() term constructor |
| `R/pffr-pcre.R` | - | pcre() term constructor |
| `tests/testthat/test-pffr.R` | - | Tests (extend this) |

### Effects API: Term Label Matching

1. Simple terms: variable name (`"xlin"` for `xlin`)
2. ff/sff/ffpc: first argument (`"X1"` for `ff(X1, xind=s)`)
3. c() wrappers: wrapped term (`"te(xte1,xte2)"` for `c(te(xte1,xte2))`)
4. s()/te(): full term string (`"s(xsmoo)"`)
5. Repeated variables: positional suffix (`"X1.1"`, `"X1.2"`)

### Integration Weights

Use rectangular rule consistent with `ff()` defaults:
```r
dx <- diff(xind)[1]
L <- matrix(dx, nrow = n, ncol = length(xind))
effect <- (L * X) %*% beta_st
```

### Noise Generation

**Gaussian**: `sigma2 = var(eta) / SNR; Y = eta + rnorm(..., sd=sqrt(sigma2))`

**Non-Gaussian**: `mu = family$linkinv(eta); Y = family$simulate(fitted.values=mu)`

### Limitations

- pffrSim uses regular grids only (same `yind`/`xind` for all observations)
- For irregular data: use `propmissing` or construct manually


--------------------------------------

# Phase 1a: Core infrastructure
  ralph-loop "Read PLAN.md. Implement Phase 1 sections 1.1-1.3 (parse_pffr_formula, preset libraries, effect computation functions). Check boxes as done. Run Rscript -e 'devtools::load_all(); devtools::test()' to verify." --completion-promise "Sections 1.1-1.3 checkboxes marked [x] and devtools::test() passes" --max-iterations 50

  # Phase 1b: Main function + backward compat
  /raralph-loop "Read PLAN.md. Continue with Phase 1 sections 1.4-1.6 (main pffrSim, backward compat, docs). Check boxes as done." --completion-promise "All Phase 1 checkboxes marked [x] and devtools::test() passes" --max-iterations 50

  # Phase 2: Tests
  ralph-loop "Read PLAN.md. Implement Phase 2 (unit tests). Check boxes as done." --completion-promise "Phase 2 checkboxes marked [x] and devtools::test() passes" --max-iterations 8

 # Phase 3: refactor
  ralph-loop "Read PLAN.md. Implement Phase 3 (refactor & test pffr + utilities). Check boxes as done." --completion-promise "Phase 3 checkboxes marked [x] and devtools::test() passes" --max-iterations 10




