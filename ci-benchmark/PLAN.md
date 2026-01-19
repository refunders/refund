# CI Benchmark Restructuring Plan

## Current Problems Identified

### 1. Evaluation is broken (pointwise instead of term-aggregated)
The current `evaluate_fit()` returns a data frame with one row per grid point per term. The `summarize_results()` does aggregate, but:
- Returns pointwise `sq_err`, `width`, `covered` in the scored data frame
- The aggregation is overly complex and happens too late
- Coverage should be: "for this term, what fraction of the pointwise CIs covered the truth?"
- Width should be: "for this term, what is the mean CI width?"

### 2. Inefficient fitting structure
Current flow:
```
for each (dgp, rep, method) combination:
    generate data  <-- wasteful! same dgp+rep generates same data
    fit model
    evaluate
```

Should be:
```
for each (dgp, rep):
    generate data ONCE
    fit pffr (base model)
    for each method:
        fit model (reusing pffr fit if needed)
        evaluate
```

### 3. Over-complex code (~2000 lines)
- `simulate_dataset()`: 200+ lines, does too much
- `run_benchmark()`: 400+ lines with nested logic
- `make_truth_functions()`: 90 lines of bespoke random function generation
- Many helper functions that obscure the data flow

### 4. Grid mismatch risk
`coef.pffr()` evaluates on its own grid (n1, n2, n3 parameters). Truth is evaluated on the same grid by calling the truth functions, but this is fragile. Better to use `predict(type="terms")` on the original data grid.

---

## Proposed Structure

### Core Design Principles
1. **Data-frame centric**: Settings, results, and metrics all flow through tibbles
2. **One dataset, all methods**: Generate data once per (dgp_id, rep_id), fit all methods
3. **Aggregate early**: Return one row per (dgp_id, rep_id, method, term_type) with metrics
4. **Simple truth**: Use `pffr_simulate()` directly with known truth, avoid custom truth generation
5. **Consistent grids**: Evaluate on same grids as DGP

### File Structure
```
ci-benchmark.R          # Main benchmark script (~300 lines)
```

Single file with sections:
1. Setup & packages
2. DGP specification (data generation)
3. Model fitting functions
4. Metric computation
5. Benchmark runner
6. Execution

---

## Detailed Implementation Plan

### Section 1: Setup (~30 lines)

```r
library(tidyverse)
library(furrr)       # for parallel map
library(progressr)   # for progress bars

devtools::load_all(".")  # load refund dev version
```

### Section 2: DGP Specification (~50 lines)

**Settings tibble**: Each row is a unique DGP configuration
```r
make_dgp_settings <- function() {

  crossing(
    n = 80L,
    nxgrid = 35L,
    nygrid = 45L,
    snr = c(5, 15, 50),
    error_dist = c("gaussian", "t6"),
    corr_type = c("iid", "ar1"),
    corr_param = NA_real_,
    hetero_type = c("none", "linear"),
    terms = list(c("ff", "concurrent", "linear", "smooth"))
  ) |>
    mutate(dgp_id = row_number())
}
```

**Methods tibble**: Which methods to fit, with feasibility constraints
```r
make_method_settings <- function() {
  tribble(
    ~method,         ~needs_pffr_fit, ~family_ok,
    "pffr",          FALSE,           c("gaussian", "gaulss", "scat"),
    "pffr_sandwich", TRUE,            c("gaussian", "gaulss", "scat"),
    "pffr_gls",      TRUE,            c("gaussian"),
    "pffr_ar",       TRUE,            c("gaussian"),
  )
}
```

**Expand to full benchmark grid**:
```r
make_benchmark_grid <- function(dgp_settings, method_settings, n_rep, seed) {
  # Cross dgp with reps
  dgp_reps <- dgp_settings |>
    crossing(rep_id = 1:n_rep) |>
    mutate(seed = seed + 1000L * dgp_id + rep_id)

  # Add methods column as list (filtered by feasibility)
  dgp_reps |>
    mutate(
      methods = map2(error_dist, hetero_type, ~{
        method_settings |>
          filter(
            # gaulss only for heteroskedastic
            !(method == "pffr" & .y == "none" & "gaulss" %in% family_ok),
            # scat only for t-distributed errors
            !(method == "pffr" & .x != "t6" & "scat" %in% family_ok)
          ) |>
          pull(method)
      })
    )
}
```

### Section 3: Data Generation (~80 lines)

**Single function to generate data with known ground truth**:
```r
generate_benchmark_data <- function(
    n, nxgrid, nygrid, snr,
    error_dist, corr_type, corr_param, hetero_type,
    terms, seed
) {
  set.seed(seed)

  # Build formula from terms
  formula_parts <- c(
    if ("ff" %in% terms) "ff(X1)" else NULL,
    if ("concurrent" %in% terms) "Xc" else NULL,
    if ("linear" %in% terms) "zlin" else NULL,
    if ("smooth" %in% terms) "s(zsmoo)" else NULL
  )
  formula <- as.formula(paste("Y ~", paste(formula_parts, collapse = " + ")))

  # Generate base data using pffr_simulate
  dat <- pffr_simulate(
    formula,
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    SNR = snr
  )

  s_grid <- attr(dat, "xindex")
  t_grid <- attr(dat, "yindex")
  truth <- attr(dat, "truth")

  # Add structured errors if needed
  if (corr_type != "iid" || hetero_type != "none" || error_dist != "gaussian") {
    dat <- add_structured_errors(dat, error_dist, corr_type, corr_param, hetero_type, snr)
  }

  # Add concurrent covariate if needed
  if ("concurrent" %in% terms) {
    dat$Xc <- I(generate_functional_covariate(n, t_grid))
    # Add concurrent effect to truth...
  }

  list(
    data = dat,
    s_grid = s_grid,
    t_grid = t_grid,
    truth = truth,
    dgp_params = list(corr_type = corr_type, hetero_type = hetero_type)
  )
}
```

### Section 4: Model Fitting (~80 lines)

**Fit all methods for one dataset, with timing**:
```r
fit_all_methods <- function(sim, methods, family = "gaussian") {
  dat <- sim$data
  t_grid <- sim$t_grid
  s_grid <- sim$s_grid

  # Build formula
  formula <- build_pffr_formula(sim)

  fits <- list()
  timings <- list()

  # Always fit pffr first (base model) - this is the pilot fit for GLS/AR
  t0 <- Sys.time()
  fits$pffr <- pffr(formula, yind = t_grid, data = dat, family = family)
  pffr_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timings$pffr <- pffr_time

  # Fit other methods, reusing pffr fit where needed
  if ("pffr_sandwich" %in% methods) {
    fits$pffr_sandwich <- fits$pffr  # same fit, different coef extraction
    timings$pffr_sandwich <- pffr_time  # same time as pffr
  }

  if ("pffr_gls" %in% methods) {
    t0 <- Sys.time()
    hatSigma <- estimate_hatSigma(fits$pffr, dat)
    fits$pffr_gls <- pffr_gls(formula, yind = t_grid, data = dat, hatSigma = hatSigma)
    gls_only_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    # Total time = pilot pffr fit + GLS fit
    timings$pffr_gls <- pffr_time + gls_only_time
    timings$pffr_gls_only <- gls_only_time  # just the GLS step
  }

  if ("pffr_ar" %in% methods) {
    t0 <- Sys.time()
    rho <- estimate_rho(fits$pffr, dat)
    fits$pffr_ar <- pffr(formula, yind = t_grid, data = dat,
                          algorithm = "bam", rho = rho)
    ar_only_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    # Total time = pilot pffr fit + AR fit
    timings$pffr_ar <- pffr_time + ar_only_time
    timings$pffr_ar_only <- ar_only_time  # just the AR step
  }

  list(fits = fits, timings = timings, pffr_pilot_time = pffr_time)
}
```

### Section 5: Metric Computation (~100 lines)

**Compute metrics for one fit, one term (aggregated, not pointwise)**:
```r
compute_term_metrics <- function(fit, truth, term_type, alpha = 0.10, use_sandwich = FALSE) {
  # Extract coefficients on appropriate grid
  coefs <- coef(fit, sandwich = use_sandwich, n1 = 50, n2 = 25, n3 = 15)

  # Find the term
  term_info <- coefs$smterms[[find_term_by_type(coefs, term_type)]]
  if (is.null(term_info)) return(NULL)

  # Get estimates and SEs
  est <- term_info$coef$value
  se <- term_info$coef$se


  # Get truth on same grid
  truth_vals <- evaluate_truth_on_grid(truth, term_type, term_info)

  # Compute CI
  z <- qnorm(1 - alpha / 2)
  lower <- est - z * se
  upper <- est + z * se
  width <- 2 * z * se

  # Aggregate metrics
  covered <- (truth_vals >= lower) & (truth_vals <= upper)

  tibble(
    term_type = term_type,
    coverage = mean(covered),
    mean_width = mean(width),
    rmse = sqrt(mean((est - truth_vals)^2)),
    n_grid = length(est)
  )
}

compute_all_metrics <- function(fit, sim, method, use_sandwich = FALSE, alpha = 0.10) {
  truth <- sim$truth
  term_types <- attr(sim$data, "truth")$term_set %||% c("ff", "concurrent", "linear", "smooth")

  map_dfr(term_types, ~{
    compute_term_metrics(fit, truth, .x, alpha, use_sandwich) |>
      mutate(method = method)
  })
}
```

### Section 6: Benchmark Runner (~80 lines)

**Process one (dgp_id, rep_id) combination**:
```r
run_one_dgp_rep <- function(row, alpha = 0.10) {
  # Generate data
  sim <- generate_benchmark_data(
    n = row$n, nxgrid = row$nxgrid, nygrid = row$nygrid,
    snr = row$snr, error_dist = row$error_dist,
    corr_type = row$corr_type, corr_param = row$corr_param,
    hetero_type = row$hetero_type, terms = row$terms[[1]],
    seed = row$seed
  )

  # Fit all applicable methods (returns fits + timings)
  methods <- row$methods[[1]]
  fit_result <- fit_all_methods(sim, methods)
  fits <- fit_result$fits
  timings <- fit_result$timings

  # Compute metrics for each method
  results <- map_dfr(names(fits), function(method) {
    use_sandwich <- (method == "pffr_sandwich")
    fit <- if (method == "pffr_sandwich") fits$pffr else fits[[method]]

    compute_all_metrics(fit, sim, method, use_sandwich, alpha) |>
      mutate(
        dgp_id = row$dgp_id,
        rep_id = row$rep_id,
        seed = row$seed,
        # Timing info
        fit_time_total = timings[[method]],
        fit_time_step_only = timings[[paste0(method, "_only")]] %||% timings[[method]],
        pffr_pilot_time = fit_result$pffr_pilot_time
      )
  })

  results
}
```

**Main runner with parallelization**:
```r
run_benchmark <- function(
    dgp_settings = make_dgp_settings(),
    n_rep = 50,
    seed = 2024,
    n_workers = parallelly::availableCores() - 1,
    output_dir = "ci-benchmark"
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Create benchmark grid
  grid <- make_benchmark_grid(dgp_settings, make_method_settings(), n_rep, seed)

  # Setup parallel backend
  plan(multisession, workers = n_workers)

  # Run with progress
  handlers(global = TRUE)
  with_progress({
    p <- progressor(nrow(grid))

    results <- future_pmap_dfr(grid, function(...) {
      row <- tibble(...)
      result <- run_one_dgp_rep(row)

      # Incremental save
      save_path <- file.path(output_dir,
                             sprintf("dgp%03d_rep%03d.rds", row$dgp_id, row$rep_id))
      saveRDS(result, save_path)

      p()
      result
    }, .options = furrr_options(seed = TRUE))
  })

  # Save combined results
  saveRDS(results, file.path(output_dir, "results_combined.rds"))

  results
}
```

### Section 7: Summarization (~40 lines)

```r
summarize_benchmark <- function(results) {
  results |>
    group_by(method, term_type, dgp_id) |>
    summarize(
      mean_coverage = mean(coverage),
      se_coverage = sd(coverage) / sqrt(n()),
      mean_width = mean(mean_width),
      mean_rmse = mean(rmse),
      # Timing summaries
      mean_fit_time_total = mean(fit_time_total),
      mean_fit_time_step = mean(fit_time_step_only),
      mean_pffr_pilot_time = mean(pffr_pilot_time),
      n_rep = n(),
      .groups = "drop"
    )
}

# Separate timing summary (aggregated across terms)
summarize_timings <- function(results) {
  results |>
    distinct(dgp_id, rep_id, method, .keep_all = TRUE) |>
    group_by(method, dgp_id) |>
    summarize(
      mean_fit_time_total = mean(fit_time_total),
      sd_fit_time_total = sd(fit_time_total),
      mean_fit_time_step = mean(fit_time_step_only),
      mean_pffr_pilot_time = mean(pffr_pilot_time),
      n_rep = n(),
      .groups = "drop"
    )
}
```

---

## Implementation Steps (Incremental)

### Step 1: Data Generation
- [ ] Implement `generate_benchmark_data()`
- [ ] Test that it produces correct truth structure
- [ ] Verify covariates are properly centered

### Step 2: Model Fitting
- [ ] Implement `fit_all_methods()`
- [ ] Test each method runs without error
- [ ] Verify pffr_gls/pffr_ar properly reuse pffr fit

### Step 3: Metric Computation
- [ ] Implement `compute_term_metrics()`
- [ ] Implement `evaluate_truth_on_grid()` - key function for grid alignment
- [ ] Test metrics on known data (e.g., coverage ~0.9 at SNR=25)

### Step 4: Benchmark Loop
- [ ] Implement `run_one_dgp_rep()`
- [ ] Implement `run_benchmark()` with parallelization
- [ ] Test on tiny grid (2 dgps, 2 reps)

### Step 5: Full Run
- [ ] Run on small grid (10 dgps, 10 reps)
- [ ] Verify results match expectations
- [ ] Run full benchmark

---

## Key Simplifications

1. **Drop custom truth generation**: Use `pffr_simulate()` with standard effects
2. **Drop visualization code**: Move to separate analysis script
3. **Drop resume logic**: Use incremental saves, manually filter completed
4. **Drop complex error structures**: Keep iid + ar1 only for initial version
5. **Single family per run**: Don't mix gaussian/gaulss/scat in same grid

---

## Expected Code Size

| Section | Lines |
|---------|-------|
| Setup | 30 |
| DGP spec | 50 |
| Data gen | 80 |
| Fitting | 80 |
| Metrics | 100 |
| Runner | 80 |
| Summary | 50 |
| **Total** | **~470** |

Compared to current ~2000 lines, this is a 75% reduction.

---

## Output Data Structure

Each row in the results tibble represents one (dgp, rep, method, term) combination:

```
# Results tibble columns:
dgp_id          : int     # DGP configuration ID
rep_id          : int     # Replication number
seed            : int     # RNG seed for reproducibility
method          : chr     # "pffr", "pffr_sandwich", "pffr_gls", "pffr_ar"
term_type       : chr     # "ff", "concurrent", "linear", "smooth"
coverage        : dbl     # Fraction of grid points where CI covered truth
mean_width      : dbl     # Mean CI width across grid
rmse            : dbl     # RMSE of estimates vs truth
n_grid          : int     # Number of grid points evaluated
fit_time_total  : dbl     # Total fit time in seconds (includes pilot for gls/ar)
fit_time_step_only: dbl   # Just this method's step (excl. pilot)
pffr_pilot_time : dbl     # Time for pilot pffr fit (same for all methods in a dgp:rep)
snr             : dbl     # Signal-to-noise ratio
corr_type       : chr     # Correlation type
hetero_type     : chr     # Heteroskedasticity type
error_dist      : chr     # Error distribution
```

This structure allows:
1. Easy filtering by method, term type, or DGP
2. Grouping for summary statistics
3. Joining with DGP settings for analysis
4. Separate timing analysis
