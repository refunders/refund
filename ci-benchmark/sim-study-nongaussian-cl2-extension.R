# ===========================================================================
# Study 1 CL2 Extension: Add cl2 sandwich results to existing Study 1
# ===========================================================================
#
# Existing Study 1 result files have 3 methods (default, hc, cluster).
# This script re-fits the same models with the same seeds and extracts
# only the cl2 sandwich metrics, saving to a separate directory.
#
# The key insight: we must re-fit the model (same seed → same data + fit)
# because cl2 requires the fitted model object for covariance correction.
#
# Usage:
#   Rscript ci-benchmark/sim-study-nongaussian-cl2-extension.R [mode]
#   mode: "smoke" (2 reps), "pilot" (10 reps), "full" (150 reps)
# ===========================================================================

# Setup -----------------------------------------------------------------------

library(tidyverse)

if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(refund)
}

library(mvtnorm)

source("ci-benchmark/benchmark-utils.R")
source("ci-benchmark/confint-benchmark.R")

# Source Study 1 for DGP functions, settings, formula builder, constants.
# The guarded main block (sys.nframe() == 0) won't execute when sourced.
source("ci-benchmark/sim-study-nongaussian-sandwich.R")

# CL2 Metric Extraction -------------------------------------------------------

#' Extract cl2-only metrics from one pffr fit
#'
#' Mirrors extract_study1_metrics() but only computes the cl2 method.
#' Pre-computes coef(fit, sandwich="cl2") once and passes to
#' compute_term_metrics() via the coefs argument.
#'
#' @param fit A fitted pffr model.
#' @param sim Simulation result from simulate_*_ff_linear().
#' @param alpha Significance level for CI.
#' @returns Tibble with metrics for cl2 × {ff, linear, E(Y)}.
extract_cl2_only_metrics <- function(fit, sim, alpha = 0.10) {
  # Pre-compute cl2 coefficients once (expensive sandwich correction)
  coefs_cl2 <- tryCatch(
    coef(fit, sandwich = "cl2"),
    error = function(e) {
      warning("coef(fit, sandwich='cl2') failed: ", conditionMessage(e))
      NULL
    }
  )

  results <- list()

  for (term_type in c("ff", "linear")) {
    metrics <- tryCatch(
      compute_term_metrics(
        fit = fit,
        truth = sim$truth,
        term_type = term_type,
        alpha = alpha,
        use_sandwich = "cl2",
        s_grid = sim$s_grid,
        t_grid = sim$t_grid,
        data = sim$data,
        coefs = coefs_cl2,
        err_struct = sim$err_struct
      ),
      error = function(e) {
        warning(sprintf(
          "compute_term_metrics failed for cl2 term=%s: %s",
          term_type,
          conditionMessage(e)
        ))
        NULL
      }
    )

    if (!is.null(metrics)) {
      metrics$method <- "cl2"
      results[[length(results) + 1]] <- metrics
    }
  }

  # E(Y) metrics — uses its own sandwich path (apply_sandwich_correction)
  ey_metrics <- tryCatch(
    compute_ey_metrics(
      fit = fit,
      truth = sim$truth,
      alpha = alpha,
      use_sandwich = "cl2",
      err_struct = sim$err_struct,
      t_grid = sim$t_grid
    ),
    error = function(e) {
      warning(sprintf(
        "compute_ey_metrics failed for cl2: %s",
        conditionMessage(e)
      ))
      NULL
    }
  )

  if (!is.null(ey_metrics)) {
    ey_metrics$method <- "cl2"
    results[[length(results) + 1]] <- ey_metrics
  }

  dplyr::bind_rows(results)
}

# Runner Functions -------------------------------------------------------------

#' Run one CL2 replicate (same DGP + fit as Study 1, cl2 extraction only)
#'
#' @param row One-row tibble/list with DGP settings.
#' @param rep_id Replicate number.
#' @param alpha Significance level.
#' @returns Tibble with cl2 metrics, or NULL on failure.
run_cl2_rep <- function(row, rep_id, alpha = 0.10) {
  if (inherits(row, "data.frame")) row <- as.list(row)

  seed <- STUDY1_BASE_SEED + 1000L * row$dgp_id + rep_id

  # Simulate data — same seed → identical DGP realization
  sim_fn <- if (row$family == "poisson") {
    simulate_poisson_ff_linear
  } else {
    simulate_binomial_ff_linear
  }

  sim <- sim_fn(
    n = row$n,
    nxgrid = row$nxgrid,
    nygrid = row$nygrid,
    corr_type = row$corr_type,
    corr_param = if (is.na(row$corr_param)) NULL else row$corr_param,
    seed = seed
  )

  # Build formula and fit — identical to run_study1_rep
  frml <- build_study1_formula(sim$s_grid)
  bs_yindex <- list(bs = "ps", k = STUDY1_K_YINDEX, m = c(2, 1))

  fam <- if (row$family == "poisson") poisson() else binomial()

  t0 <- Sys.time()
  fit <- pffr(
    frml,
    yind = sim$t_grid,
    data = sim$data,
    family = fam,
    bs.yindex = bs_yindex,
    sandwich = "none"
  )
  fit_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Extract cl2-only metrics
  metrics <- extract_cl2_only_metrics(fit, sim, alpha = alpha)

  if (nrow(metrics) == 0) {
    return(tibble(
      term_type = NA_character_,
      coverage = NA_real_,
      coverage_high_var = NA_real_,
      coverage_low_var = NA_real_,
      mean_width = NA_real_,
      rmse = NA_real_,
      bias = NA_real_,
      mean_abs_error = NA_real_,
      mean_se = NA_real_,
      median_se = NA_real_,
      z_mean = NA_real_,
      z_sd = NA_real_,
      z_kurtosis = NA_real_,
      z2_mean = NA_real_,
      n_grid = NA_integer_,
      method = "cl2",
      dgp_id = row$dgp_id,
      rep_id = rep_id,
      seed = seed,
      family = row$family,
      n = row$n,
      nxgrid = row$nxgrid,
      nygrid = row$nygrid,
      corr_type = row$corr_type,
      corr_param = row$corr_param,
      fit_time = fit_time,
      converged = FALSE,
      error_msg = "cl2 metrics extraction returned 0 rows"
    ))
  }

  # Annotate with settings
  metrics |>
    dplyr::mutate(
      dgp_id = row$dgp_id,
      rep_id = rep_id,
      seed = seed,
      family = row$family,
      n = row$n,
      nxgrid = row$nxgrid,
      nygrid = row$nygrid,
      corr_type = row$corr_type,
      corr_param = row$corr_param,
      fit_time = fit_time,
      converged = TRUE,
      error_msg = NA_character_
    )
}

#' Run CL2 extension for all Study 1 DGP cells
#'
#' @param n_rep Number of replications per DGP cell.
#' @param parallel Use parallel processing?
#' @param n_workers Number of workers.
#' @param output_dir Output directory for cl2 results.
#' @param alpha Significance level.
#' @returns Combined results tibble.
run_cl2_extension <- function(
  n_rep = 150L,
  parallel = TRUE,
  n_workers = max(1L, min(6L, parallel::detectCores() - 1L)),
  output_dir = "ci-benchmark/study1-nongaussian-cl2",
  alpha = 0.10
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  settings <- make_study1_settings()

  # Expand to full grid
  grid <- settings |>
    tidyr::crossing(rep_id = seq_len(n_rep))

  cat("Study 1 CL2 Extension\n")
  cat("=====================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("Replications:", n_rep, "\n")
  cat("Total fits:", nrow(grid), "\n")
  cat("Output dir:", output_dir, "\n\n")

  # Skip already-completed settings
  existing_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(existing_files) > 0) {
    existing_keys <- sub("\\.rds$", "", basename(existing_files))
    grid_keys <- sprintf("dgp%03d_rep%03d", grid$dgp_id, grid$rep_id)
    already_done <- grid_keys %in% existing_keys
    n_skip <- sum(already_done)
    if (n_skip > 0) {
      cat("Skipping", n_skip, "already-completed settings\n")
      grid <- grid[!already_done, , drop = FALSE]
    }
  }

  if (nrow(grid) == 0) {
    cat("All settings completed — loading from disk.\n")
    return(load_cl2_results(output_dir))
  }

  # Worker function
  run_one <- function(row_df) {
    row <- as.list(row_df)
    result <- tryCatch(
      run_cl2_rep(row, row$rep_id, alpha = alpha),
      error = function(e) {
        warning(sprintf(
          "run_cl2_rep failed for dgp=%d rep=%d: %s",
          row$dgp_id,
          row$rep_id,
          conditionMessage(e)
        ))
        tibble(
          term_type = NA_character_,
          coverage = NA_real_,
          coverage_high_var = NA_real_,
          coverage_low_var = NA_real_,
          mean_width = NA_real_,
          rmse = NA_real_,
          bias = NA_real_,
          mean_abs_error = NA_real_,
          mean_se = NA_real_,
          median_se = NA_real_,
          z_mean = NA_real_,
          z_sd = NA_real_,
          z_kurtosis = NA_real_,
          z2_mean = NA_real_,
          n_grid = NA_integer_,
          method = "cl2",
          dgp_id = row$dgp_id,
          rep_id = row$rep_id,
          seed = NA_integer_,
          family = row$family,
          n = row$n,
          nxgrid = row$nxgrid,
          nygrid = row$nygrid,
          corr_type = row$corr_type,
          corr_param = row$corr_param,
          fit_time = NA_real_,
          converged = FALSE,
          error_msg = conditionMessage(e)
        )
      }
    )

    # Always save, including failure rows
    if (!is.null(result) && nrow(result) > 0) {
      save_path <- file.path(
        output_dir,
        sprintf("dgp%03d_rep%03d.rds", row$dgp_id, row$rep_id)
      )
      atomic_saveRDS(result, save_path)
    }
    gc()
    result
  }

  if (
    parallel &&
      requireNamespace("furrr", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)
  ) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multicore, workers = n_workers)

    is_dev_mode <- file.exists("DESCRIPTION")
    pkg_dir <- normalizePath(".")
    bench_utils_path <- normalizePath("ci-benchmark/benchmark-utils.R")
    confint_path <- normalizePath("ci-benchmark/confint-benchmark.R")
    study1_path <- normalizePath(
      "ci-benchmark/sim-study-nongaussian-sandwich.R"
    )

    rows <- split(grid, seq_len(nrow(grid)))

    results <- furrr::future_map_dfr(
      rows,
      function(row_df) {
        if (is_dev_mode) {
          devtools::load_all(pkg_dir, quiet = TRUE)
        } else {
          library(refund)
        }
        source(bench_utils_path, local = TRUE)
        source(confint_path, local = TRUE)
        source(study1_path, local = TRUE)
        # Source this extension script for extract_cl2_only_metrics
        # (already in parent env when running sequentially, but workers
        # need explicit sourcing)
        run_one(row_df)
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  } else {
    results <- list()
    for (i in seq_len(nrow(grid))) {
      row_df <- grid[i, ]
      cat(sprintf(
        "\r[%d/%d] dgp=%d (family=%s, corr=%s, n=%d) rep=%d",
        i,
        nrow(grid),
        row_df$dgp_id,
        row_df$family,
        row_df$corr_type,
        row_df$n,
        row_df$rep_id
      ))
      results[[i]] <- run_one(row_df)
    }
    cat("\n")
    results <- dplyr::bind_rows(results)
  }

  # Load all results including previously completed
  all_results <- load_cl2_results(output_dir)

  atomic_saveRDS(all_results, file.path(output_dir, "results_combined.rds"))
  cat("Total cl2 results:", nrow(all_results), "rows\n")

  all_results
}

#' Load CL2 extension results from incremental saves
#'
#' @param output_dir Directory with per-rep RDS files.
#' @returns Combined tibble.
load_cl2_results <- function(
  output_dir = "ci-benchmark/study1-nongaussian-cl2"
) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(all_files) == 0) {
    warning("No result files found in ", output_dir)
    return(tibble())
  }
  dplyr::bind_rows(lapply(all_files, readRDS))
}

# Merge Function ---------------------------------------------------------------

#' Merge original Study 1 results with CL2 extension
#'
#' @param original_dir Directory with original Study 1 results.
#' @param cl2_dir Directory with CL2 extension results.
#' @param output_path Path for merged output file.
#' @returns Merged tibble.
merge_study1_cl2 <- function(
  original_dir = "ci-benchmark/study1-nongaussian",
  cl2_dir = "ci-benchmark/study1-nongaussian-cl2",
  output_path = "ci-benchmark/study1-nongaussian/results_combined_with_cl2.rds"
) {
  cat("Loading original results from:", original_dir, "\n")
  orig <- load_study1_results(original_dir)
  cat("  Rows:", nrow(orig), "\n")
  cat("  Methods:", paste(unique(orig$method), collapse = ", "), "\n")

  cat("Loading CL2 results from:", cl2_dir, "\n")
  cl2 <- load_cl2_results(cl2_dir)
  cat("  Rows:", nrow(cl2), "\n")
  cat("  Methods:", paste(unique(cl2$method), collapse = ", "), "\n")

  # Verify seed consistency: for each (dgp_id, rep_id) pair, seeds must match
  orig_seeds <- orig |>
    dplyr::filter(!is.na(seed)) |>
    dplyr::distinct(dgp_id, rep_id, seed) |>
    dplyr::rename(seed_orig = seed)

  cl2_seeds <- cl2 |>
    dplyr::filter(!is.na(seed)) |>
    dplyr::distinct(dgp_id, rep_id, seed) |>
    dplyr::rename(seed_cl2 = seed)

  seed_check <- dplyr::inner_join(
    orig_seeds,
    cl2_seeds,
    by = c("dgp_id", "rep_id")
  )

  mismatches <- seed_check |>
    dplyr::filter(seed_orig != seed_cl2)

  if (nrow(mismatches) > 0) {
    warning(
      "Seed mismatches found for ",
      nrow(mismatches),
      " (dgp_id, rep_id) pairs!"
    )
    print(head(mismatches, 10))
  } else {
    cat("Seed verification: OK (", nrow(seed_check), " pairs checked)\n")
  }

  # Merge: bind rows
  merged <- dplyr::bind_rows(orig, cl2)
  cat("Merged:", nrow(merged), "rows\n")
  cat("Methods:", paste(sort(unique(merged$method)), collapse = ", "), "\n")

  # Save
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  atomic_saveRDS(merged, output_path)
  cat("Saved to:", output_path, "\n")

  merged
}

# Main Entry Point -------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args) > 0) args[1] else "smoke"

  # Special mode: merge existing results
  if (mode == "merge") {
    merged <- merge_study1_cl2()
    cat("\nMerge complete.\n")
    quit(save = "no")
  }

  n_rep <- switch(
    mode,
    smoke = 2L,
    pilot = 10L,
    full = 150L,
    as.integer(mode)
  )
  if (is.na(n_rep)) n_rep <- 2L

  cat("Running Study 1 CL2 Extension in", mode, "mode (", n_rep, "reps)\n\n")

  results <- run_cl2_extension(
    n_rep = n_rep,
    parallel = (n_rep > 5),
    n_workers = max(1L, min(6L, parallel::detectCores() - 1L))
  )

  if (nrow(results) > 0) {
    summary_df <- results |>
      dplyr::filter(!is.na(coverage)) |>
      dplyr::group_by(family, corr_type, corr_param, n, method, term_type) |>
      dplyr::summarise(
        mean_coverage = mean(coverage, na.rm = TRUE),
        se_coverage = sd(coverage, na.rm = TRUE) / sqrt(dplyr::n()),
        mean_width = mean(mean_width, na.rm = TRUE),
        mean_z_sd = mean(z_sd, na.rm = TRUE),
        n_reps = dplyr::n(),
        .groups = "drop"
      )

    cat("\n========== CL2 EXTENSION COVERAGE SUMMARY ==========\n")
    print(
      summary_df |>
        dplyr::arrange(family, corr_type, n, term_type),
      n = 200
    )
  }

  cat("\nCL2 extension complete.\n")
}
