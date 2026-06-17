# ===========================================================================
# Study 1 Competitors Extension: t_{G-1} and simultaneous-band CIs
# ===========================================================================
#
# Re-fits each (DGP, rep) from Study 1 (same seeds, same model) and extracts
# only the competitor-CI metrics:
#
#   E1 (t_{G-1}): for use_sandwich in c("cluster", "cl2"), call
#       compute_term_metrics(..., crit_df = G-1); tag method "<sw>_t".
#
#   E2 (simultaneous): for sandwich in c("none","cluster","cl2"), compute
#       coef(fit, ci="simultaneous", level=0.90, n_sim=2000,
#            sandwich=m, sim_seed=seed+777L)
#       once per m, then extract_simultaneous_term_metrics() per term;
#       tag method "<m>_sim" (none→"default_sim").
#
# Study 1 DGP:
#   - Poisson(log) and Binomial(logit) families
#   - ff(X1) + zlin model; no "smooth" or "concurrent" terms
#   - Corr: iid, ar1(0.9), fourier_pos(0.3); n in {200, 400}; 12 DGP cells
#   - Seed scheme: STUDY1_BASE_SEED + 1000*dgp_id + rep_id (= 3001 base)
#
# Usage:
#   Rscript ci-benchmark/sim-study1-competitors-extension.R [mode]
#   mode: "smoke" (1 rep), "pilot" (10 reps), "full" (150 reps)
# ===========================================================================

# Setup -----------------------------------------------------------------------

library(tidyverse)
library(mvtnorm)

library(refund)

source("ci-benchmark/benchmark-utils.R")
source("ci-benchmark/confint-benchmark.R")

# The main guard (sys.nframe() == 0) prevents execution when sourced.
source("ci-benchmark/sim-study-nongaussian-sandwich.R")

# Constants -------------------------------------------------------------------

COMPETITORS1_OUTPUT_DIR <- "ci-benchmark/study1-competitors"
COMPETITORS1_LEVEL <- 0.90
COMPETITORS1_N_SIM <- 2000L
COMPETITORS1_SIM_SEED_OFFSET <- 777L

# Metric extraction -----------------------------------------------------------

#' Extract E1 (t_{G-1}) and E2 (simultaneous) competitor metrics for one fit
#'
#' @param fit Fitted pffr model (sandwich="none").
#' @param sim Simulation result from simulate_*_ff_linear().
#' @param seed Integer seed used for this rep.
#' @param G Number of curves (= n; df = G-1 for t-CI).
#' @param alpha Significance level (default 0.10 → 90% CIs).
#' @returns Tibble with competitor metrics, or empty tibble on complete failure.
extract_competitors_study1 <- function(fit, sim, seed, G, alpha = 0.10) {
  terms_all <- c(sim$terms, "intercept")
  results <- list()

  # --- E1: t_{G-1} critical values for cluster and cl2 ----------------------

  df_t <- G - 1L

  for (sw in c("cluster", "cl2")) {
    method_name <- paste0(sw, "_t")

    # Pre-compute coefs once per sandwich type
    coefs_t <- tryCatch(
      coef(fit, sandwich = sw, seWithMean = FALSE, n1 = 50, n2 = 25, n3 = 15),
      error = function(e) {
        warning(sprintf(
          "coef(sandwich='%s') failed in E1: %s",
          sw,
          conditionMessage(e)
        ))
        NULL
      }
    )

    for (term_type in terms_all) {
      tm <- tryCatch(
        compute_term_metrics(
          fit = fit,
          truth = sim$truth,
          term_type = term_type,
          alpha = alpha,
          use_sandwich = sw,
          s_grid = sim$s_grid,
          t_grid = sim$t_grid,
          data = sim$data,
          coefs = coefs_t,
          err_struct = sim$err_struct,
          crit_df = df_t
        ),
        error = function(e) {
          warning(sprintf(
            "compute_term_metrics E1 method=%s term=%s: %s",
            method_name,
            term_type,
            conditionMessage(e)
          ))
          NULL
        }
      )
      if (!is.null(tm)) {
        tm$method <- method_name
        results[[length(results) + 1L]] <- tm
      }
    }
  }

  # --- E2: simultaneous bands ------------------------------------------------

  sim_sandwich_map <- c(
    none = "default_sim",
    cluster = "cluster_sim",
    cl2 = "cl2_sim"
  )

  for (m in c("none", "cluster", "cl2")) {
    method_name <- unname(sim_sandwich_map[m])
    sim_seed <- seed + COMPETITORS1_SIM_SEED_OFFSET

    coefs_sim <- tryCatch(
      coef(
        fit,
        sandwich = m,
        ci = "simultaneous",
        level = COMPETITORS1_LEVEL,
        n_sim = COMPETITORS1_N_SIM,
        sim_seed = sim_seed,
        seWithMean = FALSE,
        n1 = 50,
        n2 = 25,
        n3 = 15
      ),
      error = function(e) {
        warning(sprintf(
          "coef(ci='simultaneous', sandwich='%s') failed: %s",
          m,
          conditionMessage(e)
        ))
        NULL
      }
    )

    if (is.null(coefs_sim)) next

    for (term_type in terms_all) {
      sm <- tryCatch(
        extract_simultaneous_term_metrics(
          fit = fit,
          truth = sim$truth,
          term_type = term_type,
          use_sandwich = m,
          coefs_sim = coefs_sim,
          s_grid = sim$s_grid,
          t_grid = sim$t_grid,
          data = sim$data
        ),
        error = function(e) {
          warning(sprintf(
            "extract_simultaneous_term_metrics method=%s term=%s: %s",
            method_name,
            term_type,
            conditionMessage(e)
          ))
          NULL
        }
      )
      if (!is.null(sm)) {
        sm$method <- method_name
        results[[length(results) + 1L]] <- sm
      }
    }
  }

  dplyr::bind_rows(results)
}

# Runner Functions -------------------------------------------------------------

#' Null-row skeleton for Study 1 competitor metrics
#'
#' @param row DGP settings list.
#' @param rep_id Integer replicate id.
#' @param seed Integer seed.
#' @param error_msg Character error description.
#' @returns One-row tibble.
make_null_competitor1_row <- function(
  row,
  rep_id,
  seed,
  error_msg = NA_character_
) {
  tibble(
    term_type = NA_character_,
    method = NA_character_,
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
    coverage_joint = NA_real_,
    coverage_pointwise = NA_real_,
    dgp_id = row$dgp_id,
    rep_id = rep_id,
    seed = seed,
    family = row$family,
    n = row$n,
    nxgrid = row$nxgrid,
    nygrid = row$nygrid,
    corr_type = row$corr_type,
    corr_param = row$corr_param,
    converged = FALSE,
    error_msg = error_msg
  )
}

#' Run one Study 1 competitor replicate
#'
#' Simulates data with the SAME seed as the production study, refits the model,
#' and extracts only the competitor-CI metrics (E1: t_{G-1}, E2: simultaneous).
#'
#' @param row One-row tibble/list with DGP settings.
#' @param rep_id Replicate number.
#' @param output_dir Output directory.
#' @param alpha Significance level.
#' @returns Tibble with competitor metrics, or null row on failure.
run_competitors1_rep <- function(row, rep_id, output_dir, alpha = 0.10) {
  if (inherits(row, "data.frame")) row <- as.list(row)

  seed <- STUDY1_BASE_SEED + 1000L * row$dgp_id + rep_id
  G <- row$n # G = number of curves = sample size

  save_path <- file.path(
    output_dir,
    sprintf("dgp%03d_rep%03d.rds", row$dgp_id, rep_id)
  )

  # Skip if already done
  if (file.exists(save_path)) {
    obj <- tryCatch(readRDS(save_path), error = function(e) NULL)
    if (!is.null(obj) && nrow(obj) > 0) return(obj)
  }

  # Simulate data — same seed → same DGP realization as production
  sim_fn <- if (row$family == "poisson") {
    simulate_poisson_ff_linear
  } else {
    simulate_binomial_ff_linear
  }

  sim <- tryCatch(
    sim_fn(
      n = row$n,
      nxgrid = row$nxgrid,
      nygrid = row$nygrid,
      corr_type = row$corr_type,
      corr_param = if (is.na(row$corr_param)) NULL else row$corr_param,
      seed = seed
    ),
    error = function(e) {
      warning(sprintf(
        "Simulation failed dgp=%d rep=%d: %s",
        row$dgp_id,
        rep_id,
        conditionMessage(e)
      ))
      NULL
    }
  )

  if (is.null(sim)) {
    result <- make_null_competitor1_row(row, rep_id, seed, "simulation failed")
    atomic_saveRDS(result, save_path)
    return(result)
  }

  # Build formula and fit — identical to production run_study1_rep()
  frml <- build_study1_formula(sim$s_grid)
  bs_yindex <- list(bs = "ps", k = STUDY1_K_YINDEX, m = c(2, 1))
  fam <- if (row$family == "poisson") poisson() else binomial()

  fit <- tryCatch(
    pffr(
      frml,
      yind = sim$t_grid,
      data = sim$data,
      family = fam,
      bs.yindex = bs_yindex,
      sandwich = "none"
    ),
    error = function(e) {
      warning(sprintf(
        "pffr failed dgp=%d rep=%d: %s",
        row$dgp_id,
        rep_id,
        conditionMessage(e)
      ))
      NULL
    }
  )

  if (is.null(fit)) {
    result <- make_null_competitor1_row(row, rep_id, seed, "pffr() failed")
    atomic_saveRDS(result, save_path)
    return(result)
  }

  metrics <- tryCatch(
    extract_competitors_study1(fit, sim, seed = seed, G = G, alpha = alpha),
    error = function(e) {
      warning(sprintf(
        "extract_competitors_study1 failed dgp=%d rep=%d: %s",
        row$dgp_id,
        rep_id,
        conditionMessage(e)
      ))
      tibble()
    }
  )

  if (nrow(metrics) == 0) {
    result <- make_null_competitor1_row(
      row,
      rep_id,
      seed,
      "metrics extraction returned 0 rows"
    )
  } else {
    result <- metrics |>
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
        converged = TRUE,
        error_msg = NA_character_
      )
  }

  atomic_saveRDS(result, save_path)
  gc()
  result
}

# Loader ----------------------------------------------------------------------

#' Load Study 1 competitor extension results from incremental saves
#'
#' @param output_dir Directory with per-rep RDS files.
#' @returns Combined tibble.
load_competitors1_results <- function(output_dir = COMPETITORS1_OUTPUT_DIR) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(all_files) == 0) {
    warning("No competitor result files found in ", output_dir)
    return(tibble())
  }
  dplyr::bind_rows(lapply(all_files, readRDS))
}

# Main runner -----------------------------------------------------------------

#' Run Study 1 competitors extension
#'
#' @param n_rep Number of replications per DGP cell.
#' @param parallel Use parallel processing?
#' @param n_workers Number of parallel workers.
#' @param output_dir Output directory for results.
#' @param alpha Significance level.
#' @returns Combined results tibble (invisibly).
run_competitors1 <- function(
  n_rep = 150L,
  parallel = TRUE,
  n_workers = min(3L, max(1L, parallel::detectCores() - 1L)),
  output_dir = COMPETITORS1_OUTPUT_DIR,
  alpha = 0.10
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  settings <- make_study1_settings()
  grid <- settings |> tidyr::crossing(rep_id = seq_len(n_rep))

  cat("Study 1 Competitors Extension (t_{G-1} and simultaneous bands)\n")
  cat("================================================================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("Replications:", n_rep, "\n")
  cat("Total fits:", nrow(grid), "\n")
  cat("Output dir:", output_dir, "\n\n")

  # Skip already-completed files
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
      cat("Skipping", n_skip, "already-completed reps\n")
      grid <- grid[!already_done, , drop = FALSE]
    }
  }

  if (nrow(grid) == 0) {
    cat("All reps done — loading from disk.\n")
    return(invisible(load_competitors1_results(output_dir)))
  }

  run_one <- function(row_df) {
    row <- as.list(row_df)
    tryCatch(
      run_competitors1_rep(
        row,
        row$rep_id,
        output_dir = output_dir,
        alpha = alpha
      ),
      error = function(e) {
        warning(sprintf(
          "run_competitors1_rep failed dgp=%d rep=%d: %s",
          row$dgp_id,
          row$rep_id,
          conditionMessage(e)
        ))
        make_null_competitor1_row(
          row,
          row$rep_id,
          STUDY1_BASE_SEED + 1000L * row$dgp_id + row$rep_id,
          conditionMessage(e)
        )
      }
    )
  }

  if (
    parallel &&
      n_workers > 1L &&
      requireNamespace("furrr", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)
  ) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multicore, workers = n_workers)

    pkg_dir <- normalizePath(".")
    bench_utils_path <- normalizePath("ci-benchmark/benchmark-utils.R")
    confint_path <- normalizePath("ci-benchmark/confint-benchmark.R")
    study1_path <- normalizePath(
      "ci-benchmark/sim-study-nongaussian-sandwich.R"
    )
    ext_path <- normalizePath(
      "ci-benchmark/sim-study1-competitors-extension.R"
    )

    rows <- split(grid, seq_len(nrow(grid)))

    results <- furrr::future_map_dfr(
      rows,
      function(row_df) {
        library(refund)
        library(mvtnorm)
        source(bench_utils_path, local = TRUE)
        source(confint_path, local = TRUE)
        source(study1_path, local = TRUE)
        source(ext_path, local = TRUE)
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

  all_results <- load_competitors1_results(output_dir)
  atomic_saveRDS(
    all_results,
    file.path(output_dir, "results_combined.rds")
  )
  cat("Total competitor rows:", nrow(all_results), "\n")

  invisible(all_results)
}

# Main Entry Point -------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args) > 0) args[1] else "smoke"

  n_rep <- switch(
    mode,
    smoke = 1L,
    pilot = 10L,
    full = 150L,
    as.integer(mode)
  )
  if (is.na(n_rep)) n_rep <- 1L

  n_workers <- if (mode == "smoke") 1L else
    min(3L, max(1L, parallel::detectCores() - 1L))
  do_parallel <- (n_rep > 5L) && (n_workers > 1L)

  cat(
    "Running Study 1 Competitors Extension in",
    mode,
    "mode (",
    n_rep,
    "reps,",
    n_workers,
    "workers)\n\n"
  )

  results <- run_competitors1(
    n_rep = n_rep,
    parallel = do_parallel,
    n_workers = n_workers
  )

  if (nrow(results) > 0 && any(!is.na(results$method))) {
    cat("\n========== COMPETITOR EXTENSION SUMMARY ==========\n")

    # E1: t_{G-1} pointwise coverage
    e1_summary <- results |>
      dplyr::filter(grepl("_t$", method) & !is.na(coverage)) |>
      dplyr::group_by(family, corr_type, n, method, term_type) |>
      dplyr::summarise(
        mean_coverage = mean(coverage, na.rm = TRUE),
        mean_width = mean(mean_width, na.rm = TRUE),
        n_reps = dplyr::n(),
        .groups = "drop"
      )

    if (nrow(e1_summary) > 0) {
      cat("\nE1 (t_{G-1}) coverage summary:\n")
      print(
        e1_summary |>
          dplyr::arrange(family, corr_type, n, method, term_type),
        n = 80
      )
    }

    # E2: simultaneous joint and pointwise coverage
    e2_summary <- results |>
      dplyr::filter(grepl("_sim$", method) & !is.na(coverage_joint)) |>
      dplyr::group_by(family, corr_type, n, method, term_type) |>
      dplyr::summarise(
        mean_joint = mean(coverage_joint, na.rm = TRUE),
        mean_pointwise = mean(coverage_pointwise, na.rm = TRUE),
        mean_width = mean(mean_width, na.rm = TRUE),
        n_reps = dplyr::n(),
        .groups = "drop"
      )

    if (nrow(e2_summary) > 0) {
      cat("\nE2 (simultaneous bands) coverage summary:\n")
      print(
        e2_summary |>
          dplyr::arrange(family, corr_type, n, method, term_type),
        n = 80
      )
    }
  }

  cat("\nStudy 1 Competitors Extension complete.\n")
}
