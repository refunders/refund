# ===========================================================================
# Study 2 Competitors Extension: t_{G-1} and simultaneous-band CIs
# ===========================================================================
#
# Re-fits each (DGP, rep) from Study 2 (same seeds, same model) and extracts
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
# Usage:
#   Rscript ci-benchmark/sim-study2-competitors-extension.R [mode]
#   mode: "smoke" (1 rep), "pilot" (10 reps), "full" (50 reps)
# ===========================================================================

# Setup -----------------------------------------------------------------------

library(tidyverse)

library(refund)

source("ci-benchmark/benchmark-utils.R")
source("ci-benchmark/confint-benchmark.R")

# The main guard (sys.nframe() == 0) prevents execution when sourced.
source("ci-benchmark/sim-study-grid-refinement.R")

# Constants -------------------------------------------------------------------

COMPETITORS_OUTPUT_DIR <- "ci-benchmark/study2-competitors"
COMPETITORS_LEVEL <- 0.90
COMPETITORS_N_SIM <- 2000L
COMPETITORS_SIM_SEED_OFFSET <- 777L

# Metric extraction -----------------------------------------------------------

#' Extract E1 (t_{G-1}) and E2 (simultaneous) competitor metrics for one fit
#'
#' @param fit Fitted pffr model (sandwich="none").
#' @param sim Simulation result at one grid level.
#' @param seed Integer seed used for this rep (for deterministic sim_seed).
#' @param G Number of curves (= n; used for df = G-1).
#' @param alpha Significance level for CI (default 0.10 → 90% CIs).
#' @returns Tibble with competitor metrics, or empty tibble on complete failure.
extract_competitors_grid <- function(fit, sim, seed, G, alpha = 0.10) {
  terms_all <- c(sim$terms, "intercept")
  results <- list()

  # --- E1: t_{G-1} critical values for cluster and cl2 ----------------------

  df_t <- G - 1L

  for (sw in c("cluster", "cl2")) {
    method_name <- paste0(sw, "_t")

    # Pre-compute coefs once per sandwich type (expensive sandwich correction)
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
    sim_seed <- seed + COMPETITORS_SIM_SEED_OFFSET

    coefs_sim <- tryCatch(
      coef(
        fit,
        sandwich = m,
        ci = "simultaneous",
        level = COMPETITORS_LEVEL,
        n_sim = COMPETITORS_N_SIM,
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

#' Null-row skeleton for competitor metrics
#'
#' Returns one NA row covering both E1 (coverage/mean_width) and E2
#' (coverage_joint/coverage_pointwise) metrics so a failed rep still occupies
#' a row in the combined tibble.
#'
#' @param row DGP settings list.
#' @param rep_id Integer replicate id.
#' @param seed Integer seed.
#' @param grid_label Character grid label.
#' @param grid_info List with nxgrid/nygrid.
#' @param error_msg Character error description.
#' @returns One-row tibble.
make_null_competitor_row <- function(
  row,
  rep_id,
  seed,
  grid_label,
  grid_info,
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
    corr_type = row$corr_type,
    corr_param = row$corr_param,
    n = row$n,
    snr = row$snr,
    nxgrid = grid_info$nxgrid,
    nygrid = grid_info$nygrid,
    grid_label = grid_label,
    converged = FALSE,
    error_msg = error_msg
  )
}

#' File key for one competitor result file
#'
#' @param dgp_id Integer.
#' @param n Integer sample size.
#' @param grid_label Character.
#' @param rep_id Integer.
#' @returns Character filename (no extension).
competitors2_file_key <- function(dgp_id, n, grid_label, rep_id) {
  nygrid <- study2_parse_nygrid(grid_label)
  sprintf("dgp%03d_n%03d_y%03d_rep%03d", dgp_id, n, nygrid, rep_id)
}

#' Run one (DGP, rep) pair across all grid levels — competitors extension
#'
#' Same seed and model as Study 2 main, only extracting competitor methods.
#'
#' @param row DGP settings list.
#' @param rep_id Integer replicate id.
#' @param grid_labels Character vector of grid labels to process.
#' @param output_dir Character path to output directory.
#' @param alpha Significance level.
#' @returns Tibble combining competitor metrics for all grid levels.
run_one_pair_competitors2 <- function(
  row,
  rep_id,
  grid_labels,
  output_dir,
  alpha = 0.10
) {
  if (inherits(row, "data.frame")) row <- as.list(row)

  seed <- STUDY2_BASE_SEED + 1000L * row$dgp_id + rep_id
  G <- row$n # G = number of curves

  # Generate paired data (same seed → same realization as production)
  sims <- tryCatch(
    generate_paired_data(row, seed),
    error = function(e) {
      warning(sprintf(
        "generate_paired_data failed dgp=%d rep=%d: %s",
        row$dgp_id,
        rep_id,
        conditionMessage(e)
      ))
      NULL
    }
  )
  if (is.null(sims)) {
    return(dplyr::bind_rows(lapply(grid_labels, function(gl) {
      make_null_competitor_row(
        row,
        rep_id,
        seed,
        gl,
        STUDY2_GRIDS[[gl]],
        "generate_paired_data failed"
      )
    })))
  }

  pair_results <- list()

  for (grid_label in grid_labels) {
    file_key <- competitors2_file_key(row$dgp_id, row$n, grid_label, rep_id)
    save_path <- file.path(output_dir, paste0(file_key, ".rds"))
    grid_info <- STUDY2_GRIDS[[grid_label]]

    # Skip if already done
    if (file.exists(save_path)) {
      obj <- tryCatch(readRDS(save_path), error = function(e) NULL)
      if (!is.null(obj) && nrow(obj) > 0) {
        pair_results[[grid_label]] <- obj
        next
      }
    }

    sim <- sims[[grid_label]]

    # Fit model — identical call to Study 2 main
    frml <- build_pffr_formula(sim, sim$s_grid, k_smooth = 12, k_ff = c(12, 12))
    bs_yindex <- list(bs = "ps", k = 12, m = c(2, 1))

    fit <- tryCatch(
      pffr(
        frml,
        yind = sim$t_grid,
        data = sim$data,
        bs.yindex = bs_yindex,
        sandwich = "none"
      ),
      error = function(e) {
        warning(sprintf(
          "pffr failed dgp=%d rep=%d grid=%s: %s",
          row$dgp_id,
          rep_id,
          grid_label,
          conditionMessage(e)
        ))
        NULL
      }
    )

    if (is.null(fit)) {
      result <- make_null_competitor_row(
        row,
        rep_id,
        seed,
        grid_label,
        grid_info,
        "pffr() failed"
      )
      atomic_saveRDS(result, save_path)
      pair_results[[grid_label]] <- result
      next
    }

    metrics <- tryCatch(
      extract_competitors_grid(fit, sim, seed = seed, G = G, alpha = alpha),
      error = function(e) {
        warning(sprintf(
          "extract_competitors_grid failed dgp=%d rep=%d grid=%s: %s",
          row$dgp_id,
          rep_id,
          grid_label,
          conditionMessage(e)
        ))
        tibble()
      }
    )

    if (nrow(metrics) == 0) {
      result <- make_null_competitor_row(
        row,
        rep_id,
        seed,
        grid_label,
        grid_info,
        "metrics extraction returned 0 rows"
      )
    } else {
      result <- metrics |>
        dplyr::mutate(
          dgp_id = row$dgp_id,
          rep_id = rep_id,
          seed = seed,
          corr_type = row$corr_type,
          corr_param = row$corr_param,
          n = row$n,
          snr = row$snr,
          nxgrid = grid_info$nxgrid,
          nygrid = grid_info$nygrid,
          grid_label = grid_label,
          converged = TRUE,
          error_msg = NA_character_
        )
    }

    atomic_saveRDS(result, save_path)
    pair_results[[grid_label]] <- result
    gc()
  }

  dplyr::bind_rows(pair_results)
}

# Loader ----------------------------------------------------------------------

#' Load Study 2 competitor extension results from incremental saves
#'
#' @param output_dir Directory with per-rep RDS files.
#' @returns Combined tibble.
load_competitors2_results <- function(output_dir = COMPETITORS_OUTPUT_DIR) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_n\\d+_y\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(all_files) == 0) {
    warning("No competitor result files found in ", output_dir)
    return(tibble())
  }
  dplyr::bind_rows(lapply(all_files, readRDS))
}

# Main runner -----------------------------------------------------------------

#' Run Study 2 competitors extension
#'
#' @param n_rep Number of replications per DGP cell.
#' @param grid_labels Character vector of grid labels.
#' @param parallel Use parallel processing?
#' @param n_workers Number of parallel workers.
#' @param output_dir Output directory for results.
#' @param alpha Significance level.
#' @returns Combined results tibble (invisibly).
run_competitors2 <- function(
  n_rep = 50L,
  grid_labels = names(STUDY2_GRIDS),
  parallel = TRUE,
  n_workers = min(3L, max(1L, parallel::detectCores() - 1L)),
  output_dir = COMPETITORS_OUTPUT_DIR,
  alpha = 0.10
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  settings <- make_study2_settings()
  task_grid <- settings |> tidyr::crossing(rep_id = seq_len(n_rep))

  cat("Study 2 Competitors Extension (t_{G-1} and simultaneous bands)\n")
  cat("================================================================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("Replications:", n_rep, "\n")
  cat("Grid labels:", paste(grid_labels, collapse = ", "), "\n")
  cat(
    "Total (DGP, rep) pairs:",
    nrow(task_grid),
    "× nygrid:",
    length(grid_labels),
    "=",
    nrow(task_grid) * length(grid_labels),
    "fits\n"
  )
  cat("Output dir:", output_dir, "\n\n")

  run_one <- function(row_df) {
    row <- as.list(row_df)
    tryCatch(
      run_one_pair_competitors2(
        row,
        rep_id = row$rep_id,
        grid_labels = grid_labels,
        output_dir = output_dir,
        alpha = alpha
      ),
      error = function(e) {
        warning(sprintf(
          "run_one_pair_competitors2 failed dgp=%d rep=%d: %s",
          row$dgp_id,
          row$rep_id,
          conditionMessage(e)
        ))
        NULL
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
    study2_path <- normalizePath("ci-benchmark/sim-study-grid-refinement.R")
    ext_path <- normalizePath(
      "ci-benchmark/sim-study2-competitors-extension.R"
    )

    rows <- split(task_grid, seq_len(nrow(task_grid)))

    results <- furrr::future_map_dfr(
      rows,
      function(row_df) {
        library(refund)
        source(bench_utils_path, local = TRUE)
        source(confint_path, local = TRUE)
        source(study2_path, local = TRUE)
        source(ext_path, local = TRUE)
        run_one(row_df)
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  } else {
    results <- list()
    for (i in seq_len(nrow(task_grid))) {
      row_df <- task_grid[i, ]
      cat(sprintf(
        "\r[%d/%d] dgp=%d (corr=%s, n=%d, snr=%g), rep=%d",
        i,
        nrow(task_grid),
        row_df$dgp_id,
        row_df$corr_type,
        row_df$n,
        row_df$snr,
        row_df$rep_id
      ))
      results[[i]] <- run_one(row_df)
    }
    cat("\n")
    results <- dplyr::bind_rows(results)
  }

  all_results <- load_competitors2_results(output_dir)
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
    full = 50L,
    as.integer(mode)
  )
  if (is.na(n_rep)) n_rep <- 1L

  n_workers <- if (mode == "smoke") 1L else
    min(3L, max(1L, parallel::detectCores() - 1L))
  do_parallel <- (n_rep > 5L) && (n_workers > 1L)

  cat(
    "Running Study 2 Competitors Extension in",
    mode,
    "mode (",
    n_rep,
    "reps,",
    n_workers,
    "workers)\n\n"
  )

  results <- run_competitors2(
    n_rep = n_rep,
    parallel = do_parallel,
    n_workers = n_workers
  )

  if (nrow(results) > 0 && any(!is.na(results$method))) {
    cat("\n========== COMPETITOR EXTENSION SUMMARY ==========\n")

    # E1: pointwise coverage from t-CI methods
    e1_summary <- results |>
      dplyr::filter(
        grepl("_t$", method) & !is.na(coverage)
      ) |>
      dplyr::group_by(corr_type, n, snr, grid_label, method, term_type) |>
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
          dplyr::arrange(corr_type, n, snr, grid_label, method, term_type),
        n = 80
      )
    }

    # E2: joint and pointwise coverage from simultaneous bands
    e2_summary <- results |>
      dplyr::filter(
        grepl("_sim$", method) & !is.na(coverage_joint)
      ) |>
      dplyr::group_by(corr_type, n, snr, grid_label, method, term_type) |>
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
          dplyr::arrange(corr_type, n, snr, grid_label, method, term_type),
        n = 80
      )
    }
  }

  cat("\nStudy 2 Competitors Extension complete.\n")
}
