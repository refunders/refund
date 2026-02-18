#!/usr/bin/env Rscript
# ============================================================================
# CL2 Pilot Study (Non-Gaussian): Low-Sample Calibration
# ============================================================================
#
# Focused pilot to test CL2-style cluster correction under non-Gaussian,
# non-identity-link pffr models (Poisson/log and Binomial/logit).
#
# Design:
#   - Families: poisson, binomial
#   - Model: ff + linear
#   - n in {25, 50}
#   - Correlation: iid and AR1(0.9)
#   - Methods: default, hc, cluster (CR1), cl2 (experimental)
#
# Usage:
#   Rscript ci-benchmark/cl2-pilot/sim-study-cl2-pilot-nongaussian.R [mode]
#   mode: "smoke" (1 rep), "pilot" (4 reps), "full" (10 reps)
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(refund)
}

source("ci-benchmark/benchmark-utils.R")
source("ci-benchmark/confint-benchmark.R")
source("ci-benchmark/sim-study-nongaussian-sandwich.R")
source("ci-benchmark/cl2-pilot/cl2-utils.R")

CL2_NG_BASE_SEED <- 6201L
CL2_NG_ALPHA <- 0.10
CL2_NG_TARGET_COVERAGE <- 1 - CL2_NG_ALPHA

make_cl2_ng_settings <- function() {
  tidyr::crossing(
    family = c("poisson", "binomial"),
    corr_type = c("iid", "ar1"),
    n = c(25L, 50L)
  ) |>
    dplyr::mutate(
      corr_param = dplyr::case_when(
        corr_type == "ar1" ~ 0.9,
        TRUE ~ NA_real_
      ),
      nxgrid = STUDY1_NXGRID,
      nygrid = STUDY1_NYGRID,
      dgp_id = dplyr::row_number()
    )
}

atomic_saveRDS <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) {
    warning("atomic_saveRDS: file.rename failed for ", path)
    file.copy(tmp, path, overwrite = TRUE)
    file.remove(tmp)
  }
}

make_ng_failure_row <- function(row, rep_id, seed, fit_time, method, msg) {
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
    method = method,
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
    error_msg = msg,
    cl2_capped_clusters = NA_integer_
  )
}

run_one_cl2_ng_rep <- function(row, rep_id, alpha = CL2_NG_ALPHA) {
  if (inherits(row, "data.frame")) row <- as.list(row)

  seed <- CL2_NG_BASE_SEED + 1000L * row$dgp_id + rep_id

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

  base_metrics <- tryCatch(
    extract_study1_metrics(fit, sim, alpha = alpha),
    error = function(e) {
      warning("extract_study1_metrics failed: ", conditionMessage(e))
      tibble()
    }
  )

  if (nrow(base_metrics) > 0) {
    base_metrics <- base_metrics |>
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
        error_msg = NA_character_,
        cl2_capped_clusters = NA_integer_
      )
  } else {
    base_metrics <- make_ng_failure_row(
      row = row,
      rep_id = rep_id,
      seed = seed,
      fit_time = fit_time,
      method = NA_character_,
      msg = "default/hc/cluster metrics extraction returned 0 rows"
    )
  }

  cl2_out <- tryCatch(
    cl2_compute_metrics(
      fit = fit,
      sim = sim,
      alpha = alpha,
      method_name = "cl2"
    ),
    error = function(e) e
  )

  if (inherits(cl2_out, "error")) {
    cl2_metrics <- make_ng_failure_row(
      row = row,
      rep_id = rep_id,
      seed = seed,
      fit_time = fit_time,
      method = "cl2",
      msg = paste0("cl2 failed: ", conditionMessage(cl2_out))
    )
  } else {
    cl2_metrics <- cl2_out$metrics |>
      dplyr::filter(term_type %in% c("ff", "linear", "E(Y)")) |>
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
        error_msg = NA_character_,
        cl2_capped_clusters = cl2_out$n_capped_clusters
      )

    if (nrow(cl2_metrics) == 0) {
      cl2_metrics <- make_ng_failure_row(
        row = row,
        rep_id = rep_id,
        seed = seed,
        fit_time = fit_time,
        method = "cl2",
        msg = "cl2 returned 0 rows after term filtering"
      )
    }
  }

  dplyr::bind_rows(base_metrics, cl2_metrics)
}

load_cl2_ng_results <- function(output_dir) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(all_files)) return(tibble())
  dplyr::bind_rows(lapply(all_files, readRDS))
}

run_cl2_ng_pilot <- function(
  n_rep = 4L,
  seed = CL2_NG_BASE_SEED,
  output_dir = "ci-benchmark/cl2-pilot/nongaussian/results",
  alpha = CL2_NG_ALPHA
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  settings <- make_cl2_ng_settings()
  grid <- settings |>
    tidyr::crossing(rep_id = seq_len(n_rep)) |>
    dplyr::mutate(seed = seed + 1000L * dgp_id + rep_id)

  cat("CL2 non-Gaussian pilot\n")
  cat("======================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("Replications:", n_rep, "\n")
  cat("Total fits:", nrow(grid), "\n")
  cat("Output dir:", output_dir, "\n\n")

  existing_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = FALSE
  )
  existing_keys <- sub("\\.rds$", "", existing_files)
  grid <- grid |>
    dplyr::mutate(key = sprintf("dgp%03d_rep%03d", dgp_id, rep_id))

  pending <- grid |>
    dplyr::filter(!(key %in% existing_keys))

  cat("Already completed:", nrow(grid) - nrow(pending), "\n")
  cat("Pending:", nrow(pending), "\n\n")

  if (nrow(pending) == 0) {
    return(load_cl2_ng_results(output_dir))
  }

  for (i in seq_len(nrow(pending))) {
    row <- pending[i, , drop = FALSE]
    row_list <- as.list(row)

    cat(sprintf(
      "[%d/%d] dgp=%d family=%s corr=%s n=%d rep=%d\n",
      i,
      nrow(pending),
      row_list$dgp_id,
      row_list$family,
      row_list$corr_type,
      row_list$n,
      row_list$rep_id
    ))

    res <- tryCatch(
      run_one_cl2_ng_rep(row_list, rep_id = row_list$rep_id, alpha = alpha),
      error = function(e) {
        make_ng_failure_row(
          row = row_list,
          rep_id = row_list$rep_id,
          seed = row_list$seed,
          fit_time = NA_real_,
          method = NA_character_,
          msg = conditionMessage(e)
        )
      }
    )

    save_path <- file.path(output_dir, paste0(row_list$key, ".rds"))
    atomic_saveRDS(res, save_path)
  }

  all_results <- load_cl2_ng_results(output_dir)
  atomic_saveRDS(all_results, file.path(output_dir, "results_combined.rds"))
  all_results
}

summarize_cl2_ng <- function(results, alpha = CL2_NG_ALPHA) {
  ok <- results |>
    dplyr::filter(
      converged %in% TRUE,
      !is.na(method),
      !is.na(term_type),
      !is.na(coverage)
    )

  by_cell <- ok |>
    dplyr::group_by(family, n, corr_type, corr_param, term_type, method) |>
    dplyr::summarise(
      coverage = mean(coverage, na.rm = TRUE),
      mean_width = mean(mean_width, na.rm = TRUE),
      z_sd = mean(z_sd, na.rm = TRUE),
      reps = dplyr::n_distinct(rep_id),
      n_rows = dplyr::n(),
      mean_capped_clusters = mean(cl2_capped_clusters, na.rm = TRUE),
      .groups = "drop"
    )

  compare_cluster_cl2 <- by_cell |>
    dplyr::filter(method %in% c("cluster", "cl2")) |>
    dplyr::select(
      family,
      n,
      corr_type,
      corr_param,
      term_type,
      method,
      coverage,
      mean_width,
      z_sd
    ) |>
    tidyr::pivot_wider(
      names_from = method,
      values_from = c(coverage, mean_width, z_sd)
    )

  required_cols <- c(
    "coverage_cluster",
    "coverage_cl2",
    "mean_width_cluster",
    "mean_width_cl2"
  )
  for (nm in required_cols) {
    if (!nm %in% names(compare_cluster_cl2)) {
      compare_cluster_cl2[[nm]] <- NA_real_
    }
  }

  compare_cluster_cl2 <- compare_cluster_cl2 |>
    dplyr::mutate(
      coverage_gap_cluster = abs(.data$coverage_cluster - (1 - alpha)),
      coverage_gap_cl2 = abs(.data$coverage_cl2 - (1 - alpha)),
      coverage_gap_improvement = coverage_gap_cluster - coverage_gap_cl2,
      width_ratio_cl2_vs_cluster = .data$mean_width_cl2 /
        .data$mean_width_cluster
    )

  overall <- ok |>
    dplyr::group_by(family, term_type, method) |>
    dplyr::summarise(
      coverage = mean(coverage, na.rm = TRUE),
      mean_width = mean(mean_width, na.rm = TRUE),
      z_sd = mean(z_sd, na.rm = TRUE),
      reps = dplyr::n_distinct(rep_id),
      dgp = dplyr::n_distinct(dgp_id),
      .groups = "drop"
    )

  list(
    by_cell = by_cell,
    compare_cluster_cl2 = compare_cluster_cl2,
    overall = overall,
    target_coverage = CL2_NG_TARGET_COVERAGE
  )
}

write_cl2_ng_report <- function(summary_obj, path, mode, n_rep) {
  cmp <- summary_obj$compare_cluster_cl2 |>
    dplyr::arrange(dplyr::desc(coverage_gap_improvement))

  overall <- summary_obj$overall |>
    dplyr::arrange(family, term_type, method)

  lines <- c(
    "# CL2 Non-Gaussian Pilot Report",
    "",
    sprintf("- Generated: %s", format(Sys.time(), usetz = TRUE)),
    sprintf("- Mode: %s", mode),
    sprintf("- Replications per DGP: %d", n_rep),
    sprintf("- Target coverage: %.2f", summary_obj$target_coverage),
    "",
    "## Overall",
    "```",
    capture.output(print(overall, n = Inf)),
    "```",
    "",
    "## Cluster vs CL2",
    "```",
    capture.output(print(cmp, n = Inf)),
    "```"
  )

  writeLines(lines, path)
}

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0) args[[1]] else "pilot"

n_rep <- switch(
  mode,
  smoke = 1L,
  pilot = 4L,
  full = 10L,
  as.integer(mode)
)
if (is.na(n_rep) || n_rep <= 0) {
  stop("Mode must be one of: smoke, pilot, full, or a positive integer.")
}

base_dir <- "ci-benchmark/cl2-pilot/nongaussian"
results_dir <- file.path(base_dir, "results")

cat(
  "Running non-Gaussian CL2 pilot in",
  mode,
  "mode (",
  n_rep,
  "reps per DGP)\n\n"
)

results <- run_cl2_ng_pilot(
  n_rep = n_rep,
  seed = CL2_NG_BASE_SEED,
  output_dir = results_dir,
  alpha = CL2_NG_ALPHA
)

cat("\nTotal result rows:", nrow(results), "\n")

if (nrow(results) > 0) {
  summary_obj <- summarize_cl2_ng(results, alpha = CL2_NG_ALPHA)
  atomic_saveRDS(summary_obj, file.path(base_dir, "summary.rds"))
  write.csv(
    summary_obj$by_cell,
    file.path(base_dir, "summary_by_cell.csv"),
    row.names = FALSE
  )
  write.csv(
    summary_obj$compare_cluster_cl2,
    file.path(base_dir, "summary_cluster_vs_cl2.csv"),
    row.names = FALSE
  )
  write.csv(
    summary_obj$overall,
    file.path(base_dir, "summary_overall.csv"),
    row.names = FALSE
  )
  write_cl2_ng_report(
    summary_obj,
    file.path(base_dir, "pilot_report.md"),
    mode,
    n_rep
  )
  cat("\nSaved summaries to:", base_dir, "\n")
}
