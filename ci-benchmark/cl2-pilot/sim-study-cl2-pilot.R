#!/usr/bin/env Rscript
# ============================================================================
# CL2 Pilot Study: Low-Sample Functional Data Sandwich Calibration
# ============================================================================
#
# Very small pilot to test whether a CL2-style cluster correction is likely
# worth implementing for pffr in low-cluster settings.
#
# Design:
#   - Gaussian family
#   - ff + linear model
#   - n in {20, 40}
#   - correlation structures: iid and AR1(0.9)
#   - methods: default, HC, cluster (current CR1), experimental CL2-style
#
# Usage:
#   Rscript ci-benchmark/cl2-pilot/sim-study-cl2-pilot.R [mode]
#   mode: "smoke" (2 reps), "pilot" (8 reps), "full" (20 reps)
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
source("ci-benchmark/cl2-pilot/cl2-utils.R")

CL2_BASE_SEED <- 5201L
CL2_ALPHA <- 0.10
CL2_TARGET_COVERAGE <- 1 - CL2_ALPHA

# ---------------------------------------------------------------------------
# Study design helpers
# ---------------------------------------------------------------------------

make_cl2_pilot_settings <- function() {
  corr_settings <- tibble(
    corr_type = c("iid", "ar1"),
    corr_param = c(NA_real_, 0.9)
  )

  tidyr::crossing(
    n = c(20L, 40L),
    nxgrid = 30L,
    nygrid = 40L,
    snr = 25,
    wiggliness = 5,
    error_dist = "gaussian",
    hetero_type = "none",
    hetero_param = NA_real_,
    corr_settings
  ) |>
    dplyr::mutate(
      dgp_id = dplyr::row_number(),
      terms = list(c("ff", "linear"))
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

# CL2 implementation is sourced from cl2-utils.R.

run_one_cl2_setting <- function(row, alpha = CL2_ALPHA) {
  if (inherits(row, "data.frame")) row <- as.list(row)

  out <- run_one_setting(
    row = row,
    alpha = alpha,
    methods = c("pffr", "pffr_sandwich", "pffr_hc")
  )

  base_results <- out$results
  if (!"cl2_capped_clusters" %in% names(base_results)) {
    base_results$cl2_capped_clusters <- NA_integer_
  }
  if (!is.null(out$error)) return(base_results)

  cl2_out <- tryCatch(
    cl2_compute_metrics(
      fit = out$fits$pffr,
      sim = out$sim,
      alpha = alpha,
      method_name = "pffr_cl2"
    ),
    error = function(e) e
  )

  if (inherits(cl2_out, "error")) {
    fail_row <- make_failure_metrics_row(
      row,
      paste("pffr_cl2:", cl2_out$message)
    )
    fail_row$method <- "pffr_cl2"
    fail_row$cl2_capped_clusters <- NA_integer_
    return(dplyr::bind_rows(base_results, fail_row))
  }

  cl2_metrics <- cl2_out$metrics
  if (is.null(cl2_metrics) || nrow(cl2_metrics) == 0) {
    fail_row <- make_failure_metrics_row(row, "pffr_cl2: no metrics returned")
    fail_row$method <- "pffr_cl2"
    fail_row$cl2_capped_clusters <- NA_integer_
    return(dplyr::bind_rows(base_results, fail_row))
  }

  fit_time <- out$timings$pffr %||% NA_real_
  cl2_metrics <- cl2_metrics |>
    dplyr::mutate(
      dgp_id = row$dgp_id %||% NA_integer_,
      rep_id = row$rep_id %||% NA_integer_,
      seed = row$seed %||% NA_integer_,
      n = row$n %||% NA_integer_,
      nxgrid = row$nxgrid %||% NA_integer_,
      nygrid = row$nygrid %||% NA_integer_,
      snr = row$snr %||% NA_real_,
      wiggliness = row$wiggliness %||% NA_real_,
      corr_type = row$corr_type %||% NA_character_,
      corr_param = row$corr_param %||% NA_real_,
      hetero_type = row$hetero_type %||% NA_character_,
      hetero_param = row$hetero_param %||% NA_real_,
      error_dist = row$error_dist %||% NA_character_,
      fit_time_total = fit_time,
      fit_time_step_only = 0,
      pffr_pilot_time = fit_time,
      converged = TRUE,
      error_msg = NA_character_,
      cl2_capped_clusters = cl2_out$n_capped_clusters
    )

  dplyr::bind_rows(base_results, cl2_metrics)
}

# ---------------------------------------------------------------------------
# Runner + summaries
# ---------------------------------------------------------------------------

load_cl2_results <- function(output_dir) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(all_files)) return(tibble())
  dplyr::bind_rows(lapply(all_files, readRDS))
}

run_cl2_pilot <- function(
  n_rep = 8L,
  seed = CL2_BASE_SEED,
  output_dir = "ci-benchmark/cl2-pilot/results",
  alpha = CL2_ALPHA
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  settings <- make_cl2_pilot_settings()
  grid <- settings |>
    tidyr::crossing(rep_id = seq_len(n_rep)) |>
    dplyr::mutate(seed = seed + 1000L * dgp_id + rep_id)

  cat("CL2 pilot grid:", nrow(grid), "cells (dgp x rep)\n")
  cat("Output dir:", output_dir, "\n\n")

  done_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = FALSE
  )
  done_keys <- sub("\\.rds$", "", done_files)
  grid <- grid |>
    dplyr::mutate(key = sprintf("dgp%03d_rep%03d", dgp_id, rep_id))

  pending <- grid |>
    dplyr::filter(!(key %in% done_keys))

  cat("Already completed:", nrow(grid) - nrow(pending), "\n")
  cat("Pending:", nrow(pending), "\n\n")

  if (nrow(pending) == 0) {
    return(load_cl2_results(output_dir))
  }

  for (i in seq_len(nrow(pending))) {
    row <- pending[i, , drop = FALSE]
    row_list <- as.list(row)
    save_path <- file.path(output_dir, paste0(row$key, ".rds"))

    cat(sprintf(
      "[%d/%d] dgp=%d rep=%d (n=%d, corr=%s)\n",
      i,
      nrow(pending),
      row_list$dgp_id,
      row_list$rep_id,
      row_list$n,
      row_list$corr_type
    ))

    result <- tryCatch(
      run_one_cl2_setting(row_list, alpha = alpha),
      error = function(e) {
        fail_row <- make_failure_metrics_row(row_list, conditionMessage(e))
        fail_row$cl2_capped_clusters <- NA_integer_
        fail_row
      }
    )

    atomic_saveRDS(result, save_path)
  }

  all_results <- load_cl2_results(output_dir)
  atomic_saveRDS(all_results, file.path(output_dir, "results_combined.rds"))
  all_results
}

summarize_cl2_pilot <- function(results, alpha = CL2_ALPHA) {
  ok <- results |>
    dplyr::filter(converged %in% TRUE, !is.na(method), !is.na(term_type))

  by_cell <- ok |>
    dplyr::group_by(n, corr_type, corr_param, term_type, method) |>
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
    dplyr::filter(method %in% c("pffr_sandwich", "pffr_cl2")) |>
    dplyr::select(
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
    "coverage_pffr_sandwich",
    "coverage_pffr_cl2",
    "mean_width_pffr_sandwich",
    "mean_width_pffr_cl2"
  )
  for (nm in required_cols) {
    if (!nm %in% names(compare_cluster_cl2)) {
      compare_cluster_cl2[[nm]] <- NA_real_
    }
  }

  compare_cluster_cl2 <- compare_cluster_cl2 |>
    dplyr::mutate(
      coverage_gap_cluster = abs(.data$coverage_pffr_sandwich - (1 - alpha)),
      coverage_gap_cl2 = abs(.data$coverage_pffr_cl2 - (1 - alpha)),
      coverage_gap_improvement = coverage_gap_cluster - coverage_gap_cl2,
      width_ratio_cl2_vs_cluster = .data$mean_width_pffr_cl2 /
        .data$mean_width_pffr_sandwich
    )

  overall <- ok |>
    dplyr::group_by(term_type, method) |>
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
    target_coverage = CL2_TARGET_COVERAGE
  )
}

write_cl2_report <- function(summary_obj, path, mode, n_rep) {
  cmp <- summary_obj$compare_cluster_cl2 |>
    dplyr::arrange(dplyr::desc(coverage_gap_improvement))

  overall <- summary_obj$overall |>
    dplyr::filter(
      method %in% c("pffr", "pffr_hc", "pffr_sandwich", "pffr_cl2")
    ) |>
    dplyr::arrange(term_type, method)

  lines <- c(
    "# CL2 Pilot Report",
    "",
    sprintf("- Generated: %s", format(Sys.time(), usetz = TRUE)),
    sprintf("- Mode: %s", mode),
    sprintf("- Replications per DGP: %d", n_rep),
    sprintf("- Target coverage: %.2f", summary_obj$target_coverage),
    "",
    "## Overall (all low-n cells pooled)",
    "```",
    capture.output(print(overall, n = Inf)),
    "```",
    "",
    "## Cluster vs CL2 comparison by cell",
    "```",
    capture.output(print(cmp, n = Inf)),
    "```"
  )

  writeLines(lines, path)
}

# ---------------------------------------------------------------------------
# CLI entrypoint
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0) args[[1]] else "pilot"

n_rep <- switch(
  mode,
  smoke = 2L,
  pilot = 8L,
  full = 20L,
  as.integer(mode)
)
if (is.na(n_rep) || n_rep <= 0) {
  stop("Mode must be one of: smoke, pilot, full, or a positive integer.")
}

base_dir <- "ci-benchmark/cl2-pilot"
results_dir <- file.path(base_dir, "results")

cat("Running CL2 pilot in", mode, "mode (", n_rep, "reps per DGP)\n\n")

results <- run_cl2_pilot(
  n_rep = n_rep,
  seed = CL2_BASE_SEED,
  output_dir = results_dir,
  alpha = CL2_ALPHA
)

cat("\nTotal result rows:", nrow(results), "\n")

if (nrow(results) > 0) {
  summary_obj <- summarize_cl2_pilot(results, alpha = CL2_ALPHA)

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
  write_cl2_report(
    summary_obj,
    file.path(base_dir, "pilot_report.md"),
    mode,
    n_rep
  )

  cat("\nSaved summaries to:", base_dir, "\n")
}
