#!/usr/bin/env Rscript
# ============================================================================
# Quick Test: CL2 for gaulss
# ============================================================================
#
# Small sanity check that the experimental CL2 path works for gaulss fits and
# yields plausible calibration shifts versus default/cluster.
#
# Usage:
#   Rscript ci-benchmark/cl2-pilot/quick-test-gaulss.R [mode]
#   mode: "smoke" (1 rep), "pilot" (2 reps), "full" (5 reps)
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

GAULSS_TEST_BASE_SEED <- 7301L
GAULSS_TEST_ALPHA <- 0.10

gaulss_test_settings <- function() {
  tidyr::crossing(
    n = c(30L, 60L),
    corr_type = c("iid", "ar1")
  ) |>
    dplyr::mutate(
      corr_param = dplyr::if_else(corr_type == "ar1", 0.9, NA_real_),
      dgp_id = dplyr::row_number()
    )
}

run_one_gaulss_case <- function(row, rep_id, alpha = GAULSS_TEST_ALPHA) {
  seed <- GAULSS_TEST_BASE_SEED + 1000L * row$dgp_id + rep_id

  sim <- generate_benchmark_data(
    n = row$n,
    nxgrid = 30L,
    nygrid = 40L,
    snr = 25,
    error_dist = "gaussian",
    corr_type = row$corr_type,
    corr_param = row$corr_param,
    hetero_type = "bump",
    hetero_param = 3.0,
    terms = c("ff", "linear"),
    wiggliness = 5,
    seed = seed
  )

  frml <- build_pffr_formula(sim, sim$s_grid, k_smooth = 12, k_ff = c(12, 12))
  bs_yindex <- list(bs = "ps", k = 12, m = c(2, 1))

  t0 <- Sys.time()
  fit <- pffr(
    frml,
    yind = sim$t_grid,
    data = sim$data,
    family = mgcv::gaulss(),
    bs.yindex = bs_yindex,
    sandwich = "none"
  )
  fit_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  md <- compute_all_metrics(
    fit = fit,
    sim = sim,
    method = "gaulss_default",
    use_sandwich = "none",
    alpha = alpha,
    err_struct = sim$err_struct
  )
  mc <- compute_all_metrics(
    fit = fit,
    sim = sim,
    method = "gaulss_cluster",
    use_sandwich = "cluster",
    alpha = alpha,
    err_struct = sim$err_struct
  )
  mcl2 <- cl2_compute_metrics(
    fit = fit,
    sim = sim,
    alpha = alpha,
    method_name = "gaulss_cl2"
  )

  all <- dplyr::bind_rows(md, mc, mcl2$metrics) |>
    dplyr::filter(term_type %in% c("ff", "linear", "E(Y)")) |>
    dplyr::mutate(
      dgp_id = row$dgp_id,
      rep_id = rep_id,
      seed = seed,
      n = row$n,
      corr_type = row$corr_type,
      corr_param = row$corr_param,
      fit_time = fit_time,
      converged = TRUE,
      cl2_capped_clusters = ifelse(
        method == "gaulss_cl2",
        mcl2$n_capped_clusters,
        NA_integer_
      ),
      error_msg = NA_character_
    )

  all
}

summarize_gaulss_quick <- function(results, alpha = GAULSS_TEST_ALPHA) {
  ok <- results |>
    dplyr::filter(converged %in% TRUE, !is.na(coverage), !is.na(method))

  overall <- ok |>
    dplyr::group_by(method, term_type) |>
    dplyr::summarise(
      coverage = mean(coverage, na.rm = TRUE),
      mean_width = mean(mean_width, na.rm = TRUE),
      z_sd = mean(z_sd, na.rm = TRUE),
      reps = dplyr::n_distinct(rep_id),
      dgps = dplyr::n_distinct(dgp_id),
      .groups = "drop"
    )

  compare <- ok |>
    dplyr::filter(method %in% c("gaulss_cluster", "gaulss_cl2")) |>
    dplyr::group_by(n, corr_type, term_type, method) |>
    dplyr::summarise(
      coverage = mean(coverage, na.rm = TRUE),
      mean_width = mean(mean_width, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = method,
      values_from = c(coverage, mean_width)
    ) |>
    dplyr::mutate(
      coverage_gap_cluster = abs(coverage_gaulss_cluster - (1 - alpha)),
      coverage_gap_cl2 = abs(coverage_gaulss_cl2 - (1 - alpha)),
      coverage_gap_improvement = coverage_gap_cluster - coverage_gap_cl2,
      width_ratio_cl2_vs_cluster = mean_width_gaulss_cl2 /
        mean_width_gaulss_cluster
    )

  list(overall = overall, compare = compare)
}

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0) args[[1]] else "pilot"
n_rep <- switch(
  mode,
  smoke = 1L,
  pilot = 2L,
  full = 5L,
  as.integer(mode)
)
if (is.na(n_rep) || n_rep <= 0) {
  stop("Mode must be one of: smoke, pilot, full, or a positive integer.")
}

settings <- gaulss_test_settings()
grid <- settings |>
  tidyr::crossing(rep_id = seq_len(n_rep))

cat("gaulss quick test:", nrow(settings), "cells x", n_rep, "reps\n")

out <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  row <- grid[i, , drop = FALSE]
  cat(sprintf(
    "[%d/%d] dgp=%d n=%d corr=%s rep=%d\n",
    i,
    nrow(grid),
    row$dgp_id,
    row$n,
    row$corr_type,
    row$rep_id
  ))
  out[[i]] <- tryCatch(
    run_one_gaulss_case(as.list(row), rep_id = row$rep_id),
    error = function(e) {
      tibble(
        method = NA_character_,
        term_type = NA_character_,
        coverage = NA_real_,
        mean_width = NA_real_,
        z_sd = NA_real_,
        dgp_id = row$dgp_id,
        rep_id = row$rep_id,
        seed = NA_integer_,
        n = row$n,
        corr_type = row$corr_type,
        corr_param = row$corr_param,
        fit_time = NA_real_,
        converged = FALSE,
        cl2_capped_clusters = NA_integer_,
        error_msg = conditionMessage(e)
      )
    }
  )
}

results <- dplyr::bind_rows(out)
summary_obj <- summarize_gaulss_quick(results)

base_dir <- "ci-benchmark/cl2-pilot/gaulss-quick"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(results, file.path(base_dir, "results.rds"))
write.csv(results, file.path(base_dir, "results.csv"), row.names = FALSE)
saveRDS(summary_obj, file.path(base_dir, "summary.rds"))
write.csv(
  summary_obj$overall,
  file.path(base_dir, "summary_overall.csv"),
  row.names = FALSE
)
write.csv(
  summary_obj$compare,
  file.path(base_dir, "summary_compare.csv"),
  row.names = FALSE
)

cat("\nSummary (overall):\n")
print(summary_obj$overall, n = Inf)
cat("\nSummary (cluster vs cl2):\n")
print(summary_obj$compare, n = Inf)
