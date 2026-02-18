#!/usr/bin/env Rscript
# ============================================================================
# Study 1 Follow-up: Factorial n x grid Fill-in
# ============================================================================
#
# Motivation:
#   Study 1 can be hard to interpret if n and grid were varied in paired
#   combinations only (diagonal design). This script fills missing factorial
#   n x grid cells while reusing the existing Study 1 simulation + metric code.
#
# Usage:
#   Rscript ci-benchmark/run-study1-factorial-fillin.R [mode] [design_mode]
#     mode:
#       - "smoke" (1 rep, capped task count)
#       - "pilot" (10 reps)
#       - "full"  (150 reps)
#       - or any positive integer reps
#     design_mode:
#       - "missing_diagonal" (default): run off-diagonal n x grid combinations
#         assuming diagonal pairs were already run.
#       - "missing_from_existing": infer completed n x grid pairs from an
#         existing results directory and run only missing pairs.
#       - "full_factorial": run all n x grid pairs.
#
# Notes:
#   - Reuses run_study1_rep() and summarize_study1() from
#     ci-benchmark/sim-study-nongaussian-sandwich.R.
#   - Uses incremental per-task RDS saves and resume-safe skipping.
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

STUDY1F_DEFAULT_N_LEVELS <- c(100L, 200L, 400L)
STUDY1F_DEFAULT_GRIDS <- tibble::tibble(
  grid_label = c("coarse", "medium", "fine"),
  nxgrid = c(30L, 60L, 90L),
  nygrid = c(40L, 80L, 120L)
)
STUDY1F_DEFAULT_FAMILIES <- c("poisson", "binomial")
STUDY1F_DEFAULT_CORR_TYPES <- c("iid", "ar1", "fourier_pos")
STUDY1F_DEFAULT_ALPHA <- 0.10
STUDY1F_DEFAULT_OUTPUT_DIR <- "ci-benchmark/study1-nongaussian-factorial-fillin"
STUDY1F_EXISTING_DIR <- "ci-benchmark/study1-nongaussian"
STUDY1F_SMOKE_MAX_TASKS <- 4L

build_study1_factorial_settings <- function(
  n_levels = STUDY1F_DEFAULT_N_LEVELS,
  grids = STUDY1F_DEFAULT_GRIDS,
  families = STUDY1F_DEFAULT_FAMILIES,
  corr_types = STUDY1F_DEFAULT_CORR_TYPES
) {
  tidyr::crossing(
    family = families,
    corr_type = corr_types,
    n = as.integer(n_levels),
    grids
  ) |>
    dplyr::mutate(
      corr_param = dplyr::case_when(
        corr_type == "ar1" ~ 0.9,
        corr_type == "fourier_pos" ~ 0.3,
        TRUE ~ NA_real_
      ),
      dgp_id = dplyr::row_number()
    )
}

make_diagonal_pairs <- function(n_levels, grids) {
  n_levels <- as.integer(n_levels)
  if (length(n_levels) != nrow(grids)) {
    stop(
      "missing_diagonal requires length(n_levels) == nrow(grids). ",
      "Got ",
      length(n_levels),
      " n levels and ",
      nrow(grids),
      " grids."
    )
  }

  tibble::tibble(
    n = n_levels,
    nxgrid = grids$nxgrid,
    nygrid = grids$nygrid
  )
}

extract_observed_pairs <- function(
  existing_results_dir = STUDY1F_EXISTING_DIR
) {
  existing <- tryCatch(
    load_study1_results(existing_results_dir),
    warning = function(w) tibble::tibble(),
    error = function(e) tibble::tibble()
  )

  if (!nrow(existing)) {
    return(tibble::tibble(
      n = integer(),
      nxgrid = integer(),
      nygrid = integer()
    ))
  }

  existing |>
    dplyr::filter(!is.na(n), !is.na(nxgrid), !is.na(nygrid)) |>
    dplyr::distinct(n, nxgrid, nygrid)
}

select_settings_for_mode <- function(
  full_settings,
  design_mode = c(
    "missing_diagonal",
    "missing_from_existing",
    "full_factorial"
  ),
  n_levels = STUDY1F_DEFAULT_N_LEVELS,
  grids = STUDY1F_DEFAULT_GRIDS,
  existing_results_dir = STUDY1F_EXISTING_DIR
) {
  design_mode <- match.arg(design_mode)

  if (design_mode == "full_factorial") {
    return(full_settings)
  }

  observed_pairs <- if (design_mode == "missing_diagonal") {
    make_diagonal_pairs(n_levels, grids)
  } else {
    extract_observed_pairs(existing_results_dir)
  }

  full_settings |>
    dplyr::anti_join(observed_pairs, by = c("n", "nxgrid", "nygrid"))
}

load_factorial_results <- function(output_dir) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(all_files)) return(tibble::tibble())
  dplyr::bind_rows(lapply(all_files, readRDS))
}

run_study1_factorial <- function(
  settings,
  n_rep,
  output_dir = STUDY1F_DEFAULT_OUTPUT_DIR,
  alpha = STUDY1F_DEFAULT_ALPHA,
  smoke_max_tasks = NULL
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (!nrow(settings)) {
    cat("No settings selected for this design mode.\n")
    return(load_factorial_results(output_dir))
  }

  grid <- settings |>
    tidyr::crossing(rep_id = seq_len(n_rep))

  if (!is.null(smoke_max_tasks) && nrow(grid) > smoke_max_tasks) {
    # Prefer cheapest cells for smoke tests.
    grid <- grid |>
      dplyr::mutate(smoke_cost = n * nxgrid * nygrid) |>
      dplyr::arrange(smoke_cost, family, corr_type, rep_id) |>
      dplyr::slice_head(n = smoke_max_tasks) |>
      dplyr::select(-smoke_cost)
  }

  cat("Study 1 factorial follow-up\n")
  cat("===========================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("Replications:", n_rep, "\n")
  cat("Total tasks:", nrow(grid), "\n")
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

  if (!nrow(pending)) {
    return(load_factorial_results(output_dir))
  }

  for (i in seq_len(nrow(pending))) {
    row <- pending[i, , drop = FALSE]
    row_list <- as.list(row)

    cat(sprintf(
      "[%d/%d] dgp=%d family=%s corr=%s n=%d grid=%dx%d rep=%d\n",
      i,
      nrow(pending),
      row_list$dgp_id,
      row_list$family,
      row_list$corr_type,
      row_list$n,
      row_list$nxgrid,
      row_list$nygrid,
      row_list$rep_id
    ))

    res <- tryCatch(
      run_study1_rep(row_list, rep_id = row_list$rep_id, alpha = alpha),
      error = function(e) {
        tibble::tibble(
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
          method = NA_character_,
          dgp_id = row_list$dgp_id,
          rep_id = row_list$rep_id,
          seed = NA_integer_,
          family = row_list$family,
          n = row_list$n,
          nxgrid = row_list$nxgrid,
          nygrid = row_list$nygrid,
          corr_type = row_list$corr_type,
          corr_param = row_list$corr_param,
          fit_time = NA_real_,
          converged = FALSE,
          error_msg = conditionMessage(e)
        )
      }
    )

    save_path <- file.path(output_dir, paste0(row_list$key, ".rds"))
    atomic_saveRDS(res, save_path)
  }

  all_results <- load_factorial_results(output_dir)
  atomic_saveRDS(all_results, file.path(output_dir, "results_combined.rds"))
  all_results
}

save_factorial_outputs <- function(results, settings, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (nrow(settings)) {
    settings_out <- settings |>
      dplyr::arrange(family, corr_type, n, nxgrid, nygrid)
    write.csv(
      settings_out,
      file.path(output_dir, "settings.csv"),
      row.names = FALSE
    )
  }

  if (!nrow(results)) {
    cat("No results to summarize.\n")
    return(invisible(NULL))
  }

  summary_df <- summarize_study1(results)
  write.csv(summary_df, file.path(output_dir, "summary.csv"), row.names = FALSE)
  cat("Saved summary to ", file.path(output_dir, "summary.csv"), "\n", sep = "")
}

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[[1]] else "smoke"
design_mode <- if (length(args) >= 2) args[[2]] else "missing_diagonal"

n_rep <- switch(
  mode,
  smoke = 1L,
  pilot = 10L,
  full = 150L,
  as.integer(mode)
)
if (is.na(n_rep) || n_rep <= 0) {
  stop("Mode must be one of: smoke, pilot, full, or a positive integer.")
}

if (
  !design_mode %in%
    c("missing_diagonal", "missing_from_existing", "full_factorial")
) {
  stop(
    "design_mode must be one of: missing_diagonal, missing_from_existing, full_factorial"
  )
}

full_settings <- build_study1_factorial_settings()
settings <- select_settings_for_mode(
  full_settings = full_settings,
  design_mode = design_mode,
  n_levels = STUDY1F_DEFAULT_N_LEVELS,
  grids = STUDY1F_DEFAULT_GRIDS,
  existing_results_dir = STUDY1F_EXISTING_DIR
)

cat("Mode: ", mode, " (n_rep=", n_rep, ")\n", sep = "")
cat("Design mode: ", design_mode, "\n", sep = "")
cat("Selected n x grid pairs:\n")
print(
  settings |>
    dplyr::distinct(n, nxgrid, nygrid) |>
    dplyr::arrange(n, nxgrid, nygrid),
  n = Inf
)
cat("\n")

smoke_cap <- if (mode == "smoke") STUDY1F_SMOKE_MAX_TASKS else NULL

results <- run_study1_factorial(
  settings = settings,
  n_rep = n_rep,
  output_dir = STUDY1F_DEFAULT_OUTPUT_DIR,
  alpha = STUDY1F_DEFAULT_ALPHA,
  smoke_max_tasks = smoke_cap
)

cat("\nTotal rows:", nrow(results), "\n")
save_factorial_outputs(
  results = results,
  settings = settings,
  output_dir = STUDY1F_DEFAULT_OUTPUT_DIR
)
