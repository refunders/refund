# ===========================================================================
# Bootstrap CI Extension: curve-resample bootstrap for Study 1 and Study 2
# ===========================================================================
#
# Runs pffr_coefboot() (percentile CIs, curve-resample method) on the SAME
# seeds and fits as Study 1 (non-Gaussian) and Study 2 (grid-refinement).
# One task = one DGP cell, looping over its replicates. Output goes to:
#   ci-benchmark/study1-boot/
#   ci-benchmark/study2-boot/
#
# Usage:
#   Rscript ci-benchmark/sim-study-boot-extension.R <study> <mode> [dgp_id]
#   study: "study1" or "study2"
#   mode:  "smoke" | "pilot" | "full"
#   dgp_id: integer DGP cell (optional; uses SLURM_ARRAY_TASK_ID if unset)
#
# B by mode:  smoke=5, pilot=99, full=499
# Reps:       smoke=1, pilot=10, full= (study1: 150, study2: 50)
#
# SLURM: each array task = one dgp_id, loops its reps sequentially.
# ncpus: reads NCPUS env var (set by SLURM from --cpus-per-task), default 1.
# ===========================================================================

# Setup -----------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(refund)
}

# Study 1 needs mvtnorm; load it so sourcing study1 doesn't fail
suppressPackageStartupMessages(library(mvtnorm))

source("ci-benchmark/benchmark-utils.R")
source("ci-benchmark/confint-benchmark.R")

# Source study scripts for DGP functions and constants (guarded main blocks
# won't execute because sys.nframe() > 0 when sourced).
source("ci-benchmark/sim-study-nongaussian-sandwich.R")
source("ci-benchmark/sim-study-grid-refinement.R")

# Bootstrap CI Extraction -----------------------------------------------------

#' Extract bootstrap CI coverage/width for one term
#'
#' Reuses find_term_index() and evaluate_truth_on_grid() from confint-benchmark.
#' The boot result from pffr_coefboot() has `lower`/`upper` appended to
#' smterms[[i]]$coef — exactly like coef.pffr with ci_from_coef semantics.
#'
#' @param boot_coef Return value of pffr_coefboot().
#' @param truth Truth list from simulation.
#' @param term_type One of "ff", "linear", "smooth", "concurrent", "intercept".
#' @param s_grid s evaluation grid from simulation.
#' @param t_grid t evaluation grid from simulation.
#' @param data Data frame for smooth centering.
#' @param fit Fitted pffr model (for centering decisions).
#' @returns Tibble with coverage, mean_width, n_grid — or NULL if unavailable.
extract_boot_term_metrics <- function(
  boot_coef,
  truth,
  term_type,
  s_grid = NULL,
  t_grid = NULL,
  data = NULL,
  fit = NULL
) {
  if (is.null(boot_coef) || is.null(boot_coef$smterms)) return(NULL)

  sm_names <- names(boot_coef$smterms)
  term_idx <- find_term_index(sm_names, term_type)
  if (is.null(term_idx)) return(NULL)

  term_info <- boot_coef$smterms[[term_idx]]
  if (is.null(term_info) || is.null(term_info$coef)) return(NULL)

  # pffr_coefboot appends lower/upper to smterms[[i]]$coef
  lower <- term_info$coef$lower
  upper <- term_info$coef$upper
  if (is.null(lower) || is.null(upper)) return(NULL)

  est <- term_info$coef$value

  # Evaluate truth on the same evaluation grid as the bootstrap CIs
  truth_vals <- tryCatch(
    evaluate_truth_on_grid(
      truth = truth,
      term_type = term_type,
      term_info = term_info,
      s_grid = s_grid,
      t_grid = t_grid,
      data = data,
      fit = fit,
      center_ff = "match_fit"
    ),
    error = function(e) NULL
  )
  if (is.null(truth_vals)) return(NULL)

  if (length(truth_vals) != length(lower)) {
    warning(sprintf(
      "boot extract: length mismatch for %s: lower=%d, truth=%d",
      term_type,
      length(lower),
      length(truth_vals)
    ))
    return(NULL)
  }

  covered <- (truth_vals >= lower) & (truth_vals <= upper)
  width <- upper - lower

  tibble(
    term_type = term_type,
    coverage = mean(covered, na.rm = TRUE),
    mean_width = mean(width, na.rm = TRUE),
    n_grid = sum(!is.na(covered)),
    method = "boot"
  )
}

#' Extract bootstrap metrics for all relevant terms in one fit
#'
#' @param fit Fitted pffr model.
#' @param sim Simulation result (with truth, s_grid, t_grid, data, terms).
#' @param B Number of bootstrap replicates.
#' @param ncpus Number of CPUs for parallel bootstrap.
#' @returns Tibble with per-term boot metrics, or a single-row NA tibble on error.
extract_boot_metrics <- function(fit, sim, B, ncpus = 1L) {
  # Run bootstrap — percentile CIs, curve resample, conf=0.90
  parallel_mode <- if (ncpus > 1L) "multicore" else "no"

  boot_coef <- tryCatch(
    pffr_coefboot(
      fit,
      B = B,
      n1 = 50L,
      n2 = 25L,
      n3 = 15L,
      conf = 0.90,
      type = "percent",
      method = "resample",
      parallel = parallel_mode,
      ncpus = ncpus,
      showProgress = FALSE
    ),
    error = function(e) {
      warning("pffr_coefboot failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(boot_coef)) return(NULL)

  # Determine which terms to extract: sim$terms + intercept (if present in fit)
  term_types <- unique(c(sim$terms, "intercept"))

  results <- list()
  for (tt in term_types) {
    m <- tryCatch(
      extract_boot_term_metrics(
        boot_coef = boot_coef,
        truth = sim$truth,
        term_type = tt,
        s_grid = sim$s_grid,
        t_grid = sim$t_grid,
        data = sim$data,
        fit = fit
      ),
      error = function(e) {
        warning(sprintf(
          "extract_boot_term_metrics failed for %s: %s",
          tt,
          conditionMessage(e)
        ))
        NULL
      }
    )
    if (!is.null(m)) results[[length(results) + 1]] <- m
  }

  if (length(results) == 0L) return(NULL)
  dplyr::bind_rows(results)
}

# NA row helpers ---------------------------------------------------------------

#' Build a failure row with NAs for the boot metrics columns
#'
#' @param row DGP settings list.
#' @param rep_id Rep ID.
#' @param seed Seed used.
#' @param fit_time Fitting time (or NA).
#' @param error_msg Error message string.
#' @param extra_cols Named list of extra columns (e.g. grid_label for study2).
#' @returns Single-row tibble with NA metric values.
na_boot_row <- function(
  row,
  rep_id,
  seed,
  fit_time = NA_real_,
  error_msg = NA_character_,
  extra_cols = list()
) {
  base <- tibble(
    term_type = NA_character_,
    coverage = NA_real_,
    mean_width = NA_real_,
    n_grid = NA_integer_,
    method = "boot",
    dgp_id = row$dgp_id,
    rep_id = rep_id,
    seed = seed,
    fit_time = fit_time,
    converged = FALSE,
    error_msg = error_msg
  )
  for (nm in names(extra_cols)) base[[nm]] <- extra_cols[[nm]]
  base
}

# Study 1 runner --------------------------------------------------------------

#' Run bootstrap extension for one Study 1 (dgp, rep)
#'
#' Identical seed and fit as Study 1 production; extracts boot CIs only.
#'
#' @param row DGP settings list (from make_study1_settings()).
#' @param rep_id Replicate number.
#' @param B Bootstrap replicates.
#' @param ncpus CPUs for pffr_coefboot().
#' @returns Annotated metrics tibble.
run_study1_boot_rep <- function(row, rep_id, B, ncpus = 1L) {
  if (inherits(row, "data.frame")) row <- as.list(row)

  seed <- STUDY1_BASE_SEED + 1000L * row$dgp_id + rep_id

  sim_fn <- if (row$family == "poisson") simulate_poisson_ff_linear else
    simulate_binomial_ff_linear

  sim <- sim_fn(
    n = row$n,
    nxgrid = row$nxgrid,
    nygrid = row$nygrid,
    corr_type = row$corr_type,
    corr_param = if (is.na(row$corr_param)) NULL else row$corr_param,
    seed = seed
  )

  frml <- build_study1_formula(sim$s_grid)
  fam <- if (row$family == "poisson") poisson() else binomial()

  # Store data and yind in the formula environment so pffr_coefboot can
  # resolve them. pffr_coefboot evaluates modcall$data and modcall$yind in
  # frml_env for each bootstrap replicate. sim$data/t_grid are not in frml_env.
  # bs.yindex is passed as a literal list() call so match.call() captures it
  # as an evaluable expression (not a symbol), making eval(modcall_boot) safe.
  pffr_data <- sim$data
  pffr_yind <- sim$t_grid
  frml_env <- environment(frml)
  assign("pffr_data", pffr_data, envir = frml_env)
  assign("pffr_yind", pffr_yind, envir = frml_env)

  t0 <- Sys.time()
  fit <- pffr(
    frml,
    yind = pffr_yind,
    data = pffr_data,
    family = fam,
    bs.yindex = list(bs = "ps", k = STUDY1_K_YINDEX, m = c(2, 1)),
    sandwich = "none"
  )
  fit_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  metrics <- extract_boot_metrics(fit, sim, B = B, ncpus = ncpus)

  extra <- list(
    family = row$family,
    n = row$n,
    nxgrid = row$nxgrid,
    nygrid = row$nygrid,
    corr_type = row$corr_type,
    corr_param = row$corr_param
  )

  if (is.null(metrics) || nrow(metrics) == 0L) {
    row_out <- na_boot_row(
      row,
      rep_id,
      seed,
      fit_time,
      "boot metrics extraction returned 0 rows",
      extra_cols = extra
    )
    return(row_out)
  }

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

#' Run Study 1 bootstrap extension for one DGP cell (all reps)
#'
#' Atomic per-(dgp,rep) saves; skip-existing resume.
#'
#' @param dgp_row One-row tibble from make_study1_settings().
#' @param n_rep Number of reps.
#' @param B Bootstrap replicates.
#' @param ncpus CPUs for pffr_coefboot().
#' @param output_dir Output directory.
run_study1_boot_cell <- function(
  dgp_row,
  n_rep,
  B,
  ncpus = 1L,
  output_dir = "ci-benchmark/study1-boot"
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  row <- as.list(dgp_row)

  cat(sprintf(
    "Study 1 Boot | DGP %d (family=%s, corr=%s, n=%d) | B=%d | reps=1..%d\n",
    row$dgp_id,
    row$family,
    row$corr_type,
    row$n,
    B,
    n_rep
  ))

  for (rep_id in seq_len(n_rep)) {
    save_path <- file.path(
      output_dir,
      sprintf("dgp%03d_rep%03d.rds", row$dgp_id, rep_id)
    )
    if (file.exists(save_path)) {
      cat(sprintf("  rep %d: already done, skipping\n", rep_id))
      next
    }

    cat(sprintf("  rep %d ...", rep_id))
    t_rep <- Sys.time()

    result <- tryCatch(
      run_study1_boot_rep(row, rep_id, B = B, ncpus = ncpus),
      error = function(e) {
        warning(sprintf(
          "run_study1_boot_rep failed for dgp=%d rep=%d: %s",
          row$dgp_id,
          rep_id,
          conditionMessage(e)
        ))
        extra <- list(
          family = row$family,
          n = row$n,
          nxgrid = row$nxgrid,
          nygrid = row$nygrid,
          corr_type = row$corr_type,
          corr_param = row$corr_param
        )
        seed <- STUDY1_BASE_SEED + 1000L * row$dgp_id + rep_id
        na_boot_row(row, rep_id, seed, NA_real_, conditionMessage(e), extra)
      }
    )

    atomic_saveRDS(result, save_path)
    elapsed <- as.numeric(difftime(Sys.time(), t_rep, units = "secs"))
    cat(sprintf(" done (%.1fs)\n", elapsed))
    gc()
  }
}

# Study 2 runner --------------------------------------------------------------

#' Run bootstrap extension for one Study 2 (dgp, rep) across all nygrid levels
#'
#' Generates paired data once (finest grid), subsamples to each nygrid level,
#' then runs bootstrap for each. Atomic per-(dgp,grid,rep) saves.
#'
#' @param row DGP settings list (from make_study2_settings()).
#' @param rep_id Replicate number.
#' @param B Bootstrap replicates.
#' @param ncpus CPUs for pffr_coefboot().
#' @param output_dir Output directory.
run_study2_boot_pair <- function(
  row,
  rep_id,
  B,
  ncpus = 1L,
  output_dir = "ci-benchmark/study2-boot"
) {
  if (inherits(row, "data.frame")) row <- as.list(row)

  seed <- STUDY2_BASE_SEED + 1000L * row$dgp_id + rep_id

  sims <- generate_paired_data(row, seed)

  grid_labels <- names(STUDY2_GRIDS)

  for (grid_label in grid_labels) {
    nygrid <- STUDY2_GRIDS[[grid_label]]$nygrid
    nxgrid <- STUDY2_GRIDS[[grid_label]]$nxgrid

    save_key <- sprintf(
      "dgp%03d_n%03d_y%03d_rep%03d",
      row$dgp_id,
      row$n,
      nygrid,
      rep_id
    )
    save_path <- file.path(output_dir, paste0(save_key, ".rds"))

    if (file.exists(save_path)) next

    sim <- sims[[grid_label]]

    frml <- build_pffr_formula(
      sim,
      sim$s_grid,
      k_smooth = 12,
      k_ff = c(12L, 12L)
    )

    extra <- list(
      corr_type = row$corr_type,
      corr_param = row$corr_param,
      n = row$n,
      snr = row$snr,
      nxgrid = nxgrid,
      nygrid = nygrid,
      grid_label = grid_label
    )

    # Store data and yind in the formula environment so pffr_coefboot can
    # resolve them. bs.yindex is passed as a literal list() expression so
    # that eval(modcall_boot) in pffr_coefboot can evaluate it without needing
    # a symbol lookup (list() is always available, unlike local variables).
    pffr_data <- sim$data
    pffr_yind <- sim$t_grid
    frml_env <- environment(frml)
    assign("pffr_data", pffr_data, envir = frml_env)
    assign("pffr_yind", pffr_yind, envir = frml_env)

    t0 <- Sys.time()
    fit <- tryCatch(
      pffr(
        frml,
        yind = pffr_yind,
        data = pffr_data,
        bs.yindex = list(bs = "ps", k = 12L, m = c(2L, 1L)),
        sandwich = "none"
      ),
      error = function(e) {
        warning(sprintf(
          "pffr failed for dgp=%d rep=%d grid=%s: %s",
          row$dgp_id,
          rep_id,
          grid_label,
          conditionMessage(e)
        ))
        NULL
      }
    )
    fit_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (is.null(fit)) {
      result <- na_boot_row(
        row,
        rep_id,
        seed,
        fit_time,
        "pffr fit failed",
        extra_cols = extra
      )
      atomic_saveRDS(result, save_path)
      next
    }

    metrics <- extract_boot_metrics(fit, sim, B = B, ncpus = ncpus)

    if (is.null(metrics) || nrow(metrics) == 0L) {
      result <- na_boot_row(
        row,
        rep_id,
        seed,
        fit_time,
        "boot metrics extraction returned 0 rows",
        extra_cols = extra
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
          nxgrid = nxgrid,
          nygrid = nygrid,
          grid_label = grid_label,
          fit_time = fit_time,
          converged = TRUE,
          error_msg = NA_character_
        )
    }

    atomic_saveRDS(result, save_path)
    gc()
  }
}

#' Run Study 2 bootstrap extension for one DGP cell (all reps)
#'
#' @param dgp_row One-row tibble from make_study2_settings().
#' @param n_rep Number of reps.
#' @param B Bootstrap replicates.
#' @param ncpus CPUs for pffr_coefboot().
#' @param output_dir Output directory.
run_study2_boot_cell <- function(
  dgp_row,
  n_rep,
  B,
  ncpus = 1L,
  output_dir = "ci-benchmark/study2-boot"
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  row <- as.list(dgp_row)

  cat(sprintf(
    "Study 2 Boot | DGP %d (corr=%s, n=%d, snr=%g) | B=%d | reps=1..%d\n",
    row$dgp_id,
    row$corr_type,
    row$n,
    row$snr,
    B,
    n_rep
  ))

  for (rep_id in seq_len(n_rep)) {
    # Check whether ALL grid levels for this (dgp, rep) are already done
    grid_labels <- names(STUDY2_GRIDS)
    all_done <- all(vapply(
      grid_labels,
      function(gl) {
        ny <- STUDY2_GRIDS[[gl]]$nygrid
        save_key <- sprintf(
          "dgp%03d_n%03d_y%03d_rep%03d",
          row$dgp_id,
          row$n,
          ny,
          rep_id
        )
        file.exists(file.path(output_dir, paste0(save_key, ".rds")))
      },
      logical(1)
    ))

    if (all_done) {
      cat(sprintf("  rep %d: all grids done, skipping\n", rep_id))
      next
    }

    cat(sprintf("  rep %d ...", rep_id))
    t_rep <- Sys.time()

    tryCatch(
      run_study2_boot_pair(
        row,
        rep_id,
        B = B,
        ncpus = ncpus,
        output_dir = output_dir
      ),
      error = function(e) {
        warning(sprintf(
          "run_study2_boot_pair failed for dgp=%d rep=%d: %s",
          row$dgp_id,
          rep_id,
          conditionMessage(e)
        ))
        seed <- STUDY2_BASE_SEED + 1000L * row$dgp_id + rep_id
        extra <- list(
          corr_type = row$corr_type,
          corr_param = row$corr_param,
          n = row$n,
          snr = row$snr,
          nxgrid = NA_integer_,
          nygrid = NA_integer_,
          grid_label = NA_character_
        )
        result <- na_boot_row(
          row,
          rep_id,
          seed,
          NA_real_,
          conditionMessage(e),
          extra
        )
        # Save catch-all failure to the finest-grid key so the loader
        # (which matches ^dgp\d+_n\d+_y\d+_rep\d+\.rds$) picks it up.
        ny_fine <- max(vapply(STUDY2_GRIDS, `[[`, integer(1), "nygrid"))
        save_key <- sprintf(
          "dgp%03d_n%03d_y%03d_rep%03d",
          row$dgp_id,
          row$n,
          ny_fine,
          rep_id
        )
        save_path <- file.path(output_dir, paste0(save_key, ".rds"))
        atomic_saveRDS(result, save_path)
      }
    )

    elapsed <- as.numeric(difftime(Sys.time(), t_rep, units = "secs"))
    cat(sprintf(" done (%.1fs)\n", elapsed))
  }
}

# Load Results ----------------------------------------------------------------

#' Load Study 1 boot results from incremental saves
load_study1_boot_results <- function(
  output_dir = "ci-benchmark/study1-boot"
) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(all_files) == 0L) {
    warning("No Study 1 boot result files in ", output_dir)
    return(tibble())
  }
  dplyr::bind_rows(lapply(all_files, readRDS))
}

#' Load Study 2 boot results from incremental saves
load_study2_boot_results <- function(
  output_dir = "ci-benchmark/study2-boot"
) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_n\\d+_y\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(all_files) == 0L) {
    warning("No Study 2 boot result files in ", output_dir)
    return(tibble())
  }
  dplyr::bind_rows(lapply(all_files, readRDS))
}

# Main Entry Point ------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  study <- if (length(args) >= 1L) args[1] else "study2"
  mode <- if (length(args) >= 2L) args[2] else "smoke"
  # Optional dgp_id from args[3] or SLURM_ARRAY_TASK_ID
  dgp_arg <- if (length(args) >= 3L) as.integer(args[3]) else NA_integer_
  slurm_tid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  dgp_filter <- if (!is.na(dgp_arg)) {
    dgp_arg
  } else if (nchar(slurm_tid) > 0L) {
    as.integer(slurm_tid)
  } else {
    NA_integer_ # run all dgp cells (sequential, for local smoke)
  }

  # Number of CPUs: from NCPUS env or SLURM_CPUS_PER_TASK, default 1
  ncpus_env <- Sys.getenv("NCPUS", unset = "")
  ncpus_slurm <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  ncpus <- if (nchar(ncpus_env) > 0L) {
    max(1L, as.integer(ncpus_env))
  } else if (nchar(ncpus_slurm) > 0L) {
    max(1L, as.integer(ncpus_slurm))
  } else {
    1L
  }

  # B and n_rep by mode
  B <- switch(mode, smoke = 5L, pilot = 99L, full = 499L, {
    warning("Unknown mode '", mode, "', defaulting to smoke")
    5L
  })
  if (study == "study1") {
    n_rep <- switch(mode, smoke = 1L, pilot = 10L, full = 150L, 1L)
    output_dir <- "ci-benchmark/study1-boot"
  } else if (study == "study2") {
    n_rep <- switch(mode, smoke = 1L, pilot = 10L, full = 50L, 1L)
    output_dir <- "ci-benchmark/study2-boot"
  } else {
    stop("study must be 'study1' or 'study2', got: ", study)
  }

  cat("=================================================\n")
  cat("Bootstrap CI Extension\n")
  cat("  study      :", study, "\n")
  cat("  mode       :", mode, "\n")
  cat("  B          :", B, "\n")
  cat("  n_rep      :", n_rep, "\n")
  cat("  ncpus      :", ncpus, "\n")
  cat("  dgp_filter :", if (is.na(dgp_filter)) "all" else dgp_filter, "\n")
  cat("  output_dir :", output_dir, "\n")
  cat("=================================================\n\n")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Dispatch
  if (study == "study1") {
    settings <- make_study1_settings()

    if (!is.na(dgp_filter)) {
      settings <- dplyr::filter(settings, dgp_id == dgp_filter)
      if (nrow(settings) == 0L) {
        stop("No Study 1 DGP found with dgp_id = ", dgp_filter)
      }
    }

    for (i in seq_len(nrow(settings))) {
      run_study1_boot_cell(
        dgp_row = settings[i, ],
        n_rep = n_rep,
        B = B,
        ncpus = ncpus,
        output_dir = output_dir
      )
    }

    # Print summary if in smoke/pilot mode
    if (mode %in% c("smoke", "pilot")) {
      results <- load_study1_boot_results(output_dir)
      if (nrow(results) > 0L) {
        cat("\n========== STUDY 1 BOOT COVERAGE SUMMARY ==========\n")
        summary_df <- results |>
          dplyr::filter(!is.na(coverage)) |>
          dplyr::group_by(family, corr_type, n, term_type) |>
          dplyr::summarise(
            mean_coverage = mean(coverage, na.rm = TRUE),
            mean_width = mean(mean_width, na.rm = TRUE),
            n_reps = dplyr::n(),
            .groups = "drop"
          )
        print(summary_df, n = 100)
      }
    }
  } else {
    # study2
    settings <- make_study2_settings()

    if (!is.na(dgp_filter)) {
      settings <- dplyr::filter(settings, dgp_id == dgp_filter)
      if (nrow(settings) == 0L) {
        stop("No Study 2 DGP found with dgp_id = ", dgp_filter)
      }
    }

    for (i in seq_len(nrow(settings))) {
      run_study2_boot_cell(
        dgp_row = settings[i, ],
        n_rep = n_rep,
        B = B,
        ncpus = ncpus,
        output_dir = output_dir
      )
    }

    if (mode %in% c("smoke", "pilot")) {
      results <- load_study2_boot_results(output_dir)
      if (nrow(results) > 0L) {
        cat("\n========== STUDY 2 BOOT COVERAGE SUMMARY ==========\n")
        summary_df <- results |>
          dplyr::filter(!is.na(coverage)) |>
          dplyr::group_by(corr_type, n, snr, grid_label, term_type) |>
          dplyr::summarise(
            mean_coverage = mean(coverage, na.rm = TRUE),
            mean_width = mean(mean_width, na.rm = TRUE),
            n_reps = dplyr::n(),
            .groups = "drop"
          )
        print(summary_df, n = 100)
      }
    }
  }

  cat("\nBootstrap extension complete.\n")
}
