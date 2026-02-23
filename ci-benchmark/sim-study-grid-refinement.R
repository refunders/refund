# ===========================================================================
# Study 2: Grid Refinement and Sandwich Covariance Quality
# ===========================================================================
#
# Quantifies how CI quality and runtime change with sample size (n) and
# response/covariate observation grid resolution, including CL2 correction.
#
# Design:
#   - Gaussian family (fixed)
#   - 3 correlation structures: IID, AR1(0.9), fourier_pos(0.3)
#   - 3 sample sizes: n = 20, 40, 80
#   - 2 SNR levels: 25 (high), 3 (low)
#   - Full model: ff + linear + smooth + concurrent
#   - Grid: nxgrid in {30,60,90} x nygrid in {40,80,120} (9 combos, factorial)
#   - Factorial DGP design: correlation x n x snr; grid factorial and paired
#   - Paired within each (correlation, n, snr, rep): generate finest grid once,
#     deterministic subsampling for coarser grids
#
# Phases:
#   Pilot: 30 reps per DGP x grid for timing + crude coverage
#   Main:  50 reps per DGP x grid (full factorial comparison)
#
# Methods: default (none), HC, cluster, and CL2 sandwich.
#
# Usage:
#   Rscript ci-benchmark/sim-study-grid-refinement.R [mode]
#   mode: "smoke" (2 reps), "pilot" (30 reps), "main" (120 reps)
# ===========================================================================

# Setup -----------------------------------------------------------------------

library(tidyverse)

if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(refund)
}

source("ci-benchmark/benchmark-utils.R")
source("ci-benchmark/confint-benchmark.R")

# Constants -------------------------------------------------------------------

STUDY2_BASE_SEED <- 4001L
STUDY2_N_LEVELS <- c(20L, 40L, 80L)
STUDY2_SNR_LEVELS <- c(25, 3)
STUDY2_WIGGLINESS <- 5
STUDY2_TERMS <- c("ff", "linear", "smooth", "concurrent")

# Grid levels: full factorial nxgrid x nygrid, subsample from finest
STUDY2_NXGRID_LEVELS <- c(30L, 60L, 90L)
STUDY2_NYGRID_LEVELS <- c(40L, 80L, 120L)
STUDY2_FINE_NX <- max(STUDY2_NXGRID_LEVELS)
STUDY2_FINE_NY <- max(STUDY2_NYGRID_LEVELS)

STUDY2_GRIDS <- local({
  combos <- expand.grid(
    nxgrid = STUDY2_NXGRID_LEVELS,
    nygrid = STUDY2_NYGRID_LEVELS
  )
  grids <- lapply(seq_len(nrow(combos)), function(i) {
    list(nxgrid = combos$nxgrid[i], nygrid = combos$nygrid[i])
  })
  names(grids) <- sprintf("%dx%d", combos$nxgrid, combos$nygrid)
  grids
})

# Grid Subsampling ------------------------------------------------------------

#' Deterministically subsample functional data from fine grid to coarser grid
#'
#' Selects evenly-spaced indices from the fine grid to produce the target
#' resolution. Ensures truly paired data across grid levels.
#'
#' @param sim Simulation result from generate_benchmark_data() on fine grid.
#' @param target_nxgrid Target covariate grid size.
#' @param target_nygrid Target response grid size.
#' @returns Modified sim list with subsampled data and truth.
subsample_grid <- function(sim, target_nxgrid, target_nygrid) {
  orig_nx <- length(sim$s_grid)
  orig_ny <- length(sim$t_grid)

  stopifnot(target_nxgrid <= orig_nx, target_nygrid <= orig_ny)

  # Select evenly-spaced indices
  s_idx <- round(seq(1, orig_nx, length.out = target_nxgrid))
  t_idx <- round(seq(1, orig_ny, length.out = target_nygrid))

  # Ensure unique indices (can happen with rounding)
  s_idx <- unique(s_idx)
  t_idx <- unique(t_idx)
  stopifnot(length(s_idx) == target_nxgrid, length(t_idx) == target_nygrid)

  new_s_grid <- sim$s_grid[s_idx]
  new_t_grid <- sim$t_grid[t_idx]

  # Subsample data matrices
  new_data <- sim$data
  new_data$Y <- I(as.matrix(sim$data$Y)[, t_idx])
  new_data$X1 <- I(as.matrix(sim$data$X1)[, s_idx])
  if (!is.null(sim$data$Xconc)) {
    new_data$Xconc <- I(as.matrix(sim$data$Xconc)[, t_idx])
  }
  # Scalar covariates (zlin, zsmoo) unchanged

  # Subsample truth
  new_truth <- sim$truth

  # Beta components
  beta <- new_truth$beta
  if (!is.null(beta$`ff(X1)`) && is.matrix(beta$`ff(X1)`)) {
    new_truth$beta$`ff(X1)` <- beta$`ff(X1)`[s_idx, t_idx]
  }
  if (!is.null(beta$zlin) && !is.function(beta$zlin)) {
    new_truth$beta$zlin <- beta$zlin[t_idx]
  }
  if (!is.null(beta$Xconc) && !is.function(beta$Xconc)) {
    new_truth$beta$Xconc <- beta$Xconc[t_idx]
  }
  if (!is.null(beta$intercept) && !is.function(beta$intercept)) {
    new_truth$beta$intercept <- beta$intercept[t_idx]
  }
  # s(zsmoo) truth is a function — no subsampling needed

  # Eta and etaTerms
  if (!is.null(new_truth$eta) && is.matrix(new_truth$eta)) {
    new_truth$eta <- new_truth$eta[, t_idx]
  }
  if (!is.null(new_truth$etaTerms)) {
    new_truth$etaTerms <- lapply(new_truth$etaTerms, function(m) {
      if (is.matrix(m)) m[, t_idx, drop = FALSE] else m
    })
  }
  if (!is.null(new_truth$epsilon) && is.matrix(new_truth$epsilon)) {
    new_truth$epsilon <- new_truth$epsilon[, t_idx]
  }

  # Rebuild error structure for new t_grid
  new_err_struct <- make_error_structure(
    t_grid = new_t_grid,
    corr_type = sim$err_struct$corr_type,
    corr_param = sim$err_struct$corr_param,
    hetero_type = sim$err_struct$hetero_type,
    hetero_param = sim$err_struct$hetero_param
  )

  list(
    data = new_data,
    s_grid = new_s_grid,
    t_grid = new_t_grid,
    truth = new_truth,
    err_struct = new_err_struct,
    terms = sim$terms,
    # Metadata for validation
    s_idx = s_idx,
    t_idx = t_idx
  )
}

# Covariance Quality Metrics --------------------------------------------------

#' Extract coefficient vector, SE vector, and truth from coef.pffr for one term
#'
#' @param fit Fitted pffr model.
#' @param term_type One of "ff", "linear", "smooth", "concurrent", "intercept".
#' @param sandwich_type One of "none", "hc", "cluster", or "cl2".
#' @param truth Truth list from simulation (for truth centering in cov quality).
#' @param s_grid s grid from simulation.
#' @param t_grid t grid from simulation.
#' @param data Model data frame (for concurrent term grid evaluation).
#' @returns List with est, se, truth_vals vectors (or NULL).
extract_coef_and_se <- function(
  fit,
  term_type,
  sandwich_type = "none",
  truth = NULL,
  s_grid = NULL,
  t_grid = NULL,
  data = NULL
) {
  sandwich_type <- normalize_sandwich_type(sandwich_type)

  cc <- tryCatch(
    coef(
      fit,
      sandwich = sandwich_type,
      seWithMean = FALSE,
      n1 = 50,
      n2 = 25,
      n3 = 15
    ),
    error = function(e) NULL
  )
  if (is.null(cc) || is.null(cc$smterms)) return(NULL)

  sm_names <- names(cc$smterms)
  term_idx <- find_term_index(sm_names, term_type)
  if (is.null(term_idx)) return(NULL)

  term_info <- cc$smterms[[term_idx]]
  if (is.null(term_info) || is.null(term_info$coef)) return(NULL)

  # Evaluate truth on same grid for centering (R1-H3)
  truth_vals <- NULL
  if (!is.null(truth) && !is.null(s_grid) && !is.null(t_grid)) {
    truth_vals <- tryCatch(
      as.vector(evaluate_truth_on_grid(
        truth,
        term_type,
        term_info,
        s_grid,
        t_grid,
        data = data,
        fit = fit
      )),
      error = function(e) NULL
    )
  }

  list(
    est = term_info$coef$value,
    se = term_info$coef$se,
    truth_vals = truth_vals,
    x = term_info$x,
    y = term_info$y
  )
}

#' Compute covariance quality metrics from collected vectors
#'
#' Compares model-estimated SEs (averaged across reps) to MC covariance of
#' estimation errors. When truth varies across reps (random truth functions),
#' raw estimate variance includes both estimation + truth variation. Centering
#' by truth isolates the estimation variance that sandwich targets. (R1-H3)
#'
#' @param est_matrix Matrix of coefficient estimates (n_rep x n_coef).
#' @param se_matrix Matrix of model SEs (n_rep x n_coef).
#' @param term_type Term identifier for labeling.
#' @param truth_matrix Optional matrix of truth values (n_rep x n_coef).
#'   If provided, MC variance is computed from (est - truth) instead of est.
#' @returns Tibble with covariance quality metrics.
compute_cov_quality_metrics <- function(
  est_matrix,
  se_matrix,
  term_type,
  truth_matrix = NULL
) {
  n_rep <- nrow(est_matrix)
  n_coef <- ncol(est_matrix)

  if (n_rep < 10 || n_coef < 2) {
    return(tibble(
      term_type = term_type,
      n_coef = n_coef,
      n_rep = n_rep,
      diag_ratio_median = NA_real_,
      diag_ratio_iqr_low = NA_real_,
      diag_ratio_iqr_high = NA_real_,
      frobenius_rel_error = NA_real_,
      leading_eigenvalue_ratio = NA_real_,
      truth_centered = !is.null(truth_matrix)
    ))
  }

  # Center by truth if available (R1-H3: removes truth variation from MC var)
  error_matrix <- if (
    !is.null(truth_matrix) &&
      identical(dim(truth_matrix), dim(est_matrix))
  ) {
    est_matrix - truth_matrix
  } else {
    est_matrix
  }

  # Diagonal comparison: mean(se^2) vs var(est - truth)
  mean_model_var <- colMeans(se_matrix^2, na.rm = TRUE)
  mc_var <- apply(error_matrix, 2, var, na.rm = TRUE)

  # Avoid division by zero
  ok <- mc_var > 0 & is.finite(mean_model_var)
  diag_ratio <- rep(NA_real_, n_coef)
  diag_ratio[ok] <- mean_model_var[ok] / mc_var[ok]

  diag_quantiles <- quantile(
    diag_ratio,
    probs = c(0.25, 0.5, 0.75),
    na.rm = TRUE
  )

  # Diagonal-based covariance comparison for 1D terms (R1-M4: renamed to clarify
  # this compares diag(mean_se^2) vs MC covariance, not full model covariance)
  frobenius_rel <- NA_real_
  leading_eig_ratio <- NA_real_

  if (n_coef <= 100 && n_rep > n_coef + 10) {
    Sigma_MC <- cov(error_matrix, use = "pairwise.complete.obs")

    # Diagonal approximation from mean(se^2)
    # (full model Sigma_hat would require lpmatrix extraction)
    Sigma_hat_diag <- diag(mean_model_var, n_coef)

    # Frobenius relative error (diag approximation vs MC)
    frob_diff <- norm(Sigma_hat_diag - Sigma_MC, "F")
    frob_mc <- norm(Sigma_MC, "F")
    if (frob_mc > 0) frobenius_rel <- frob_diff / frob_mc

    # Leading eigenvalue ratio
    eig_mc <- eigen(Sigma_MC, symmetric = TRUE, only.values = TRUE)$values
    eig_hat <- sort(mean_model_var, decreasing = TRUE) # diagonal eigenvalues
    if (length(eig_mc) > 0 && eig_mc[1] > 0 && length(eig_hat) > 0) {
      leading_eig_ratio <- eig_hat[1] / eig_mc[1]
    }
  }

  tibble(
    term_type = term_type,
    n_coef = n_coef,
    n_rep = n_rep,
    diag_ratio_median = unname(diag_quantiles[2]),
    diag_ratio_iqr_low = unname(diag_quantiles[1]),
    diag_ratio_iqr_high = unname(diag_quantiles[3]),
    frobenius_rel_error = frobenius_rel,
    leading_eigenvalue_ratio = leading_eig_ratio,
    truth_centered = !is.null(truth_matrix)
  )
}

# DGP Settings ----------------------------------------------------------------

#' Create Study 2 DGP settings
#'
#' Full factorial in correlation and sample size. Grid levels handled in runner.
#'
#' @param n_levels Integer vector of sample sizes.
#' @returns Tibble with DGP configurations.
make_study2_settings <- function(
  n_levels = STUDY2_N_LEVELS,
  snr_levels = STUDY2_SNR_LEVELS
) {
  corr_settings <- tibble(
    corr_type = c("iid", "ar1", "fourier_pos"),
    corr_param = c(NA_real_, 0.9, 0.3)
  )

  tidyr::crossing(
    corr_settings,
    n = as.integer(n_levels),
    snr = as.numeric(snr_levels)
  ) |>
    dplyr::mutate(
      wiggliness = STUDY2_WIGGLINESS,
      error_dist = "gaussian",
      hetero_type = "none",
      hetero_param = NA_real_,
      dgp_id = dplyr::row_number()
    )
}

# Runner Functions ------------------------------------------------------------

#' Generate data on fine grid and create paired subsamples
#'
#' @param row DGP settings row.
#' @param seed RNG seed.
#' @returns List of simulations keyed by grid label.
generate_paired_data <- function(row, seed) {
  if (inherits(row, "data.frame")) row <- as.list(row)

  # Generate on finest grid
  sim_fine <- generate_benchmark_data(
    n = row$n,
    nxgrid = STUDY2_FINE_NX,
    nygrid = STUDY2_FINE_NY,
    snr = row$snr,
    error_dist = row$error_dist,
    corr_type = row$corr_type,
    corr_param = row$corr_param,
    hetero_type = row$hetero_type,
    hetero_param = row$hetero_param,
    terms = STUDY2_TERMS,
    wiggliness = row$wiggliness,
    seed = seed
  )

  # Subsample to all factorial grid combinations
  sims <- list()
  for (grid_label in names(STUDY2_GRIDS)) {
    grid_info <- STUDY2_GRIDS[[grid_label]]
    if (
      grid_info$nxgrid == STUDY2_FINE_NX &&
        grid_info$nygrid == STUDY2_FINE_NY
    ) {
      sims[[grid_label]] <- sim_fine
    } else {
      sims[[grid_label]] <- subsample_grid(
        sim_fine,
        grid_info$nxgrid,
        grid_info$nygrid
      )
    }
  }
  sims
}

#' Fit and extract metrics for one grid level
#'
#' @param sim Simulation result (at a given grid level).
#' @param alpha Significance level.
#' @param collect_vectors If TRUE, also collect coefficient vectors for
#'   covariance quality metrics.
#' @returns List with metrics tibble and optional coef_data.
fit_and_extract_grid <- function(sim, alpha = 0.10, collect_vectors = TRUE) {
  # Build formula
  frml <- build_pffr_formula(sim, sim$s_grid, k_smooth = 12, k_ff = c(12, 12))
  bs_yindex <- list(bs = "ps", k = 12, m = c(2, 1))

  # Fit model
  t0 <- Sys.time()
  fit <- pffr(
    frml,
    yind = sim$t_grid,
    data = sim$data,
    bs.yindex = bs_yindex,
    sandwich = "none"
  )
  fit_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  methods <- list(
    default = "none",
    hc = "hc",
    cluster = "cluster",
    cl2 = "cl2"
  )
  metrics_list <- list()
  coef_data_list <- list()

  for (method_name in names(methods)) {
    sw_type <- methods[[method_name]]

    # Term-level metrics
    for (term_type in c(sim$terms, "intercept")) {
      tm <- tryCatch(
        compute_term_metrics(
          fit,
          sim$truth,
          term_type,
          alpha = alpha,
          use_sandwich = sw_type,
          s_grid = sim$s_grid,
          t_grid = sim$t_grid,
          data = sim$data,
          err_struct = sim$err_struct
        ),
        error = function(e) {
          warning(sprintf(
            "compute_term_metrics failed for method=%s term=%s: %s",
            method_name,
            term_type,
            conditionMessage(e)
          ))
          NULL
        }
      )
      if (!is.null(tm)) {
        tm$method <- method_name
        tm$fit_time <- fit_time
        metrics_list[[length(metrics_list) + 1]] <- tm
      }

      # Collect coefficient vectors + truth for covariance metrics (R1-H3)
      if (collect_vectors) {
        cv <- tryCatch(
          extract_coef_and_se(
            fit,
            term_type,
            sw_type,
            truth = sim$truth,
            s_grid = sim$s_grid,
            t_grid = sim$t_grid,
            data = sim$data
          ),
          error = function(e) NULL
        )
        if (!is.null(cv)) {
          coef_data_list[[paste(method_name, term_type, sep = ".")]] <- cv
        }
      }
    }

    # E(Y) metrics
    ey_tm <- tryCatch(
      compute_ey_metrics(
        fit,
        sim$truth,
        alpha,
        use_sandwich = sw_type,
        err_struct = sim$err_struct,
        t_grid = sim$t_grid
      ),
      error = function(e) NULL
    )
    if (!is.null(ey_tm)) {
      ey_tm$method <- method_name
      ey_tm$fit_time <- fit_time
      metrics_list[[length(metrics_list) + 1]] <- ey_tm
    }
  }

  list(
    metrics = dplyr::bind_rows(metrics_list),
    fit_time = fit_time,
    coef_data = if (collect_vectors) coef_data_list else NULL
  )
}

#' Atomic save
atomic_saveRDS <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds.tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) {
    warning("atomic_saveRDS: file.rename failed for ", path)
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
}

study2_main_file_key <- function(dgp_id, n, grid_label, rep_id) {
  sprintf("dgp%03d_n%03d_grid%s_rep%03d", dgp_id, n, grid_label, rep_id)
}

#' List existing Study 2 main result files (new factorial filename pattern)
#'
#' @param output_dir Study 2 main output directory.
#' @returns Tibble with parsed file keys.
list_study2_main_files <- function(output_dir) {
  files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_n\\d+_grid\\d+x\\d+_rep\\d+\\.rds$",
    full.names = FALSE
  )
  if (length(files) == 0) {
    return(tibble(
      file = character(),
      dgp_id = integer(),
      n = integer(),
      grid_label = character(),
      rep_id = integer()
    ))
  }

  parsed <- stringr::str_match(
    files,
    "^dgp(\\d+)_n(\\d+)_grid(\\d+x\\d+)_rep(\\d+)\\.rds$"
  )

  tibble(
    file = files,
    dgp_id = as.integer(parsed[, 2]),
    n = as.integer(parsed[, 3]),
    grid_label = parsed[, 4],
    rep_id = as.integer(parsed[, 5])
  )
}

#' Build pending Study 2 main tasks for resume-friendly execution
#'
#' Determines which (dgp, rep) pairs are already complete across the requested
#' grid labels and skips them up front. Partially complete pairs are retained and
#' completed using the per-grid file checks in run_one_pair().
#'
#' @param settings DGP settings tibble from make_study2_settings().
#' @param n_rep Number of reps requested.
#' @param grid_labels Grid labels requested for this run.
#' @param output_dir Study 2 main output directory.
#' @returns List with task_grid and resume summary counts.
plan_study2_main_tasks <- function(settings, n_rep, grid_labels, output_dir) {
  grid_labels <- unique(as.character(grid_labels))

  pair_grid <- settings |>
    tidyr::crossing(rep_id = seq_len(n_rep))

  target_grid_tasks <- pair_grid |>
    dplyr::select(dgp_id, n, rep_id) |>
    tidyr::crossing(grid_label = grid_labels)

  existing_grid_tasks <- list_study2_main_files(output_dir) |>
    dplyr::filter(grid_label %in% grid_labels) |>
    dplyr::distinct(dgp_id, n, rep_id, grid_label)

  pair_status <- target_grid_tasks |>
    dplyr::count(dgp_id, n, rep_id, name = "n_grid_target") |>
    dplyr::left_join(
      existing_grid_tasks |>
        dplyr::count(dgp_id, n, rep_id, name = "n_grid_done"),
      by = c("dgp_id", "n", "rep_id")
    ) |>
    dplyr::mutate(
      n_grid_done = tidyr::replace_na(n_grid_done, 0L),
      pair_complete = n_grid_done >= n_grid_target
    )

  complete_pairs <- pair_status |>
    dplyr::filter(pair_complete) |>
    dplyr::select(dgp_id, n, rep_id)

  task_grid <- pair_grid |>
    dplyr::anti_join(complete_pairs, by = c("dgp_id", "n", "rep_id"))

  pending_grid_tasks <- target_grid_tasks |>
    dplyr::anti_join(
      existing_grid_tasks,
      by = c("dgp_id", "n", "rep_id", "grid_label")
    )

  list(
    task_grid = task_grid,
    pair_status = pair_status,
    total_pairs = nrow(pair_grid),
    complete_pairs = sum(pair_status$pair_complete),
    partial_pairs = sum(
      pair_status$n_grid_done > 0 & !pair_status$pair_complete
    ),
    pending_pairs = nrow(task_grid),
    total_grid_tasks = nrow(target_grid_tasks),
    existing_grid_tasks = nrow(existing_grid_tasks),
    pending_grid_tasks = nrow(pending_grid_tasks)
  )
}

#' Run Study 2 pilot phase
#'
#' @param n_rep Number of reps per DGP x grid.
#' @param parallel Use parallel processing?
#' @param n_workers Number of workers.
#' @param output_dir Output directory.
#' @returns Pilot summary tibble.
run_study2_pilot <- function(
  n_rep = 1L,
  parallel = TRUE,
  n_workers = max(1L, min(6L, parallel::detectCores() - 1L)),
  output_dir = "ci-benchmark/study2-grid-refinement"
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  settings <- make_study2_settings()
  grid_labels <- names(STUDY2_GRIDS)

  cat("Study 2 Pilot: Grid Refinement\n")
  cat("==============================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("n levels:", paste(sort(unique(settings$n)), collapse = ", "), "\n")
  cat("SNR levels:", paste(sort(unique(settings$snr)), collapse = ", "), "\n")
  cat("Grid levels:", paste(grid_labels, collapse = ", "), "\n")
  cat("Reps per cell:", n_rep, "\n")
  cat("Parallel:", parallel, "workers:", n_workers, "\n\n")

  # Build task grid: each task is one (DGP, rep) combination
  # (all 3 grid levels are fitted within each task for paired comparison)
  task_grid <- settings |> tidyr::crossing(rep_id = seq_len(n_rep))

  # Worker function: one (DGP, rep) → fit all grid levels
  run_pilot_task <- function(row) {
    if (inherits(row, "data.frame")) row <- as.list(row)
    seed <- STUDY2_BASE_SEED + 1000L * row$dgp_id + row$rep_id

    sims <- generate_paired_data(row, seed)

    task_results <- list()
    task_timings <- list()

    for (grid_label in grid_labels) {
      sim <- sims[[grid_label]]
      grid_info <- STUDY2_GRIDS[[grid_label]]

      res <- fit_and_extract_grid(sim, collect_vectors = FALSE)

      if (nrow(res$metrics) > 0) {
        result <- res$metrics |>
          dplyr::mutate(
            dgp_id = row$dgp_id,
            rep_id = row$rep_id,
            seed = seed,
            corr_type = row$corr_type,
            corr_param = row$corr_param,
            n = row$n,
            snr = row$snr,
            nxgrid = grid_info$nxgrid,
            nygrid = grid_info$nygrid,
            grid_label = grid_label
          )
        task_results[[grid_label]] <- result
        task_timings[[grid_label]] <- tibble(
          dgp_id = row$dgp_id,
          corr_type = row$corr_type,
          snr = row$snr,
          grid_label = grid_label,
          nxgrid = grid_info$nxgrid,
          nygrid = grid_info$nygrid,
          rep_id = row$rep_id,
          fit_time = result$fit_time[1]
        )
      } else {
        # Record failure row instead of silently dropping (R1-H2)
        task_results[[grid_label]] <- tibble(
          term_type = NA_character_,
          coverage = NA_real_,
          method = NA_character_,
          fit_time = res$fit_time,
          converged = FALSE,
          dgp_id = row$dgp_id,
          rep_id = row$rep_id,
          seed = seed,
          corr_type = row$corr_type,
          corr_param = row$corr_param,
          n = row$n,
          snr = row$snr,
          nxgrid = grid_info$nxgrid,
          nygrid = grid_info$nygrid,
          grid_label = grid_label,
          error_msg = "metrics extraction returned 0 rows"
        )
        task_timings[[grid_label]] <- tibble(
          dgp_id = row$dgp_id,
          corr_type = row$corr_type,
          snr = row$snr,
          grid_label = grid_label,
          nxgrid = grid_info$nxgrid,
          nygrid = grid_info$nygrid,
          rep_id = row$rep_id,
          fit_time = res$fit_time
        )
      }
    }

    list(
      results = dplyr::bind_rows(task_results),
      timings = dplyr::bind_rows(task_timings)
    )
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

    rows <- split(task_grid, seq_len(nrow(task_grid)))

    out_list <- furrr::future_map(
      rows,
      function(row_df) {
        if (is_dev_mode) {
          devtools::load_all(pkg_dir, quiet = TRUE)
        } else {
          library(refund)
        }
        source(normalizePath("ci-benchmark/benchmark-utils.R"), local = TRUE)
        source(normalizePath("ci-benchmark/confint-benchmark.R"), local = TRUE)
        source(
          normalizePath("ci-benchmark/sim-study-grid-refinement.R"),
          local = TRUE
        )
        tryCatch(
          run_pilot_task(row_df),
          error = function(e) {
            message(
              "Pilot error dgp=",
              row_df$dgp_id,
              " rep=",
              row_df$rep_id,
              ": ",
              e$message
            )
            list(results = tibble(), timings = tibble())
          }
        )
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )

    all_results <- lapply(out_list, `[[`, "results")
    timing_results <- lapply(out_list, `[[`, "timings")
  } else {
    all_results <- list()
    timing_results <- list()

    for (i in seq_len(nrow(task_grid))) {
      row_df <- task_grid[i, ]
      cat(sprintf(
        "\r  [%d/%d] DGP %d (corr=%s, n=%d, snr=%g), rep %d",
        i,
        nrow(task_grid),
        row_df$dgp_id,
        row_df$corr_type,
        row_df$n,
        row_df$snr,
        row_df$rep_id
      ))

      out <- tryCatch(
        run_pilot_task(row_df),
        error = function(e) {
          message("\n  Error: ", e$message)
          list(results = tibble(), timings = tibble())
        }
      )
      all_results[[i]] <- out$results
      timing_results[[i]] <- out$timings
    }
    cat("\n")
  }

  results <- dplyr::bind_rows(all_results)
  timings <- dplyr::bind_rows(timing_results)
  if (nrow(timings) == 0) {
    timings <- tibble(
      n = integer(),
      snr = numeric(),
      corr_type = character(),
      grid_label = character(),
      nxgrid = integer(),
      nygrid = integer(),
      rep_id = integer(),
      fit_time = numeric()
    )
  }

  timing_summary <- if (nrow(timings) > 0) {
    timings |>
      dplyr::group_by(n, snr, corr_type, grid_label, nxgrid, nygrid) |>
      dplyr::summarise(
        median_time = median(fit_time, na.rm = TRUE),
        mean_time = mean(fit_time, na.rm = TRUE),
        max_time = max(fit_time, na.rm = TRUE),
        n_fits = dplyr::n(),
        .groups = "drop"
      )
  } else {
    tibble(
      n = integer(),
      snr = numeric(),
      corr_type = character(),
      grid_label = character(),
      nxgrid = integer(),
      nygrid = integer(),
      median_time = numeric(),
      mean_time = numeric(),
      max_time = numeric(),
      n_fits = integer()
    )
  }

  timing_scaling <- if (nrow(timing_summary) > 0) {
    baseline_label <- "30x40"
    timing_summary |>
      dplyr::group_by(n, snr, corr_type) |>
      dplyr::mutate(
        baseline_median = median_time[grid_label == baseline_label][1],
        time_ratio_vs_baseline = median_time / baseline_median
      ) |>
      dplyr::ungroup()
  } else {
    tibble(
      n = integer(),
      snr = numeric(),
      corr_type = character(),
      grid_label = character(),
      nxgrid = integer(),
      nygrid = integer(),
      median_time = numeric(),
      mean_time = numeric(),
      max_time = numeric(),
      n_fits = integer(),
      baseline_median = numeric(),
      time_ratio_vs_baseline = numeric()
    )
  }

  cat("\nPilot Timing Summary:\n")
  print(timing_summary)
  cat(
    "\nPilot Timing Scaling (median_time / 30x40 baseline within n,snr,corr):\n"
  )
  print(
    timing_scaling |>
      dplyr::select(
        n,
        snr,
        corr_type,
        grid_label,
        median_time,
        time_ratio_vs_baseline
      ) |>
      dplyr::arrange(n, snr, corr_type, grid_label)
  )

  # Coverage summary by grid
  if (nrow(results) > 0) {
    cov_summary <- results |>
      dplyr::filter(!is.na(coverage)) |>
      dplyr::group_by(n, snr, grid_label, corr_type, method, term_type) |>
      dplyr::summarise(
        mean_coverage = mean(coverage, na.rm = TRUE),
        mean_z_sd = mean(z_sd, na.rm = TRUE),
        n_rep = dplyr::n(),
        .groups = "drop"
      )

    cat("\nPilot Coverage Summary:\n")
    print(
      cov_summary |>
        dplyr::filter(term_type %in% c("ff", "linear")) |>
        dplyr::arrange(n, snr, corr_type, grid_label, method, term_type),
      n = 100
    )
  }

  pilot_output <- list(
    results = results,
    timing_summary = timing_summary,
    timing_scaling = timing_scaling
  )
  atomic_saveRDS(pilot_output, file.path(output_dir, "pilot_summary.rds"))

  pilot_output
}

#' Run Study 2 main phase
#'
#' @param n_rep Number of reps per DGP x grid cell.
#' @param grid_labels Character vector of grid labels to run. Default uses all
#'   9 factorial nxgrid x nygrid combinations (e.g. "30x40", "60x120").
#' @param parallel Use parallel processing?
#' @param n_workers Number of workers.
#' @param output_dir Output directory.
#' @param alpha Significance level.
#' @returns Combined results tibble.
run_study2_main <- function(
  n_rep = 50L,
  grid_labels = names(STUDY2_GRIDS),
  parallel = TRUE,
  n_workers = max(1L, min(6L, parallel::detectCores() - 1L)),
  output_dir = "ci-benchmark/study2-grid-refinement/main",
  alpha = 0.10
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  settings <- make_study2_settings()

  # Validate grid labels
  valid_labels <- names(STUDY2_GRIDS)
  bad <- setdiff(grid_labels, valid_labels)
  if (length(bad) > 0) {
    warning(
      "Unknown grid labels: ",
      paste(bad, collapse = ", "),
      ". Using all grids instead."
    )
    grid_labels <- valid_labels
  }

  cat("Study 2 Main: Grid Refinement (factorial n x snr x nxgrid x nygrid)\n")
  cat("===================================================================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("n levels:", paste(sort(unique(settings$n)), collapse = ", "), "\n")
  cat("SNR levels:", paste(sort(unique(settings$snr)), collapse = ", "), "\n")
  cat("Grid levels:", paste(grid_labels, collapse = " vs "), "\n")
  cat("Reps per cell:", n_rep, "\n")
  cat(
    "Total fits (planned):",
    nrow(settings) * length(grid_labels) * n_rep,
    "\n"
  )
  cat("Output dir:", output_dir, "\n\n")

  resume_plan <- plan_study2_main_tasks(
    settings = settings,
    n_rep = n_rep,
    grid_labels = grid_labels,
    output_dir = output_dir
  )
  task_grid <- resume_plan$task_grid

  cat("Resume scan:\n")
  cat("  Pair tasks total:", resume_plan$total_pairs, "\n")
  cat("  Pair tasks complete (skipped):", resume_plan$complete_pairs, "\n")
  cat("  Pair tasks partial (resume):", resume_plan$partial_pairs, "\n")
  cat("  Pair tasks pending:", resume_plan$pending_pairs, "\n")
  cat("  Grid tasks existing:", resume_plan$existing_grid_tasks, "\n")
  cat("  Grid tasks pending:", resume_plan$pending_grid_tasks, "\n\n")

  # Worker function for one (DGP, rep) pair across all grid levels
  run_one_pair <- function(row, rep_id, alpha) {
    seed <- STUDY2_BASE_SEED + 1000L * row$dgp_id + rep_id

    # Generate paired data (same seed for all grids)
    sims <- generate_paired_data(row, seed)

    pair_results <- list()
    for (grid_label in grid_labels) {
      # Check if already done
      save_key <- study2_main_file_key(
        dgp_id = row$dgp_id,
        n = row$n,
        grid_label = grid_label,
        rep_id = rep_id
      )
      save_path <- file.path(output_dir, paste0(save_key, ".rds"))

      if (file.exists(save_path)) {
        obj <- readRDS(save_path)
        # Saved objects are list(metrics=..., coef_data=...) — extract metrics
        pair_results[[grid_label]] <- if (
          is.list(obj) && "metrics" %in% names(obj)
        ) {
          obj$metrics
        } else {
          obj
        }
        next
      }

      sim <- sims[[grid_label]]
      grid_info <- STUDY2_GRIDS[[grid_label]]

      res <- fit_and_extract_grid(sim, alpha = alpha, collect_vectors = TRUE)

      if (nrow(res$metrics) > 0) {
        result <- res$metrics |>
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
            grid_label = grid_label
          )

        # Save per-rep result (with coef vectors for covariance quality)
        save_obj <- list(metrics = result, coef_data = res$coef_data)
        atomic_saveRDS(save_obj, save_path)
        pair_results[[grid_label]] <- result
      } else {
        # Record failure row instead of silently dropping (R1-H2)
        result <- tibble(
          term_type = NA_character_,
          coverage = NA_real_,
          method = NA_character_,
          fit_time = res$fit_time,
          converged = FALSE,
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
          error_msg = "metrics extraction returned 0 rows"
        )
        save_obj <- list(metrics = result, coef_data = NULL)
        atomic_saveRDS(save_obj, save_path)
        pair_results[[grid_label]] <- result
      }
    }

    dplyr::bind_rows(pair_results)
  }

  # Run all (DGP, rep) pairs
  if (nrow(task_grid) == 0) {
    cat("No pending Study 2 pair tasks. Loading existing results only.\n")
    results <- tibble()
  } else if (parallel && requireNamespace("furrr", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multicore, workers = n_workers)

    is_dev_mode <- file.exists("DESCRIPTION")
    pkg_dir <- normalizePath(".")

    rows <- split(task_grid, seq_len(nrow(task_grid)))

    results <- furrr::future_map_dfr(
      rows,
      function(row_df) {
        row <- as.list(row_df)
        if (is_dev_mode) {
          devtools::load_all(pkg_dir, quiet = TRUE)
        } else {
          library(refund)
        }
        source(normalizePath("ci-benchmark/benchmark-utils.R"), local = TRUE)
        source(normalizePath("ci-benchmark/confint-benchmark.R"), local = TRUE)
        source(
          normalizePath("ci-benchmark/sim-study-grid-refinement.R"),
          local = TRUE
        )

        tryCatch(
          run_one_pair(row, row$rep_id, alpha = alpha),
          error = function(e) {
            message(
              "Error dgp=",
              row$dgp_id,
              " rep=",
              row$rep_id,
              ": ",
              e$message
            )
            NULL
          }
        )
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  } else {
    # Sequential
    results <- list()

    for (i in seq_len(nrow(task_grid))) {
      row <- as.list(task_grid[i, ])
      cat(sprintf(
        "\r[%d/%d] dgp=%d (corr=%s, n=%d, snr=%g), rep=%d",
        i,
        nrow(task_grid),
        row$dgp_id,
        row$corr_type,
        row$n,
        row$snr,
        row$rep_id
      ))

      res <- tryCatch(
        run_one_pair(row, row$rep_id, alpha = alpha),
        error = function(e) {
          message("\nError: ", e$message)
          NULL
        }
      )
      if (!is.null(res)) results[[i]] <- res
    }
    cat("\n")
    results <- dplyr::bind_rows(results)
  }

  # Load all results
  all_results <- load_study2_results(output_dir)

  # Compute covariance quality metrics
  cov_quality <- compute_study2_cov_quality(output_dir, settings, grid_labels)

  # Save combined
  atomic_saveRDS(
    all_results,
    file.path(dirname(output_dir), "main_results_combined.rds")
  )
  atomic_saveRDS(cov_quality, file.path(dirname(output_dir), "cov_quality.rds"))

  cat("Total results:", nrow(all_results), "rows\n")
  list(results = all_results, cov_quality = cov_quality)
}

#' Load Study 2 main results from incremental saves
#'
#' @param output_dir Directory with per-rep RDS files.
#' @returns Combined metrics tibble.
load_study2_results <- function(output_dir) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_n\\d+_grid\\w+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(all_files) == 0) return(tibble())

  results_list <- lapply(all_files, function(f) {
    obj <- readRDS(f)
    if (is.list(obj) && "metrics" %in% names(obj)) obj$metrics else obj
  })
  dplyr::bind_rows(results_list)
}

#' Compute Study 2 covariance quality metrics across reps
#'
#' @param output_dir Main results directory.
#' @param settings DGP settings tibble.
#' @param grid_labels Grid labels to process.
#' @returns Tibble with covariance quality metrics.
compute_study2_cov_quality <- function(output_dir, settings, grid_labels) {
  cov_results <- list()

  for (s in seq_len(nrow(settings))) {
    row <- as.list(settings[s, ])

    for (grid_label in grid_labels) {
      for (method_name in c("default", "hc", "cluster", "cl2")) {
        for (term_type in c(
          "linear",
          "smooth",
          "concurrent",
          "ff",
          "intercept"
        )) {
          key <- paste(method_name, term_type, sep = ".")

          # Collect coefficient vectors and truth across reps (R1-H3)
          est_list <- list()
          se_list <- list()
          truth_list <- list()

          files <- list.files(
            output_dir,
            pattern = sprintf(
              "^dgp%03d_n%03d_grid%s_rep\\d+\\.rds$",
              row$dgp_id,
              row$n,
              grid_label
            ),
            full.names = TRUE
          )

          for (f in files) {
            obj <- tryCatch(readRDS(f), error = function(e) NULL)
            if (is.null(obj) || is.null(obj$coef_data)) next

            cv <- obj$coef_data[[key]]
            if (is.null(cv)) next

            est_list[[length(est_list) + 1]] <- cv$est
            se_list[[length(se_list) + 1]] <- cv$se
            # Collect truth for centering (R1-H3)
            if (!is.null(cv$truth_vals)) {
              truth_list[[length(truth_list) + 1]] <- cv$truth_vals
            }
          }

          if (length(est_list) < 10) next

          # Align vectors (should all be same length)
          n_coef <- length(est_list[[1]])
          ok <- vapply(est_list, function(x) length(x) == n_coef, logical(1))
          est_list <- est_list[ok]
          se_list <- se_list[ok]

          # Also align truth (and check SE alignment per R1-L3)
          truth_mat <- NULL
          if (length(truth_list) == length(est_list)) {
            ok_truth <- vapply(
              truth_list,
              function(x) length(x) == n_coef,
              logical(1)
            )
            if (all(ok_truth)) {
              truth_mat <- do.call(rbind, truth_list)
            }
          }

          if (length(est_list) < 10) next

          est_mat <- do.call(rbind, est_list)
          se_mat <- do.call(rbind, se_list)

          cov_metrics <- compute_cov_quality_metrics(
            est_mat,
            se_mat,
            term_type,
            truth_matrix = truth_mat
          )

          cov_metrics$dgp_id <- row$dgp_id
          cov_metrics$corr_type <- row$corr_type
          cov_metrics$n <- row$n
          cov_metrics$snr <- row$snr
          cov_metrics$grid_label <- grid_label
          cov_metrics$method <- method_name

          cov_results[[length(cov_results) + 1]] <- cov_metrics
        }
      }
    }
  }

  dplyr::bind_rows(cov_results)
}

# Summarization ---------------------------------------------------------------

#' Summarize Study 2 results
#'
#' @param results Results tibble.
#' @returns Summary tibble.
summarize_study2 <- function(results) {
  results |>
    dplyr::filter(!is.na(coverage)) |>
    dplyr::group_by(
      n,
      snr,
      corr_type,
      grid_label,
      nxgrid,
      nygrid,
      method,
      term_type
    ) |>
    dplyr::summarise(
      mean_coverage = mean(coverage, na.rm = TRUE),
      se_coverage = sd(coverage, na.rm = TRUE) / sqrt(dplyr::n()),
      mean_width = mean(mean_width, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_z_sd = mean(z_sd, na.rm = TRUE),
      mean_fit_time = mean(fit_time, na.rm = TRUE),
      n_reps = dplyr::n(),
      .groups = "drop"
    )
}

# Validation Checks -----------------------------------------------------------

#' Validate paired grid subsampling
#'
#' Verify same seeds + subsample produce identical truth functions.
#'
#' @param seed RNG seed.
#' @returns List with pass and max_difference.
validate_paired_grids <- function(seed = 5555) {
  settings <- make_study2_settings()
  row <- as.list(settings[1, ]) # IID case

  sims <- generate_paired_data(row, seed)

  # Check: truth on coarse grid should match subsampled fine truth
  fine_truth_ff <- sims$fine$truth$beta$`ff(X1)`
  coarse_truth_ff <- sims$coarse$truth$beta$`ff(X1)`

  # Subsample fine truth to coarse indices
  s_idx <- sims$coarse$s_idx
  t_idx <- sims$coarse$t_idx
  fine_subsampled <- fine_truth_ff[s_idx, t_idx]

  max_diff <- max(abs(coarse_truth_ff - fine_subsampled))

  # Linear term
  fine_truth_lin <- sims$fine$truth$beta$zlin
  coarse_truth_lin <- sims$coarse$truth$beta$zlin
  fine_sub_lin <- fine_truth_lin[t_idx]
  max_diff_lin <- max(abs(coarse_truth_lin - fine_sub_lin))

  max_diff_total <- max(max_diff, max_diff_lin)

  list(
    pass = max_diff_total < 1e-10,
    max_difference_ff = max_diff,
    max_difference_lin = max_diff_lin,
    max_difference = max_diff_total
  )
}

# Main Entry Point -------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args) > 0) args[1] else "smoke"

  base_dir <- "ci-benchmark/study2-grid-refinement"

  if (mode == "smoke") {
    cat("Study 2: Smoke test (2 reps)\n\n")

    # Validation checks
    cat("Paired grid validation: ")
    v <- validate_paired_grids()
    cat(
      if (v$pass) "PASS" else "FAIL",
      sprintf("(max diff = %.2e)\n", v$max_difference)
    )

    pilot <- run_study2_pilot(n_rep = 2L, output_dir = base_dir)
  } else if (mode == "pilot") {
    cat("Study 2: Pilot phase (30 reps)\n\n")
    pilot <- run_study2_pilot(n_rep = 30L, output_dir = base_dir)
  } else if (mode == "main") {
    cat("Using full factorial design with all nxgrid x nygrid combinations\n")

    main_results <- run_study2_main(
      n_rep = 50L,
      grid_labels = names(STUDY2_GRIDS), # all 9 nxgrid x nygrid combos
      parallel = TRUE,
      n_workers = max(1L, min(10L, parallel::detectCores() - 1L)),
      output_dir = file.path(base_dir, "main")
    )

    # Summary
    if (nrow(main_results$results) > 0) {
      summary_df <- summarize_study2(main_results$results)

      cat("\n========== STUDY 2 COVERAGE SUMMARY ==========\n")
      print(
        summary_df |>
          dplyr::filter(term_type %in% c("ff", "linear", "E(Y)")) |>
          dplyr::arrange(n, snr, corr_type, grid_label, method, term_type),
        n = 200
      )

      cat("\n========== COVARIANCE QUALITY ==========\n")
      print(
        main_results$cov_quality |>
          dplyr::arrange(n, snr, corr_type, grid_label, method, term_type),
        n = 100
      )

      write.csv(
        summary_df,
        file.path(base_dir, "main_summary.csv"),
        row.names = FALSE
      )
    }
  } else {
    n_rep <- as.integer(mode)
    if (is.na(n_rep)) n_rep <- 2L
    pilot <- run_study2_pilot(n_rep = n_rep, output_dir = base_dir)
  }

  cat("\nStudy 2 complete.\n")
}
