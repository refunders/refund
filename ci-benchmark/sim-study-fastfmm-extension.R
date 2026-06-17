# ===========================================================================
# fastFMM (FUI) Extension: Add E4 competitor to Study 2 (and optionally Study 1)
# ===========================================================================
#
# Re-runs the same seeds as the production study (Study 2 base 4001; Study 1
# base 3001 via command-line arg), regenerates identical data, fits
# fastFMM::fui() on the terms it can estimate (intercept + linear scalar;
# concurrent if present), and extracts:
#
#   method="fastfmm"     — pointwise CI coverage + mean width
#   method="fastfmm_sim" — simultaneous/joint band coverage + band width
#
# SCOPE LIMITATIONS (fastFMM cannot fit):
#   - ff (function-on-function): ∫X(s)β(s,t)ds is not supported by FUI.
#   - s(zsmoo): fastFMM treats all fixed effects linearly (step 1 = pointwise
#     LME with fixed scalar predictors). The nonlinear smooth effect in pffr
#     maps to a linear term in fastFMM — not comparable.
#   - Concurrent (Xconc) is ONLY included if concurrent=TRUE is passed to fui()
#     AND the data was generated with that term. For our Study 2 DGP,
#     concurrent=TRUE is used since Xconc is evaluated on the response grid.
#
# Terms compared (Study 2 Gaussian):
#   intercept, linear (zlin). NOT ff, NOT s(zsmoo) (linear zsmoo included
#   in the fui model but not compared — truth is nonlinear).
#   If Xconc present: concurrent is also compared.
#
# Study 1 (non-Gaussian: Poisson/Binomial):
#   intercept, linear (zlin). analytic=FALSE, boot=500. IMPLEMENT but
#   do NOT run locally (LRZ cluster only).
#
# Output dirs:
#   ci-benchmark/study2-fastfmm/   (Study 2, one rds per dgp/grid/rep)
#   ci-benchmark/study1-fastfmm/   (Study 1, one rds per dgp/rep)
#
# File naming mirrors study2-competitors: dgp%03d_n%03d_y%03d_rep%03d.rds
#
# Usage:
#   Rscript ci-benchmark/sim-study-fastfmm-extension.R [mode] [study]
#   mode:  "smoke" (1 rep, 1 dgp), "pilot" (10 reps), "full" (50 reps)
#   study: "study2" (default, Gaussian) or "study1" (non-Gaussian, LRZ only)
# ===========================================================================

# Setup -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(fastFMM)
})

if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(refund)
}

source("ci-benchmark/benchmark-utils.R")
source("ci-benchmark/confint-benchmark.R")

# Source Study 2 for DGP functions, settings, constants. The guarded main
# block (sys.nframe() == 0) prevents execution when sourced.
source("ci-benchmark/sim-study-grid-refinement.R")

# Also source Study 1 DGP for the non-Gaussian path.
source("ci-benchmark/sim-study-nongaussian-sandwich.R")

# Constants -------------------------------------------------------------------

FASTFMM_ALPHA <- 0.10 # 90% CIs (nominal 1-alpha)

# Study 2 output dir
FASTFMM_STUDY2_DIR <- "ci-benchmark/study2-fastfmm"

# Study 1 output dir
FASTFMM_STUDY1_DIR <- "ci-benchmark/study1-fastfmm"

# Terms that fastFMM can estimate and that we compare.
# "linear" = zlin (linear scalar effect, apples-to-apples).
# "intercept" = intercept function.
# NOTE: "concurrent" included for Study 2 (Xconc) when present in DGP.
# NOTE: "zsmoo" is fitted by fui linearly but NOT compared — truth is
#   nonlinear s(zsmoo) so it is not apples-to-apples.
FASTFMM_COMPARABLE_TERMS_STUDY2 <- c("intercept", "linear", "concurrent")
FASTFMM_COMPARABLE_TERMS_STUDY1 <- c("intercept", "linear")

# Atomic save -----------------------------------------------------------------

atomic_saveRDS_ff <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds.tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) {
    warning("atomic_saveRDS_ff: file.rename failed for ", path)
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
}

# FUI data preparation --------------------------------------------------------

#' Build a fastFMM-compatible data frame from a simulation result
#'
#' FUI requires:
#'   - Response: an N x L matrix stored as a single column named with the
#'     LHS of the formula. We store it as a regular matrix (not I(...)).
#'   - Scalar covariates: standard columns.
#'   - Subject ID: a unique integer column "id" (one curve per unit).
#'
#' NOTE: In our DGP, each subject contributes exactly one functional
#' observation (no repeated measures). This makes the random intercept
#' (1|id) degenerate. We use override_zero_var=TRUE in fui() to proceed.
#' The re variance estimate will be ~0, which is correct: no genuine
#' between-subject random variation beyond the fixed effects.
#'
#' @param sim Simulation result from generate_benchmark_data() or
#'   subsample_grid().
#' @returns A data frame suitable for fui(), with response matrix column
#'   named "Y", scalar cols "zlin", "zsmoo", "id", and "Xconc" if present.
build_fastfmm_data <- function(sim) {
  dat <- sim$data
  n <- nrow(dat)

  Y_mat <- as.matrix(dat$Y) # n x nygrid

  fmm_df <- data.frame(
    id = seq_len(n),
    zlin = as.numeric(dat$zlin)
  )

  # zsmoo only present in Study 2 (pffr smooth term). Study 1 has no zsmoo.
  # fastFMM cannot estimate the nonlinear smooth anyway (linear approximation),
  # so we include it for completeness in Study 2 but skip in Study 1.
  has_zsmoo <- !is.null(dat$zsmoo) && length(dat$zsmoo) == n
  if (has_zsmoo) {
    fmm_df$zsmoo <- as.numeric(dat$zsmoo)
  }

  # Add concurrent covariate if present: Xconc is n x nygrid in sim$data
  has_concurrent <- "concurrent" %in%
    sim$terms &&
    !is.null(dat$Xconc) &&
    is.matrix(as.matrix(dat$Xconc))

  if (has_concurrent) {
    fmm_df$Xconc <- I(as.matrix(dat$Xconc))
  }

  # Attach Y as a plain matrix column (fui detects by column-name pattern)
  fmm_df$Y <- I(Y_mat)

  fmm_df
}

# FUI fitting -----------------------------------------------------------------

#' Fit fui() on a simulation result and return the model object
#'
#' For Gaussian family (Study 2): analytic=TRUE, parallel=FALSE.
#' For non-Gaussian (Study 1): analytic=FALSE, n_boots=500, boot_type="cluster".
#'
#' Handles the degenerate-RE case (one curve per subject) with
#' override_zero_var=TRUE.
#'
#' @param sim Simulation result.
#' @param family Character family: "gaussian", "poisson", or "binomial".
#' @param n_boots Number of bootstrap samples for non-Gaussian (default 500).
#' @param fui_seed Seed for FUI bootstrap (passed to fui seed argument).
#' @returns A fui model object or NULL on failure.
fit_fui <- function(sim, family = "gaussian", n_boots = 500L, fui_seed = 1L) {
  fmm_df <- build_fastfmm_data(sim)
  has_concurrent <- !is.null(fmm_df$Xconc)
  has_zsmoo <- !is.null(fmm_df$zsmoo)

  # Build formula:
  # - zlin: linear scalar (apples-to-apples with pffr's zlin term)
  # - zsmoo: included as linear ONLY when present in the data (Study 2 only).
  #   NOT compared — truth is nonlinear; fui estimates it linearly (not comparable).
  # - Xconc: concurrent functional covariate if present in data
  # - (1|id): random intercept (degenerate but required by FUI for variance)
  rhs_terms <- c(
    "zlin",
    if (has_zsmoo) "zsmoo",
    if (has_concurrent) "Xconc",
    "(1 | id)"
  )
  frml <- as.formula(paste("Y ~", paste(rhs_terms, collapse = " + ")))

  is_gaussian <- (family == "gaussian")

  fit <- tryCatch(
    suppressWarnings(
      fui(
        formula = frml,
        data = fmm_df,
        family = family,
        analytic = is_gaussian,
        var = TRUE,
        parallel = FALSE,
        silent = TRUE,
        n_boots = if (!is_gaussian) n_boots else 500L,
        boot_type = if (!is_gaussian) "cluster" else NULL,
        seed = fui_seed,
        # concurrent=TRUE tells fui that Xconc is a functional covariate
        # evaluated at each response grid point (same as the Y domain).
        # This enables the pointwise concurrent X(t)beta(t) model.
        concurrent = has_concurrent,
        override_zero_var = TRUE
      )
    ),
    error = function(e) {
      warning(sprintf("fui() failed: %s", conditionMessage(e)))
      NULL
    }
  )
  fit
}

# Metric extraction -----------------------------------------------------------

#' Compute joint CI critical value at given alpha from betaHat_var
#'
#' var_analytic() hardcodes qn as the 0.95 quantile of max|Z(t)|, giving 95%
#' simultaneous bands. We recompute at any alpha using the same Gaussian
#' multiplier approach (same formula, different quantile).
#'
#' @param Sigma L x L covariance matrix for one coefficient function.
#' @param alpha Significance level (default 0.10 for 90% band).
#' @param N Number of Monte Carlo samples (default 10000 as in fui).
#' @param seed Monte Carlo seed.
#' @returns Scalar critical value qn at level 1-alpha.
compute_qn_at_alpha <- function(Sigma, alpha = 0.10, N = 10000L, seed = 1L) {
  L <- nrow(Sigma)
  sqrt_diag <- sqrt(diag(Sigma))
  # Avoid division by zero
  safe_diag <- pmax(sqrt_diag, 1e-14)
  # Standardize correlation matrix
  S_scl <- diag(1 / safe_diag, L)
  Sigma_corr <- as.matrix(S_scl %*% Sigma %*% S_scl)
  # Ensure positive definite
  Sigma_corr <- 0.5 * (Sigma_corr + t(Sigma_corr))
  diag(Sigma_corr) <- 1
  set.seed(seed)
  x_sample <- tryCatch(
    abs(mvtnorm::rmvnorm(N, rep(0, L), Sigma_corr)),
    error = function(e) NULL
  )
  if (is.null(x_sample)) return(NA_real_)
  un <- apply(x_sample, 1, max)
  stats::quantile(un, 1 - alpha)
}

#' Extract per-term coverage and width metrics from a fui fit
#'
#' Maps fui coefficient curves to truth on sim$t_grid and computes:
#'   pointwise coverage = fraction of t-grid points where CI contains truth
#'   joint coverage     = indicator that ALL t-grid points are inside the
#'                        simultaneous band
#'   mean_width (pointwise) and mean_width_joint (simultaneous band)
#'
#' Term mapping (row names of fui$betaHat):
#'   "(Intercept)"  -> term "intercept" -> truth$beta$intercept
#'   "zlin"         -> term "linear"    -> truth$beta$zlin
#'   "Xconc"        -> term "concurrent"-> truth$beta$Xconc  (if concurrent=TRUE)
#'   "zsmoo"        -> EXCLUDED (truth nonlinear; not comparable)
#'
#' FUI argvals: internally 1:L (indices). We map index k to t_grid[k] via
#' approx(). The L entries in betaHat correspond to the L columns of Y (the
#' functional response), i.e., betaHat[r, k] estimates beta_r(t_grid[k]).
#'
#' Alpha: 90% CIs (alpha=0.10).
#'   Pointwise: beta +/- z_{0.95} * sqrt(diag(betaHat_var[,,r]))
#'   Simultaneous: note fui's stored qn always targets 95% (hardcoded).
#'     We recompute a 90% joint qn via compute_qn_at_alpha() using the
#'     same Gaussian multiplier method as fui's var_analytic.R.
#'
#' @param fit fui model object.
#' @param sim Simulation result (for truth and t_grid).
#' @param alpha Significance level (default 0.10).
#' @param comparable_terms Character vector of term types to extract.
#' @param qn_seed Seed for Monte Carlo critical value simulation.
#' @returns Tibble with columns: term_type, coverage, mean_width,
#'   coverage_joint, mean_width_joint. One row per comparable term found.
#'   Returns empty tibble if fit is NULL or no terms found.
extract_fui_metrics <- function(
  fit,
  sim,
  alpha = FASTFMM_ALPHA,
  comparable_terms = FASTFMM_COMPARABLE_TERMS_STUDY2,
  qn_seed = 42L
) {
  if (is.null(fit) || is.null(fit$betaHat)) return(tibble())

  t_grid <- sim$t_grid
  L <- length(t_grid)
  truth <- sim$truth

  # Row names of betaHat identify terms
  term_rows <- rownames(fit$betaHat)
  if (is.null(term_rows)) return(tibble())

  # argvals is 1:L inside fui (analytic path). Map index -> t_grid value.
  # For the analytic path fui$argvals is always 1:L regardless of t_grid.
  fui_argvals <- fit$argvals # 1:L
  if (length(fui_argvals) != L) {
    warning(sprintf(
      "extract_fui_metrics: fui argvals length (%d) != t_grid length (%d); skipping",
      length(fui_argvals),
      L
    ))
    return(tibble())
  }

  # Mapping from fui row name -> our term type + truth extractor
  # (zsmoo is excluded: truth is nonlinear but fui models it linearly)
  term_map <- list(
    "(Intercept)" = list(
      type = "intercept",
      truth = function(b) {
        tv <- b$intercept
        if (is.null(tv)) return(NULL)
        # interpolate truth onto t_grid (should already be on t_grid)
        stats::approx(t_grid, tv, xout = t_grid, rule = 2)$y
      }
    ),
    "zlin" = list(
      type = "linear",
      truth = function(b) {
        tv <- b$zlin
        if (is.null(tv)) return(NULL)
        stats::approx(t_grid, tv, xout = t_grid, rule = 2)$y
      }
    ),
    "Xconc" = list(
      type = "concurrent",
      truth = function(b) {
        tv <- b$Xconc
        if (is.null(tv)) return(NULL)
        stats::approx(t_grid, tv, xout = t_grid, rule = 2)$y
      }
    )
  )

  crit_pw <- qnorm(1 - alpha / 2) # z_{0.95} ≈ 1.645 for 90% CI

  # For concurrent models, fui names the Xconc coefficient "Xconc_1" (using
  # the first argval index). Use prefix matching for the Xconc term.
  find_row <- function(row_nm, term_rows) {
    # Exact match first
    r <- which(term_rows == row_nm)
    if (length(r) > 0) return(r[1])
    # Prefix match (handles Xconc -> Xconc_1, Xconc_2, etc.)
    r <- which(startsWith(term_rows, paste0(row_nm, "_")))
    if (length(r) > 0) return(r[1])
    NULL
  }

  results <- list()
  for (row_nm in names(term_map)) {
    info <- term_map[[row_nm]]
    if (!(info$type %in% comparable_terms)) next
    r <- find_row(row_nm, term_rows)
    if (is.null(r)) next # term not in model

    # Point estimate (1 x L row of betaHat)
    est <- as.numeric(fit$betaHat[r, ])

    # Truth on t_grid
    truth_vals <- tryCatch(info$truth(truth$beta), error = function(e) NULL)
    if (is.null(truth_vals) || length(truth_vals) != L) next

    # ---- Pointwise CI -------------------------------------------------------
    # Requires betaHat_var: L x L x p array; diagonal gives pointwise SE^2
    if (!is.null(fit$betaHat_var)) {
      var_diag <- diag(fit$betaHat_var[,, r]) # length L
      se_pw <- sqrt(pmax(var_diag, 0))
      lower_pw <- est - crit_pw * se_pw
      upper_pw <- est + crit_pw * se_pw
      covered_pw <- (truth_vals >= lower_pw) & (truth_vals <= upper_pw)
      cov_pw <- mean(covered_pw, na.rm = TRUE)
      width_pw <- mean(upper_pw - lower_pw, na.rm = TRUE)
    } else {
      cov_pw <- NA_real_
      width_pw <- NA_real_
    }

    # ---- Simultaneous/joint band -------------------------------------------
    # fui's stored qn is always the 0.95 quantile (hardcoded in var_analytic.R
    # line 342). To compare at our nominal 1-alpha=0.90, we recompute qn at
    # the 0.90 quantile using the same Gaussian multiplier approach.
    if (!is.null(fit$betaHat_var)) {
      Sigma_r <- fit$betaHat_var[,, r]
      qn_90 <- tryCatch(
        compute_qn_at_alpha(Sigma_r, alpha = alpha, N = 10000L, seed = qn_seed),
        error = function(e) NA_real_
      )
      if (is.finite(qn_90)) {
        se_pw_j <- sqrt(pmax(diag(Sigma_r), 0))
        lower_j <- est - qn_90 * se_pw_j
        upper_j <- est + qn_90 * se_pw_j
        covered_j <- (truth_vals >= lower_j) & (truth_vals <= upper_j)
        # Joint coverage = 1 if ALL points covered
        cov_joint <- as.numeric(all(covered_j, na.rm = TRUE))
        width_joint <- mean(upper_j - lower_j, na.rm = TRUE)
      } else {
        cov_joint <- NA_real_
        width_joint <- NA_real_
      }
    } else {
      cov_joint <- NA_real_
      width_joint <- NA_real_
    }

    results[[length(results) + 1]] <- tibble(
      term_type = info$type,
      coverage = cov_pw, # pointwise coverage
      mean_width = width_pw, # pointwise mean width
      coverage_joint = cov_joint, # simultaneous joint coverage (0/1)
      mean_width_joint = width_joint # simultaneous band mean width
    )
  }

  dplyr::bind_rows(results)
}

# Null-row helpers ------------------------------------------------------------

null_row_study2 <- function(
  row,
  rep_id,
  seed,
  grid_label,
  grid_info,
  error_msg = NA_character_
) {
  term_types <- FASTFMM_COMPARABLE_TERMS_STUDY2
  tibble(
    term_type = term_types,
    coverage = NA_real_,
    mean_width = NA_real_,
    coverage_joint = NA_real_,
    mean_width_joint = NA_real_,
    method = NA_character_,
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
    fit_time = NA_real_,
    converged = FALSE,
    error_msg = error_msg
  )
}

null_row_study1 <- function(row, rep_id, seed, error_msg = NA_character_) {
  term_types <- FASTFMM_COMPARABLE_TERMS_STUDY1
  tibble(
    term_type = term_types,
    coverage = NA_real_,
    mean_width = NA_real_,
    coverage_joint = NA_real_,
    mean_width_joint = NA_real_,
    method = NA_character_,
    dgp_id = row$dgp_id,
    rep_id = rep_id,
    seed = seed,
    family = row$family,
    n = row$n,
    nxgrid = row$nxgrid,
    nygrid = row$nygrid,
    corr_type = row$corr_type,
    corr_param = row$corr_param,
    fit_time = NA_real_,
    converged = FALSE,
    error_msg = error_msg
  )
}

# Study 2 runner --------------------------------------------------------------

#' Run one (DGP, rep) pair across all grid levels — fastFMM extension
#'
#' Same seed and data as Study 2 production run (STUDY2_BASE_SEED + ...).
#'
#' @param row DGP settings list.
#' @param rep_id Integer replicate id.
#' @param grid_labels Character vector of grid labels to process.
#' @param output_dir Character path to output directory.
#' @param alpha Significance level.
#' @returns Tibble combining fastFMM metrics for all grid levels.
run_one_pair_fastfmm2 <- function(
  row,
  rep_id,
  grid_labels,
  output_dir,
  alpha = FASTFMM_ALPHA
) {
  if (inherits(row, "data.frame")) row <- as.list(row)
  seed <- STUDY2_BASE_SEED + 1000L * row$dgp_id + rep_id

  # Generate paired data — identical to Study 2 production
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
      null_row_study2(
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
    nygrid <- study2_parse_nygrid(grid_label)
    file_key <- sprintf(
      "dgp%03d_n%03d_y%03d_rep%03d",
      row$dgp_id,
      row$n,
      nygrid,
      rep_id
    )
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

    # Determine comparable terms (concurrent only if Xconc present in DGP)
    has_concurrent <- "concurrent" %in% sim$terms && !is.null(sim$data$Xconc)
    comp_terms <- c(
      "intercept",
      "linear",
      if (has_concurrent) "concurrent" else NULL
    )

    # Fit fui
    t0 <- Sys.time()
    fit <- fit_fui(sim, family = "gaussian")
    fit_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (is.null(fit)) {
      result <- null_row_study2(
        row,
        rep_id,
        seed,
        grid_label,
        grid_info,
        "fui() returned NULL"
      )
      result$fit_time <- fit_time
      atomic_saveRDS_ff(result, save_path)
      pair_results[[grid_label]] <- result
      next
    }

    # Extract metrics
    metrics <- tryCatch(
      extract_fui_metrics(
        fit,
        sim,
        alpha = alpha,
        comparable_terms = comp_terms
      ),
      error = function(e) {
        warning(sprintf(
          "extract_fui_metrics failed dgp=%d rep=%d grid=%s: %s",
          row$dgp_id,
          rep_id,
          grid_label,
          conditionMessage(e)
        ))
        tibble()
      }
    )

    if (nrow(metrics) == 0) {
      result <- null_row_study2(
        row,
        rep_id,
        seed,
        grid_label,
        grid_info,
        "extract_fui_metrics returned 0 rows"
      )
      result$fit_time <- fit_time
    } else {
      # Each row contains both pointwise (coverage/mean_width) and
      # simultaneous (coverage_joint/mean_width_joint) metrics.
      # Use method="fastfmm" for pointwise, method_joint="fastfmm_sim".
      # This mirrors the competitors extension schema (coverage + coverage_joint
      # coexist in the same row; method identifies the primary CI type).
      result <- metrics |>
        dplyr::mutate(
          method = "fastfmm",
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
          fit_time = fit_time,
          converged = TRUE,
          error_msg = NA_character_
        )
    }

    atomic_saveRDS_ff(result, save_path)
    pair_results[[grid_label]] <- result
    gc()
  }

  dplyr::bind_rows(pair_results)
}

# Study 1 runner (non-Gaussian) -----------------------------------------------

#' Run one (DGP, rep) — fastFMM extension for Study 1 (Poisson/Binomial)
#'
#' WARNING: analytic=FALSE with 500 bootstraps is expensive.
#' Do NOT run locally; use LRZ cluster.
#'
#' @param row DGP settings list (from make_study1_settings()).
#' @param rep_id Integer replicate id.
#' @param output_dir Output directory.
#' @param alpha Significance level.
#' @returns Tibble with fastFMM metrics for this (DGP, rep).
run_one_rep_fastfmm1 <- function(
  row,
  rep_id,
  output_dir,
  alpha = FASTFMM_ALPHA
) {
  if (inherits(row, "data.frame")) row <- as.list(row)
  seed <- STUDY1_BASE_SEED + 1000L * row$dgp_id + rep_id

  file_key <- sprintf("dgp%03d_rep%03d", row$dgp_id, rep_id)
  save_path <- file.path(output_dir, paste0(file_key, ".rds"))

  if (file.exists(save_path)) {
    obj <- tryCatch(readRDS(save_path), error = function(e) NULL)
    if (!is.null(obj) && nrow(obj) > 0) return(obj)
  }

  # Simulate — same DGP as Study 1 production
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
        "simulate failed dgp=%d rep=%d: %s",
        row$dgp_id,
        rep_id,
        conditionMessage(e)
      ))
      NULL
    }
  )
  if (is.null(sim)) {
    result <- null_row_study1(row, rep_id, seed, "simulate failed")
    atomic_saveRDS_ff(result, save_path)
    return(result)
  }

  t0 <- Sys.time()
  fit <- fit_fui(
    sim,
    family = row$family,
    n_boots = 500L,
    fui_seed = seed %% 1000L + 1L
  )
  fit_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(fit)) {
    result <- null_row_study1(row, rep_id, seed, "fui() returned NULL")
    result$fit_time <- fit_time
    atomic_saveRDS_ff(result, save_path)
    return(result)
  }

  metrics <- tryCatch(
    extract_fui_metrics(
      fit,
      sim,
      alpha = alpha,
      comparable_terms = FASTFMM_COMPARABLE_TERMS_STUDY1
    ),
    error = function(e) {
      warning(sprintf(
        "extract_fui_metrics Study1 dgp=%d rep=%d: %s",
        row$dgp_id,
        rep_id,
        conditionMessage(e)
      ))
      tibble()
    }
  )

  if (nrow(metrics) == 0) {
    result <- null_row_study1(
      row,
      rep_id,
      seed,
      "extract_fui_metrics returned 0 rows"
    )
    result$fit_time <- fit_time
  } else {
    result <- metrics |>
      dplyr::mutate(
        method = "fastfmm",
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

  atomic_saveRDS_ff(result, save_path)
  gc()
  result
}

# Main runner functions --------------------------------------------------------

#' Run fastFMM Study 2 extension
#'
#' @param n_rep Number of replicates per DGP cell.
#' @param grid_labels Grid labels to run (default: all STUDY2_GRIDS).
#' @param output_dir Output directory.
#' @param alpha Significance level.
#' @param dgp_ids_limit Optional integer vector to restrict DGP ids (smoke mode).
#' @returns Combined results tibble.
run_fastfmm_study2 <- function(
  n_rep = 50L,
  grid_labels = names(STUDY2_GRIDS),
  output_dir = FASTFMM_STUDY2_DIR,
  alpha = FASTFMM_ALPHA,
  dgp_ids_limit = NULL
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  settings <- make_study2_settings()
  if (!is.null(dgp_ids_limit)) {
    settings <- settings |> dplyr::filter(dgp_id %in% dgp_ids_limit)
  }

  task_grid <- settings |> tidyr::crossing(rep_id = seq_len(n_rep))

  # Skip already-complete (DGP, rep) pairs
  nygrid_first <- study2_parse_nygrid(grid_labels[1])
  existing_keys <- sub(
    "\\.rds$",
    "",
    list.files(output_dir, pattern = "^dgp\\d+.*\\.rds$")
  )
  grid_keys_first <- sprintf(
    "dgp%03d_n%03d_y%03d_rep%03d",
    task_grid$dgp_id,
    task_grid$n,
    nygrid_first,
    task_grid$rep_id
  )
  already_done <- grid_keys_first %in% existing_keys
  if (any(already_done)) {
    cat("  Skipping", sum(already_done), "already-completed (dgp, rep) pairs\n")
    task_grid <- task_grid[!already_done, , drop = FALSE]
  }

  cat("fastFMM Study 2 Extension\n")
  cat("=========================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("Reps per cell:", n_rep, "\n")
  cat("Grid labels:", paste(grid_labels, collapse = ", "), "\n")
  cat("Pending pairs:", nrow(task_grid), "\n")
  cat("Output dir:", output_dir, "\n\n")

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
      run_one_pair_fastfmm2(
        row,
        row$rep_id,
        grid_labels,
        output_dir,
        alpha = alpha
      ),
      error = function(e) {
        message("\nError: ", e$message)
        NULL
      }
    )
    if (!is.null(res)) results[[i]] <- res
  }
  cat("\n")

  dplyr::bind_rows(results)
}

#' Run fastFMM Study 1 extension (non-Gaussian)
#'
#' EXPENSIVE: uses bootstrap (n_boots=500) at each location × rep × DGP.
#' Intended for LRZ cluster only.
#'
#' @param n_rep Replicates per DGP cell.
#' @param output_dir Output directory.
#' @param alpha Significance level.
#' @param dgp_ids_limit Optional integer vector to restrict DGPs (smoke mode).
#' @returns Combined results tibble.
run_fastfmm_study1 <- function(
  n_rep = 150L,
  output_dir = FASTFMM_STUDY1_DIR,
  alpha = FASTFMM_ALPHA,
  dgp_ids_limit = NULL
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  settings <- make_study1_settings()
  if (!is.null(dgp_ids_limit)) {
    settings <- settings |> dplyr::filter(dgp_id %in% dgp_ids_limit)
  }

  grid <- settings |> tidyr::crossing(rep_id = seq_len(n_rep))
  existing_keys <- sub(
    "\\.rds$",
    "",
    list.files(output_dir, pattern = "^dgp\\d+_rep\\d+\\.rds$")
  )
  grid_keys <- sprintf("dgp%03d_rep%03d", grid$dgp_id, grid$rep_id)
  already_done <- grid_keys %in% existing_keys
  if (any(already_done)) {
    cat("  Skipping", sum(already_done), "already-completed reps\n")
    grid <- grid[!already_done, , drop = FALSE]
  }

  cat("fastFMM Study 1 Extension (non-Gaussian — intended for LRZ)\n")
  cat("===========================================================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("Reps per cell:", n_rep, "\n")
  cat("Pending reps:", nrow(grid), "\n")
  cat("Output dir:", output_dir, "\n\n")

  results <- list()
  for (i in seq_len(nrow(grid))) {
    row <- as.list(grid[i, ])
    cat(sprintf(
      "\r[%d/%d] dgp=%d (family=%s, corr=%s, n=%d) rep=%d",
      i,
      nrow(grid),
      row$dgp_id,
      row$family,
      row$corr_type,
      row$n,
      row$rep_id
    ))
    res <- tryCatch(
      run_one_rep_fastfmm1(row, row$rep_id, output_dir, alpha = alpha),
      error = function(e) {
        message("\nError: ", e$message)
        NULL
      }
    )
    if (!is.null(res)) results[[i]] <- res
  }
  cat("\n")
  dplyr::bind_rows(results)
}

# Summary helper --------------------------------------------------------------

summarize_fastfmm <- function(results) {
  results |>
    dplyr::filter(!is.na(coverage) | !is.na(coverage_joint)) |>
    dplyr::group_by(term_type, n, snr, corr_type, grid_label, nygrid) |>
    dplyr::summarise(
      mean_coverage = mean(coverage, na.rm = TRUE),
      mean_width_pw = mean(mean_width, na.rm = TRUE),
      mean_coverage_joint = mean(coverage_joint, na.rm = TRUE),
      mean_width_joint = mean(mean_width_joint, na.rm = TRUE),
      n_reps = dplyr::n(),
      n_converged = sum(converged, na.rm = TRUE),
      .groups = "drop"
    )
}

# Main Entry Point ------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args) >= 1) args[1] else "smoke"
  study <- if (length(args) >= 2) args[2] else "study2"

  study <- match.arg(study, c("study2", "study1"))

  n_rep <- switch(
    mode,
    smoke = 1L,
    pilot = 10L,
    full = 50L,
    as.integer(mode)
  )
  if (is.na(n_rep)) n_rep <- 1L

  cat(sprintf(
    "fastFMM Extension: mode=%s, study=%s, n_rep=%d\n\n",
    mode,
    study,
    n_rep
  ))

  if (study == "study2") {
    # Smoke: single smallest DGP (dgp_id=1: iid, n=20, snr=25)
    dgp_limit <- if (mode == "smoke") 1L else NULL
    grid_labs <- if (mode == "smoke") names(STUDY2_GRIDS)[1] else
      names(STUDY2_GRIDS)

    results <- run_fastfmm_study2(
      n_rep = n_rep,
      grid_labels = grid_labs,
      dgp_ids_limit = dgp_limit
    )

    if (nrow(results) > 0) {
      cat("\n========== fastFMM Study 2 RESULTS PREVIEW ==========\n")
      pw_cols <- c(
        "term_type",
        "coverage",
        "mean_width",
        "converged",
        "grid_label",
        "corr_type",
        "n",
        "snr"
      )
      joint_cols <- c(
        "term_type",
        "coverage_joint",
        "mean_width_joint",
        "grid_label",
        "corr_type",
        "n",
        "snr"
      )

      avail_pw <- intersect(pw_cols, names(results))
      avail_joint <- intersect(joint_cols, names(results))

      cat("Pointwise (fastfmm):\n")
      print(
        results |>
          dplyr::select(dplyr::any_of(avail_pw)) |>
          dplyr::filter(!is.na(coverage)),
        n = 30
      )

      cat("\nSimultaneous/joint (fastfmm_sim):\n")
      print(
        results |>
          dplyr::select(dplyr::any_of(avail_joint)) |>
          dplyr::filter(!is.na(coverage_joint)),
        n = 30
      )

      # Validation checks
      cat("\n--- Validation ---\n")
      cat(
        "coverage in [0,1]:",
        all(results$coverage >= 0 & results$coverage <= 1, na.rm = TRUE),
        "\n"
      )
      cat(
        "joint in [0,1]:",
        all(
          results$coverage_joint >= 0 & results$coverage_joint <= 1,
          na.rm = TRUE
        ),
        "\n"
      )
      cat("mean_width_joint >= mean_width (joint >= pointwise width):\n")
      width_check <- results |>
        dplyr::filter(!is.na(mean_width) & !is.na(mean_width_joint)) |>
        dplyr::summarise(all_ok = all(mean_width_joint >= mean_width * 0.999))
      cat(" ", width_check$all_ok, "\n")
      cat(
        "Terms in output:",
        paste(sort(unique(results$term_type)), collapse = ", "),
        "\n"
      )
      cat("Files written to:", FASTFMM_STUDY2_DIR, "\n")
      cat(
        "Files:",
        length(list.files(FASTFMM_STUDY2_DIR, pattern = "\\.rds$")),
        "\n"
      )
    } else {
      cat("No results returned.\n")
    }
  } else if (study == "study1") {
    # Study 1 non-Gaussian: IMPLEMENT but do NOT run in smoke here locally
    # (too slow without bootstrap=500). Only run in full mode on LRZ.
    if (mode == "smoke") {
      cat("Study 1 non-Gaussian fastFMM: smoke mode runs 1 rep for 1 DGP.\n")
      cat("WARNING: bootstrap (n_boots=500) is very slow locally.\n")
      cat("For full runs, use LRZ cluster.\n\n")
      dgp_limit <- 1L
    } else {
      dgp_limit <- NULL
    }
    results <- run_fastfmm_study1(
      n_rep = n_rep,
      dgp_ids_limit = dgp_limit
    )
    if (nrow(results) > 0) {
      cat("\n========== fastFMM Study 1 RESULTS PREVIEW ==========\n")
      print(
        results |>
          dplyr::select(
            dplyr::any_of(c(
              "term_type",
              "coverage",
              "coverage_joint",
              "mean_width",
              "converged",
              "family",
              "corr_type",
              "n"
            ))
          ),
        n = 20
      )
    }
  }

  cat("\nfastFMM extension complete.\n")
}
