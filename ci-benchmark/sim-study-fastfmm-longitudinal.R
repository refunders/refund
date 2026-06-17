# ===========================================================================
# E4: fastFMM / FUI comparison — Repeated-measures (longitudinal) sub-study
# ===========================================================================
#
# PAPER CONTEXT
# Study 1/2 DGPs have one curve per subject (independent curves); the subject
# random intercept (1|id) is unidentifiable there and fastFMM/FUI cannot be
# applied. This sub-study uses a *longitudinal* DGP with J repeated
# functional observations per subject, the setting FUI was designed for.
#
# ADEMP
#   Aims: Compare pffr (Bayesian/sandwich CIs) vs fastFMM::fui (analytic LME
#         CIs) for coefficient function coverage on longitudinal functional
#         data with known truth.
#   Data-generating mechanisms (DGP):
#     Y_ij(t) = beta0(t) + zlin_ij * beta1(t) + b_i + eps_ij(t)
#     - N subjects, J visits each.
#     - zlin_ij ~ N(0,1), then CENTERED per rep (mean subtracted).
#       Centering ensures mgcv's mean-centering of the design matrix
#       does not shift the intercept, so truth alignment is clean.
#     - b_i ~ N(0, sigma_b^2): SCALAR random intercept per subject,
#       constant over t. This is the direct analogue of (1|id) in FUI,
#       and is what pffr's c(s(id, bs="re")) models.
#       NOTE: We do NOT use a functional b_i(t) because (a) FUI only models
#       a scalar RE (1|id), and (b) functional RE requires s(id, bs="re")
#       in pffr (no c() wrapper) which creates a tensor product and causes
#       centering identifiability issues. Using a scalar b_i gives the
#       apples-to-apples comparison the study aims for.
#     - eps_ij(t) ~ N(0, sigma_eps^2 * I_L): iid Gaussian pointwise noise.
#     - Response grid L points on [0,1]; Gaussian family.
#     - beta0(t) and beta1(t) are random B-spline functions, CENTERED to
#       mean zero over t. Centering is essential: pffr's Intercept(yindex)
#       smooth has a sum-to-zero-like constraint (in the sense that the mean
#       goes into the parametric intercept, which varies per rep/sample). By
#       centering beta0, the Intercept(yindex) smooth = beta0(t) directly.
#       Similarly for beta1 and zlin centering.
#   Estimands: beta0(t), beta1(t) on the response grid (known truth).
#   Methods:
#     pffr  — scalar random intercept via c(s(id, bs="re")); Bayesian
#             default CIs + cluster-sandwich + cl2 CIs (cluster = by id).
#             The c() wrapper makes the RE constant over t (scalar shift per
#             subject), directly comparable to (1|id) in FUI.
#             Full intercept estimate = pterms[1,1] + Intercept(yindex).
#     fastFMM — fui(Y ~ zlin + (1|id), analytic=TRUE); pointwise CIs from
#             betaHat_var diagonal; simultaneous via compute_qn_at_alpha().
#   Performance measures:
#     Pointwise coverage + mean CI width at 90% nominal (alpha=0.10).
#     Simultaneous joint coverage + band width.
#     Monte Carlo SEs: sqrt(p(1-p)/R).
#
# DGP CELLS (factorial):
#   N in {50, 100}, J in {3, 5} → 4 cells
#   Fixed: L=60, sigma_b=1, sigma_eps=0.5
#   (Full run would expand N/J range and add more reps.)
#
# SEED SCHEME
#   Base: STUDY4_BASE_SEED = 5001
#   Per (cell, rep): seed = STUDY4_BASE_SEED + 1000 * cell_id + rep_id
#   cell_id: row index of make_study4_settings() (1-based).
#
# PFFR SYNTAX (verified):
#   pffr(Y ~ zlin + c(s(id, bs="re")), yind=t_grid, data=long_df)
#   - c() wrapper: makes the RE constant over t (= scalar per-subject shift).
#   - Appears in coef() as smterms: Intercept(yindex), zlin(yindex), id(yindex).
#   - Full intercept = pterms[1,1] (scalar) + smterms[[1]] (functional smooth).
#     We compare only the functional smooth against centered beta0(t).
#     (When beta0 and zlin are centered, the scalar pterms ≈ 0 is the grand mean.)
#
# PFFR SANDWICH CIs — IMPORTANT LIMITATION:
#   pffr's built-in cluster/cl2 sandwich ALWAYS clusters by curve (i.e.,
#   each of the N*J observation rows = one cluster). In the longitudinal
#   setting the correct unit is the SUBJECT (N clusters of J curves each),
#   but pffr has no built-in way to cluster by subject. As a result:
#   - "pffr_default" / "pffr_cluster": cluster-by-curve SE ≈ Bayesian SE
#     (near-identical because clustering a single observation = no sandwich
#     correction). This is INCORRECT for the longitudinal setting (SE too
#     small by ~3x for the linear term).
#   - "pffr_cl2": bias-corrected cluster-by-curve; still underestimates but
#     less severely (~2x).
#   This is a KEY SCIENTIFIC FINDING of this sub-study: pffr's current
#   sandwich infrastructure is not designed for multi-observation-per-curve
#   (repeated-measures) data. Custom subject-level clustering would require
#   modifying build_cluster_id() to group by subject rather than by curve.
#   This is recorded as a finding, not fixed here.
#
# PFFR INTERCEPT TRUTH (centering note):
#   pffr's Intercept(yindex) smooth is the "shape" of beta0(t) around the
#   per-rep grand mean of Y. When beta0(t) has mean 0 and zlin is centered,
#   the Intercept(yindex) estimate ≈ beta0(t) directly (the scalar parametric
#   intercept ≈ 0). We compare against sim$truth$beta$intercept (which is
#   already centered = mean zero by construction).
#
# FASTFMM TERM MAPPING:
#   fui row "(Intercept)" -> "intercept" -> truth$beta$intercept = beta0(t)
#   fui row "zlin"        -> "linear"    -> truth$beta$zlin = beta1(t)
#   (1|id) is the scalar random effect — not an estimand.
#
# DATA FORMAT:
#   Long format: N*J rows, each row one (subject i, visit j) observation.
#   Columns: id (factor), visit (integer), zlin (numeric, centered), Y (1×L).
#   pffr stacks this as N*J curves with the scalar RE term.
#   fui receives the same long data frame (its native format).
#
# REUSED HELPERS (from sim-study-fastfmm-extension.R):
#   extract_fui_metrics(), compute_qn_at_alpha(), atomic_saveRDS_ff().
#
# OUTPUT:
#   ci-benchmark/study4-longitudinal/
#   File naming: cell%03d_rep%03d.rds (per-cell, per-rep atomic saves).
#   One RDS = tibble with method × term_type rows for one (cell, rep).
#
# USAGE:
#   Rscript ci-benchmark/sim-study-fastfmm-longitudinal.R [mode]
#   mode: "smoke" (2 reps, 1 cell), "pilot" (10 reps, all cells),
#         "full" (50 reps, all cells — do NOT run locally)
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

# Reuse the atomic save and FUI helpers from the existing fastFMM extension.
# Its main block only runs when it is the entry point (sys.nframe() == 0),
# so sourcing is safe here.
source("ci-benchmark/sim-study-fastfmm-extension.R")

# Constants -------------------------------------------------------------------

STUDY4_BASE_SEED <- 5001L
STUDY4_ALPHA <- 0.10 # 90% CIs (nominal 1-alpha)
STUDY4_L <- 60L # Response grid length on [0,1]
STUDY4_SIGMA_B <- 1.0 # Subject-level scalar RE SD
STUDY4_SIGMA_EPS <- 0.5 # Pointwise noise SD
STUDY4_K_BETA <- 6L # B-spline basis dimension for beta0, beta1
STUDY4_OUTPUT_DIR <- "ci-benchmark/study4-longitudinal"

# Comparable terms for E4 (only population-level fixed-effect coefficients)
STUDY4_COMPARABLE_TERMS <- c("intercept", "linear")

# DGP settings ----------------------------------------------------------------

#' Build factorial DGP settings for Study 4
#'
#' 2×2 factorial: N in {50, 100}, J in {3, 5}.
#' Each row is one DGP cell identified by cell_id.
#'
#' @returns Tibble with columns: cell_id, N, J.
make_study4_settings <- function() {
  grid <- tidyr::crossing(N = c(50L, 100L), J = c(3L, 5L)) |>
    dplyr::arrange(N, J) |>
    dplyr::mutate(cell_id = dplyr::row_number())
  grid
}

# DGP simulation --------------------------------------------------------------

#' Generate zero-mean smooth 1D function via random B-spline coefficients
#'
#' Generates a random function then subtracts its mean over the grid so
#' the result has mean zero. This centering is critical: pffr's
#' Intercept(yindex) smooth must be compared against a zero-mean truth;
#' the non-zero mean is absorbed into the scalar parametric intercept.
#' Similarly, a zero-mean beta1 ensures that centering zlin does not shift
#' the effective intercept.
#'
#' @param grid Evaluation grid (numeric vector).
#' @param k Number of B-spline basis functions.
#' @param amplitude Target SD of the resulting function.
#' @returns Numeric vector, length(grid), with mean ≈ 0.
generate_smooth_1d_zeromean <- function(
  grid,
  k = STUDY4_K_BETA,
  amplitude = 0.3
) {
  B <- splines::bs(grid, df = k, intercept = TRUE)
  coefs <- stats::rnorm(k)
  f <- as.vector(B %*% coefs)
  # Center to mean zero
  f <- f - mean(f)
  sd_f <- stats::sd(f)
  if (sd_f > 0) f <- f / sd_f * amplitude
  f
}

#' Simulate one longitudinal functional dataset (Study 4 DGP)
#'
#' DGP: Y_ij(t) = beta0(t) + zlin_ij * beta1(t) + b_i + eps_ij(t)
#' where:
#'   - beta0(t), beta1(t): zero-mean B-spline functions (centering design).
#'   - b_i ~ N(0, sigma_b^2): scalar per-subject random intercept.
#'   - zlin_ij ~ N(0,1), centered across all N*J observations per rep.
#'   - eps_ij(t) iid ~ N(0, sigma_eps^2).
#'
#' Centering beta0, beta1, and zlin ensures:
#'   (a) pffr's Intercept(yindex) smooth ≈ beta0(t) (zero-mean truth).
#'   (b) pffr's zlin(yindex) smooth ≈ beta1(t).
#'   (c) The scalar parametric intercept ≈ grand_mean(Y) ≈ 0.
#'
#' @param N Number of subjects.
#' @param J Number of visits per subject.
#' @param L Response grid length.
#' @param sigma_b SD of scalar subject-level RE.
#' @param sigma_eps Pointwise noise SD.
#' @param seed RNG seed (integer).
#' @returns List with: data (long data.frame), t_grid, truth (list with
#'   $beta = list(intercept, zlin)), N, J, L, sigma_b.
simulate_longitudinal <- function(
  N,
  J,
  L = STUDY4_L,
  sigma_b = STUDY4_SIGMA_B,
  sigma_eps = STUDY4_SIGMA_EPS,
  seed
) {
  set.seed(seed)
  t_grid <- seq(0, 1, length.out = L)

  # --- Truth: zero-mean coefficient functions (vary per rep via seed) ------
  # Both are centered to mean zero so pffr intercept/linear smterms ≈ truth.
  beta0 <- generate_smooth_1d_zeromean(
    t_grid,
    k = STUDY4_K_BETA,
    amplitude = 0.5
  )
  beta1 <- generate_smooth_1d_zeromean(
    t_grid,
    k = STUDY4_K_BETA,
    amplitude = 0.4
  )

  # --- Subject-level scalar random effects (constant over t) ---------------
  b_i <- stats::rnorm(N, mean = 0, sd = sigma_b)

  # --- Observation-level data in long format --------------------------------
  n_obs <- N * J
  id_vec <- rep(seq_len(N), each = J)
  visit_vec <- rep(seq_len(J), times = N)
  zlin_raw <- stats::rnorm(n_obs)
  # Center zlin across all observations: removes mgcv's automatic centering
  # which otherwise absorbs mean(zlin)*beta1(t) into the intercept.
  zlin_vec <- zlin_raw - mean(zlin_raw)

  # Assemble Y matrix (n_obs x L)
  Y_mat <- matrix(0, nrow = n_obs, ncol = L)
  for (obs in seq_len(n_obs)) {
    i <- id_vec[obs]
    mu_ij <- beta0 + zlin_vec[obs] * beta1 + b_i[i]
    Y_mat[obs, ] <- mu_ij + stats::rnorm(L, mean = 0, sd = sigma_eps)
  }

  long_df <- data.frame(
    id = factor(id_vec),
    visit = visit_vec,
    zlin = zlin_vec
  )
  long_df$Y <- I(Y_mat)

  list(
    data = long_df,
    t_grid = t_grid,
    truth = list(
      beta = list(
        intercept = beta0, # zero-mean; compared against Intercept(yindex) smooth
        zlin = beta1 # zero-mean; compared against zlin(yindex) smooth
      )
    ),
    N = N,
    J = J,
    L = L,
    sigma_b = sigma_b
  )
}

# pffr fitting and metric extraction ------------------------------------------

#' Fit pffr with scalar constant random intercept on longitudinal data
#'
#' Model: Y ~ zlin + c(s(id, bs="re")), yind=t_grid.
#' The c() wrapper makes s(id, bs="re") constant over t: each subject gets
#' a scalar random shift b_i (BLUP), not a functional curve. This is the
#' correct pffr analogue of (1|id) in FUI.
#'
#' NOTE on identifiability: the RE is identified from within-subject
#' replication (J visits per subject). With J >= 2, the scalar b_i is
#' identifiable separately from the fixed-effect parameters.
#'
#' Sandwich types extracted:
#'   "none"    — default Bayesian/REML SEs (Wood 2006)
#'   "cluster" — cluster-robust sandwich SE (cluster = id)
#'   "cl2"     — bias-corrected cluster-robust (Bell-McCaffrey)
#'
#' @param sim Simulation result from simulate_longitudinal().
#' @returns pffr fit object or NULL on error.
fit_pffr_longitudinal <- function(sim) {
  df <- sim$data
  t_grid <- sim$t_grid

  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      pffr(
        Y ~ zlin + c(s(id, bs = "re")),
        yind = t_grid,
        data = df
      )
    )),
    error = function(e) {
      warning(sprintf("pffr() failed: %s", conditionMessage(e)))
      NULL
    }
  )
  fit
}

#' Extract pffr CI metrics for intercept or linear term
#'
#' For the intercept: pffr reports a scalar parametric intercept (pterms)
#' PLUS a functional smooth Intercept(yindex). The estimand is beta0(t) which,
#' by construction in simulate_longitudinal(), has mean zero. With mean-zero
#' beta0 and centered zlin, the Intercept(yindex) smooth ≈ beta0(t) directly
#' and pterms ≈ grand_mean(Y) ≈ 0.
#'
#' DESIGN CHOICE: We compare pffr's Intercept(yindex) smooth against the
#' centered truth (beta0, already mean-zero). The scalar pterms (which
#' estimates the grand mean of Y) is NOT included in the comparison because:
#'   (a) it conflates the fixed-effect intercept with the random effect mean
#'       (sum of b_i / N is not guaranteed zero in finite samples), and
#'   (b) fui's (Intercept) also estimates the population intercept function,
#'       not the sample grand mean — so the comparison is symmetric.
#'
#' For the linear term: pffr's zlin(yindex) smooth vs truth beta1(t).
#'
#' The SEs from coef.pffr include the seWithMean=TRUE uncertainty (which
#' incorporates the uncertainty in the mean), but with mean-zero functions
#' this makes negligible difference.
#'
#' @param fit pffr fit object.
#' @param sim Simulation result (for truth, t_grid).
#' @param term_type One of "intercept", "linear".
#' @param sandwich_type One of "none", "cluster", "cl2".
#' @param alpha Significance level (default STUDY4_ALPHA).
#' @returns One-row tibble with coverage, mean_width, etc., or NULL.
extract_pffr_metrics_long <- function(
  fit,
  sim,
  term_type,
  sandwich_type = "none",
  alpha = STUDY4_ALPHA
) {
  if (is.null(fit)) return(NULL)

  t_grid <- sim$t_grid

  coefs <- tryCatch(
    coef(
      fit,
      sandwich = sandwich_type,
      seWithMean = FALSE,
      n1 = 50,
      n2 = 25,
      n3 = 15
    ),
    error = function(e) {
      warning(sprintf(
        "coef.pffr(%s, sandwich=%s): %s",
        term_type,
        sandwich_type,
        conditionMessage(e)
      ))
      NULL
    }
  )
  if (is.null(coefs) || is.null(coefs$smterms)) return(NULL)

  sm_names <- names(coefs$smterms)

  if (term_type == "intercept") {
    # Intercept(yindex) smooth — compared against zero-mean beta0(t)
    idx <- which(grepl("^intercept\\(", tolower(sm_names)))
    if (length(idx) == 0) return(NULL)
    info <- coefs$smterms[[idx[1]]]
    truth_vals <- sim$truth$beta$intercept # zero-mean by construction
  } else if (term_type == "linear") {
    # zlin(yindex) smooth — compared against zero-mean beta1(t)
    idx <- which(grepl("zlin", tolower(sm_names)))
    if (length(idx) == 0) return(NULL)
    info <- coefs$smterms[[idx[1]]]
    truth_vals <- sim$truth$beta$zlin # zero-mean by construction
  } else {
    return(NULL)
  }

  est <- info$coef$value # on 100-pt coef.pffr grid
  se <- info$coef$se
  x_eval <- info$x # 100 evaluation points on [0,1]
  if (is.null(est) || is.null(se) || is.null(x_eval)) return(NULL)

  # Interpolate truth onto pffr's evaluation grid
  truth_on_grid <- stats::approx(t_grid, truth_vals, xout = x_eval, rule = 2)$y
  if (length(truth_on_grid) != length(est)) return(NULL)

  crit <- stats::qnorm(1 - alpha / 2)
  lower <- est - crit * se
  upper <- est + crit * se
  covered <- (truth_on_grid >= lower) & (truth_on_grid <= upper)

  tibble(
    term_type = term_type,
    coverage = mean(covered, na.rm = TRUE),
    mean_width = mean(upper - lower, na.rm = TRUE),
    coverage_joint = NA_real_, # pffr simultaneous CIs handled in E2
    mean_width_joint = NA_real_,
    rmse = sqrt(mean((est - truth_on_grid)^2, na.rm = TRUE)),
    bias = mean(est - truth_on_grid, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE),
    n_grid = sum(!is.na(covered))
  )
}

# fui fitting and metric extraction -------------------------------------------

#' Build fui-compatible data frame from longitudinal simulation result
#'
#' fui requires a data.frame with:
#'   - Response: N*J × L matrix in a single column named "Y".
#'   - id: integer grouping variable identifying subjects.
#'   - Scalar covariates: standard columns (zlin here).
#'
#' @param sim Simulation result from simulate_longitudinal().
#' @returns Data frame suitable for fui().
build_fui_data_long <- function(sim) {
  df <- sim$data
  # fui needs a plain (not I()-wrapped) matrix column
  df$Y <- as.matrix(df$Y)
  df$id <- as.integer(df$id)
  df
}

#' Fit fui() on longitudinal simulation result
#'
#' Formula: Y ~ zlin + (1|id)
#' analytic=TRUE: Gaussian-analytic path (Cui et al. 2022 Algorithm 1).
#' var=TRUE: compute betaHat_var for CI construction.
#'
#' @param sim Simulation result from simulate_longitudinal().
#' @returns fui fit object or NULL on error.
fit_fui_longitudinal <- function(sim) {
  fui_df <- build_fui_data_long(sim)
  frml <- Y ~ zlin + (1 | id)

  fit <- tryCatch(
    suppressWarnings(
      fui(
        formula = frml,
        data = fui_df,
        family = "gaussian",
        analytic = TRUE,
        var = TRUE,
        parallel = FALSE,
        silent = TRUE
      )
    ),
    error = function(e) {
      warning(sprintf("fui() failed: %s", conditionMessage(e)))
      NULL
    }
  )
  fit
}

#' Extract fui CI metrics for Study 4 (longitudinal)
#'
#' Adapts extract_fui_metrics() (from sim-study-fastfmm-extension.R) to the
#' Study 4 DGP. The truth$beta keys match extract_fui_metrics() conventions:
#'   "intercept" -> truth$beta$intercept = beta0(t)
#'   "zlin"      -> truth$beta$zlin      = beta1(t) (mapped to type "linear")
#'
#' @param fit fui fit object.
#' @param sim Simulation result.
#' @param alpha Significance level.
#' @returns Tibble (term_type, coverage, mean_width, coverage_joint, mean_width_joint).
extract_fui_metrics_long <- function(fit, sim, alpha = STUDY4_ALPHA) {
  extract_fui_metrics(
    fit = fit,
    sim = sim,
    alpha = alpha,
    comparable_terms = STUDY4_COMPARABLE_TERMS,
    qn_seed = 42L
  )
}

# Null-row helper -------------------------------------------------------------

#' Build a null result row for one (cell, rep) when a fit fails
#'
#' @param cell DGP settings list.
#' @param rep_id Integer replicate id.
#' @param seed Seed used.
#' @param method Character method label.
#' @param error_msg Character error message.
#' @returns Tibble with NA metric columns and metadata.
null_row_study4 <- function(
  cell,
  rep_id,
  seed,
  method,
  error_msg = NA_character_
) {
  tibble(
    term_type = STUDY4_COMPARABLE_TERMS,
    coverage = NA_real_,
    mean_width = NA_real_,
    coverage_joint = NA_real_,
    mean_width_joint = NA_real_,
    rmse = NA_real_,
    bias = NA_real_,
    mean_se = NA_real_,
    n_grid = NA_integer_,
    method = method,
    cell_id = cell$cell_id,
    rep_id = rep_id,
    seed = seed,
    N = cell$N,
    J = cell$J,
    n_obs = cell$N * cell$J,
    fit_time = NA_real_,
    converged = FALSE,
    error_msg = error_msg
  )
}

# Per-(cell, rep) runner ------------------------------------------------------

#' Run one (cell, rep) for Study 4 — both pffr and fui
#'
#' Generates data once, fits pffr (3 sandwich variants) and fui, extracts
#' coverage metrics for each method. Saves atomically to output_dir.
#' Skips if the file already exists and is non-empty.
#'
#' @param cell Named list: cell_id, N, J (from make_study4_settings()).
#' @param rep_id Integer replicate id.
#' @param output_dir Character path.
#' @param alpha Significance level.
#' @returns Tibble combining all method rows for this (cell, rep).
run_one_rep_study4 <- function(
  cell,
  rep_id,
  output_dir = STUDY4_OUTPUT_DIR,
  alpha = STUDY4_ALPHA
) {
  if (inherits(cell, "data.frame")) cell <- as.list(cell)
  seed <- STUDY4_BASE_SEED + 1000L * cell$cell_id + rep_id

  file_key <- sprintf("cell%03d_rep%03d", cell$cell_id, rep_id)
  save_path <- file.path(output_dir, paste0(file_key, ".rds"))

  # Resume: skip if already completed
  if (file.exists(save_path)) {
    obj <- tryCatch(readRDS(save_path), error = function(e) NULL)
    if (!is.null(obj) && nrow(obj) > 0) return(obj)
  }

  all_methods <- c("pffr_default", "pffr_cluster", "pffr_cl2", "fastfmm")

  # --- Simulate data --------------------------------------------------------
  sim <- tryCatch(
    simulate_longitudinal(
      N = cell$N,
      J = cell$J,
      L = STUDY4_L,
      sigma_b = STUDY4_SIGMA_B,
      sigma_eps = STUDY4_SIGMA_EPS,
      seed = seed
    ),
    error = function(e) {
      warning(sprintf(
        "simulate_longitudinal failed cell=%d rep=%d: %s",
        cell$cell_id,
        rep_id,
        conditionMessage(e)
      ))
      NULL
    }
  )

  if (is.null(sim)) {
    result <- dplyr::bind_rows(lapply(all_methods, function(m) {
      null_row_study4(cell, rep_id, seed, m, "simulate_longitudinal failed")
    }))
    atomic_saveRDS_ff(result, save_path)
    return(result)
  }

  results <- list()

  # --- pffr fit (shared across 3 sandwich variants) -------------------------
  t0_pffr <- Sys.time()
  pffr_fit <- fit_pffr_longitudinal(sim)
  pffr_time <- as.numeric(difftime(Sys.time(), t0_pffr, units = "secs"))

  sandwich_types <- list(
    pffr_default = "none",
    pffr_cluster = "cluster",
    pffr_cl2 = "cl2"
  )

  for (method_nm in names(sandwich_types)) {
    sw <- sandwich_types[[method_nm]]

    if (is.null(pffr_fit)) {
      r <- null_row_study4(
        cell,
        rep_id,
        seed,
        method_nm,
        "pffr() returned NULL"
      )
      r$fit_time <- pffr_time
      results[[method_nm]] <- r
      next
    }

    # Extract metrics for each comparable term
    term_rows <- lapply(STUDY4_COMPARABLE_TERMS, function(tt) {
      m <- extract_pffr_metrics_long(
        pffr_fit,
        sim,
        tt,
        sandwich_type = sw,
        alpha = alpha
      )
      if (is.null(m)) {
        return(tibble(
          term_type = tt,
          coverage = NA_real_,
          mean_width = NA_real_,
          coverage_joint = NA_real_,
          mean_width_joint = NA_real_,
          rmse = NA_real_,
          bias = NA_real_,
          mean_se = NA_real_,
          n_grid = NA_integer_
        ))
      }
      m
    })

    pffr_metrics <- dplyr::bind_rows(term_rows)
    results[[method_nm]] <- pffr_metrics |>
      dplyr::mutate(
        method = method_nm,
        cell_id = cell$cell_id,
        rep_id = rep_id,
        seed = seed,
        N = cell$N,
        J = cell$J,
        n_obs = cell$N * cell$J,
        fit_time = pffr_time,
        converged = TRUE,
        error_msg = NA_character_
      )
  }

  # --- fui fit --------------------------------------------------------------
  t0_fui <- Sys.time()
  fui_fit <- fit_fui_longitudinal(sim)
  fui_time <- as.numeric(difftime(Sys.time(), t0_fui, units = "secs"))

  if (is.null(fui_fit)) {
    r <- null_row_study4(cell, rep_id, seed, "fastfmm", "fui() returned NULL")
    r$fit_time <- fui_time
    results[["fastfmm"]] <- r
  } else {
    fui_metrics <- tryCatch(
      extract_fui_metrics_long(fui_fit, sim, alpha = alpha),
      error = function(e) {
        warning(sprintf(
          "extract_fui_metrics_long cell=%d rep=%d: %s",
          cell$cell_id,
          rep_id,
          conditionMessage(e)
        ))
        tibble()
      }
    )

    if (nrow(fui_metrics) == 0) {
      r <- null_row_study4(
        cell,
        rep_id,
        seed,
        "fastfmm",
        "extract_fui_metrics returned 0 rows"
      )
      r$fit_time <- fui_time
      results[["fastfmm"]] <- r
    } else {
      # add rmse/bias/mean_se/n_grid columns matching pffr schema if absent
      if (!"rmse" %in% names(fui_metrics)) fui_metrics$rmse <- NA_real_
      if (!"bias" %in% names(fui_metrics)) fui_metrics$bias <- NA_real_
      if (!"mean_se" %in% names(fui_metrics)) fui_metrics$mean_se <- NA_real_
      if (!"n_grid" %in% names(fui_metrics)) fui_metrics$n_grid <- NA_integer_

      results[["fastfmm"]] <- fui_metrics |>
        dplyr::mutate(
          method = "fastfmm",
          cell_id = cell$cell_id,
          rep_id = rep_id,
          seed = seed,
          N = cell$N,
          J = cell$J,
          n_obs = cell$N * cell$J,
          fit_time = fui_time,
          converged = TRUE,
          error_msg = NA_character_
        )
    }
  }

  result <- dplyr::bind_rows(results)
  atomic_saveRDS_ff(result, save_path)
  gc()
  result
}

# Main runner -----------------------------------------------------------------

#' Run Study 4 (E4 longitudinal sub-study)
#'
#' Iterates over all DGP cells × reps; saves per-(cell,rep) atomically;
#' resumes if output files already exist.
#'
#' @param n_rep Number of replicates per DGP cell.
#' @param output_dir Output directory path.
#' @param alpha Significance level.
#' @param cell_ids_limit Optional integer vector to restrict cells (smoke mode).
#' @returns Combined results tibble (all methods × terms × cells × reps).
run_study4 <- function(
  n_rep = 10L,
  output_dir = STUDY4_OUTPUT_DIR,
  alpha = STUDY4_ALPHA,
  cell_ids_limit = NULL
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  settings <- make_study4_settings()
  if (!is.null(cell_ids_limit)) {
    settings <- settings |> dplyr::filter(cell_id %in% cell_ids_limit)
  }

  task_grid <- settings |> tidyr::crossing(rep_id = seq_len(n_rep))

  # Skip already-complete (cell, rep) pairs
  existing_keys <- sub(
    "\\.rds$",
    "",
    list.files(output_dir, pattern = "^cell\\d+_rep\\d+\\.rds$")
  )
  task_keys <- sprintf("cell%03d_rep%03d", task_grid$cell_id, task_grid$rep_id)
  already_done <- task_keys %in% existing_keys
  if (any(already_done)) {
    cat(
      "  Skipping",
      sum(already_done),
      "already-completed (cell, rep) pairs\n"
    )
    task_grid <- task_grid[!already_done, , drop = FALSE]
  }

  cat("E4 Longitudinal Study (pffr vs fastFMM)\n")
  cat("========================================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("Reps per cell:", n_rep, "\n")
  cat("Pending pairs:", nrow(task_grid), "\n")
  cat("Output dir:", output_dir, "\n\n")
  cat(
    "Settings:\n",
    "  N:",
    paste(unique(settings$N), collapse = ", "),
    "\n",
    "  J:",
    paste(unique(settings$J), collapse = ", "),
    "\n",
    "  L:",
    STUDY4_L,
    " sigma_b:",
    STUDY4_SIGMA_B,
    " sigma_eps:",
    STUDY4_SIGMA_EPS,
    "\n\n"
  )

  results <- list()
  for (i in seq_len(nrow(task_grid))) {
    row <- as.list(task_grid[i, ])
    cat(sprintf(
      "\r[%d/%d] cell=%d (N=%d, J=%d), rep=%d   ",
      i,
      nrow(task_grid),
      row$cell_id,
      row$N,
      row$J,
      row$rep_id
    ))
    res <- tryCatch(
      run_one_rep_study4(row, row$rep_id, output_dir, alpha = alpha),
      error = function(e) {
        message("\nError in run_one_rep_study4: ", e$message)
        NULL
      }
    )
    if (!is.null(res)) results[[i]] <- res
  }
  cat("\n")
  dplyr::bind_rows(results)
}

# Summary helper --------------------------------------------------------------

#' Summarise Study 4 results over reps
#'
#' @param results Tibble from run_study4().
#' @returns Summary tibble with mean coverage, width, MC SE.
summarize_study4 <- function(results) {
  results |>
    dplyr::filter(!is.na(coverage)) |>
    dplyr::group_by(method, term_type, N, J) |>
    dplyr::summarise(
      mean_coverage = mean(coverage, na.rm = TRUE),
      mc_se_coverage = sqrt(
        mean(coverage, na.rm = TRUE) *
          (1 - mean(coverage, na.rm = TRUE)) /
          sum(!is.na(coverage))
      ),
      mean_width = mean(mean_width, na.rm = TRUE),
      mean_coverage_joint = mean(coverage_joint, na.rm = TRUE),
      mean_width_joint = mean(mean_width_joint, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      n_reps = dplyr::n(),
      n_converged = sum(converged, na.rm = TRUE),
      .groups = "drop"
    )
}

# Validation checks -----------------------------------------------------------

#' Print smoke/pilot validation checks on results
#'
#' @param results Tibble from run_study4().
#' @param nominal Nominal coverage (default 1 - STUDY4_ALPHA).
validate_study4 <- function(results, nominal = 1 - STUDY4_ALPHA) {
  cat("\n--- Validation (Study 4 E4 Longitudinal) ---\n")
  cat(
    "coverage in [0,1]:",
    all(results$coverage >= 0 & results$coverage <= 1, na.rm = TRUE),
    "\n"
  )
  cat(
    "mean_width > 0 (where non-NA):",
    all(results$mean_width[!is.na(results$mean_width)] > 0),
    "\n"
  )
  fui_rows <- dplyr::filter(
    results,
    method == "fastfmm",
    !is.na(mean_width_joint)
  )
  if (nrow(fui_rows) > 0) {
    cat(
      "fui joint width >= pointwise width:",
      all(
        fui_rows$mean_width_joint >= fui_rows$mean_width * 0.999,
        na.rm = TRUE
      ),
      "\n"
    )
  }
  cat(
    "Methods in output:",
    paste(sort(unique(results$method)), collapse = ", "),
    "\n"
  )
  cat(
    "Terms in output:",
    paste(sort(unique(results$term_type)), collapse = ", "),
    "\n"
  )
  n_na_cells <- results |>
    dplyr::group_by(method, term_type, cell_id) |>
    dplyr::summarise(all_na = all(is.na(coverage)), .groups = "drop") |>
    dplyr::filter(all_na) |>
    nrow()
  cat("Fully-NA (method, term, cell) combos:", n_na_cells, "\n")
}

# Main Entry Point ------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args) >= 1) args[1] else "smoke"
  mode <- match.arg(mode, c("smoke", "pilot", "full"))

  n_rep <- switch(mode, smoke = 2L, pilot = 10L, full = 50L)
  cell_limit <- if (mode == "smoke") 1L else NULL

  # Remove any stale smoke files before re-running (since we rewrote the DGP)
  if (mode == "smoke") {
    stale <- list.files(
      STUDY4_OUTPUT_DIR,
      pattern = "^cell001_rep.*\\.rds$",
      full.names = TRUE
    )
    if (length(stale) > 0) {
      file.remove(stale)
      cat(sprintf(
        "  Removed %d stale file(s) from %s\n",
        length(stale),
        STUDY4_OUTPUT_DIR
      ))
    }
  }

  cat(sprintf(
    "\nE4 Longitudinal: mode=%s, n_rep=%d, cells=%s\n\n",
    mode,
    n_rep,
    if (is.null(cell_limit)) "all" else paste(cell_limit, collapse = ",")
  ))

  results <- run_study4(n_rep = n_rep, cell_ids_limit = cell_limit)

  if (!is.null(results) && nrow(results) > 0) {
    cat("\n========== E4 RESULTS PREVIEW ==========\n")
    print(
      results |>
        dplyr::select(dplyr::any_of(c(
          "method",
          "term_type",
          "coverage",
          "mean_width",
          "coverage_joint",
          "mean_width_joint",
          "bias",
          "N",
          "J",
          "rep_id",
          "fit_time",
          "converged"
        ))),
      n = 40
    )

    if (n_rep >= 2L) {
      cat("\n--- Summary across reps ---\n")
      print(summarize_study4(results), n = 40)
    }

    validate_study4(results)

    cat("\nFiles written to:", STUDY4_OUTPUT_DIR, "\n")
    cat(
      "Files:",
      length(list.files(STUDY4_OUTPUT_DIR, pattern = "\\.rds$")),
      "\n"
    )
  } else {
    cat("No results returned.\n")
  }

  cat("\nE4 longitudinal study complete.\n")
}
