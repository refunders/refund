# ===========================================================================
# Sandwich Benchmark Part 4: gaulss (Gaussian location-scale) family
# ===========================================================================
#
# DGP: Function-on-scalar, Gaussian location-scale
#   Y_i(t) = beta_0(t) + x_i * beta_1(t) + sigma(t) * eps_i(t)
#   sigma(t) = exp(gamma_0 + gamma_1 * t)   [hetero=TRUE]
#            = exp(gamma_0)                  [hetero=FALSE]
#   eps_i(t) ~ MVN(0, R(rho))               [AR(1)]
#
# Settings (small -- gaulss is slow):
#   n in {50, 100}, rho in {0, 0.6}, hetero in {FALSE, TRUE}
#   B_rep = 50, T_grid = 30
#   4 settings x 50 reps = 200 fits
#
# Methods: nosandwich, hc, cluster_freq, cluster_bayes
#
# Usage:
#   Rscript ci-benchmark/sandwich-benchmark-gaulss.R
# ===========================================================================

devtools::load_all(".")
library(mgcv)
library(dplyr)
library(tidyr)

source("ci-benchmark/benchmark-utils.R")

# ---------------------------------------------------------------------------
# Shared helpers (same pattern as sandwich-benchmark-2.R)
# ---------------------------------------------------------------------------

#' Fit pffr gaulss and extract coefficient CIs for all 4 sandwich methods
#'
#' @param m A fitted pffr model (gaulss family).
#' @param term_name Grep pattern for coef.pffr smterms name.
#' @param truth_vec Numeric vector of true values on the eval grid.
#' @returns Data frame with method, idx, est, se, truth, covered, width.
extract_all_methods <- function(m, term_name, truth_vec) {
  methods <- list(
    nosandwich = list(sandwich = FALSE, freq = FALSE),
    hc = list(sandwich = "hc", freq = FALSE),
    cluster_freq = list(sandwich = TRUE, freq = TRUE),
    cluster_bayes = list(sandwich = TRUE, freq = FALSE)
  )

  results <- list()
  for (mname in names(methods)) {
    args <- methods[[mname]]
    cc <- tryCatch(
      coef(
        m,
        sandwich = args$sandwich,
        freq = args$freq,
        seWithMean = FALSE,
        n1 = 50,
        n2 = 25
      ),
      error = function(e) NULL
    )
    if (is.null(cc)) next

    sm_names <- names(cc$smterms)
    idx <- grep(term_name, sm_names, ignore.case = TRUE)
    if (length(idx) == 0) next
    ti <- cc$smterms[[idx[1]]]
    if (is.null(ti) || is.null(ti$coef)) next

    est <- ti$coef$value
    se <- ti$coef$se
    if (is.null(se) || length(est) != length(truth_vec)) next

    results[[mname]] <- data.frame(
      method = mname,
      idx = seq_along(est),
      est = est,
      se = se,
      truth = truth_vec,
      lower = est - 1.96 * se,
      upper = est + 1.96 * se,
      stringsAsFactors = FALSE
    ) |>
      mutate(
        covered = truth >= lower & truth <= upper,
        width = upper - lower
      )
  }
  bind_rows(results)
}


# ---------------------------------------------------------------------------
# DGP
# ---------------------------------------------------------------------------

generate_gaulss_data <- function(n, T_grid, rho, hetero, seed) {
  set.seed(seed)

  t_grid <- seq(0, 1, length.out = T_grid)

  # Covariates
  x <- rnorm(n)

  # True coefficient functions
  beta0 <- sin(2 * pi * t_grid)
  beta1 <- cos(2 * pi * t_grid)

  # True heteroskedasticity
  gamma0 <- log(0.5)
  gamma1 <- if (hetero) 1.5 else 0
  log_sigma_t <- gamma0 + gamma1 * t_grid
  sigma_t <- exp(log_sigma_t)

  # Mean: n x T
  mu_mat <- outer(rep(1, n), beta0) + outer(x, beta1)

  # Correlated errors
  if (rho > 0) {
    dist_mat <- abs(outer(t_grid, t_grid, "-"))
    R <- rho^(dist_mat / min(diff(t_grid)))
    L <- chol(R)
    eps <- matrix(rnorm(n * T_grid), n, T_grid) %*% L
  } else {
    eps <- matrix(rnorm(n * T_grid), n, T_grid)
  }

  # Scale by sigma(t)
  eps <- sweep(eps, 2, sigma_t, "*")

  Y <- mu_mat + eps

  dat <- data.frame(x = x)
  dat$Y <- I(Y)

  list(
    data = dat,
    t_grid = t_grid,
    beta0 = beta0,
    beta1 = beta1,
    sigma_t = sigma_t
  )
}


# ---------------------------------------------------------------------------
# Scenario grid
# ---------------------------------------------------------------------------

settings <- expand.grid(
  n = c(50L, 100L),
  rho = c(0, 0.6),
  hetero = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

B_rep <- 50L
T_grid <- 30L
base_seed <- 7700L

cat("gaulss Sandwich Benchmark\n")
cat("=========================\n")
cat("Settings:", nrow(settings), "\n")
cat("Reps per setting:", B_rep, "\n")
cat("T_grid:", T_grid, "\n\n")

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

all_results <- list()
counter <- 0L

for (si in seq_len(nrow(settings))) {
  cfg <- settings[si, ]
  cat(sprintf(
    "\n--- Setting %d/%d: n=%d, rho=%.1f, hetero=%s ---\n",
    si,
    nrow(settings),
    cfg$n,
    cfg$rho,
    cfg$hetero
  ))

  for (rep in seq_len(B_rep)) {
    counter <- counter + 1L
    seed <- base_seed + 1000L * si + rep

    sim <- tryCatch(
      generate_gaulss_data(
        n = cfg$n,
        T_grid = T_grid,
        rho = cfg$rho,
        hetero = cfg$hetero,
        seed = seed
      ),
      error = function(e) {
        cat(sprintf("  [%d] DGP error: %s\n", rep, e$message))
        NULL
      }
    )
    if (is.null(sim)) next

    # Fit gaulss model
    m <- tryCatch(
      pffr(
        Y ~ x,
        yind = sim$t_grid,
        data = sim$data,
        family = gaulss()
      ),
      error = function(e) {
        cat(sprintf("  [%d] Fit error: %s\n", rep, e$message))
        NULL
      }
    )
    if (is.null(m)) next

    # Get truth on coef grid
    cc_ref <- tryCatch(
      coef(m, sandwich = FALSE, seWithMean = FALSE, n1 = 50),
      error = function(e) NULL
    )
    if (is.null(cc_ref)) next

    # Find the x term (beta1)
    sm_names <- names(cc_ref$smterms)
    x_idx <- grep(
      "x\\(yindex\\)|yindex\\.vec.*:x",
      sm_names,
      ignore.case = TRUE
    )
    if (length(x_idx) == 0) {
      cat(sprintf(
        "  [%d] Cannot find x term in: %s\n",
        rep,
        paste(sm_names, collapse = ", ")
      ))
      next
    }

    ti <- cc_ref$smterms[[x_idx[1]]]
    t_eval <- ti$x
    truth_on_grid <- cos(2 * pi * t_eval)

    # Extract all methods
    res <- tryCatch(
      extract_all_methods(m, "x\\(yindex\\)|yindex\\.vec.*:x", truth_on_grid),
      error = function(e) {
        cat(sprintf("  [%d] Extract error: %s\n", rep, e$message))
        NULL
      }
    )
    if (is.null(res) || nrow(res) == 0) next

    res$setting_id <- si
    res$rep <- rep
    res$n <- cfg$n
    res$rho <- cfg$rho
    res$hetero <- cfg$hetero

    all_results[[counter]] <- res

    if (rep %% 10 == 0) {
      cat(sprintf("  [%d/%d] done\n", rep, B_rep))
    }
  }
}

results <- bind_rows(all_results)

# ---------------------------------------------------------------------------
# Summarize
# ---------------------------------------------------------------------------

summary_cov <- results |>
  group_by(n, rho, hetero, method) |>
  summarize(
    coverage = mean(covered, na.rm = TRUE),
    mean_width = mean(width, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE),
    n_rep = n_distinct(rep),
    .groups = "drop"
  ) |>
  arrange(n, rho, hetero, method)

cat("\n\n========== COVERAGE SUMMARY ==========\n")
print(summary_cov, n = 100)

# SE calibration
cal <- results |>
  group_by(n, rho, hetero, method, rep) |>
  summarize(mean_se = mean(se), .groups = "drop") |>
  group_by(n, rho, hetero, method) |>
  summarize(
    mean_se = mean(mean_se),
    emp_sd_est = sd(mean_se),
    .groups = "drop"
  )

cat("\n\n========== SE CALIBRATION ==========\n")
print(cal, n = 100)

# Save results
dir.create("ci-benchmark/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, "ci-benchmark/results/gaulss_results.rds")
saveRDS(summary_cov, "ci-benchmark/results/gaulss_summary.rds")
cat("\nResults saved to ci-benchmark/results/gaulss_*.rds\n")
