# explore-dgp.R
#
# Interactive DGP/Model Explorer for pffr Benchmark
#
# Purpose: Run single replicates of benchmark settings, visualize data and all
# terms, show estimate + CI vs truth for every term type, and diagnose term
# contributions and collinearity.
#
# Usage:
#   source("ci-benchmark/explore-dgp.R")
#   exp <- explore_setting(snr = 15, wiggliness = 0.1, seed = 42)
#   exp <- explore_setting(snr = 3, wiggliness = 1, seed = 42, corr_type = "ar1")
#   plot_response(exp)
#   plot_all_terms(exp, method = "pffr")
#   summarize_term_contributions(exp)
#   summarize_coverage(exp, method = "pffr")
#   summarize_coverage(exp, method = "pffr_ar")
#   summarize_coverage(exp, method = "pffr_sandwich")
# Setup -----------------------------------------------------------------------

library(tidyverse)
library(patchwork)

if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else if (file.exists("../DESCRIPTION")) {
  devtools::load_all("..", quiet = TRUE)
} else {
  library(refund)
}

source("ci-benchmark/benchmark-utils.R")
source("ci-benchmark/confint-benchmark.R")

theme_set(theme_minimal(base_size = 12))

# Main Explorer Function ------------------------------------------------------

#' Explore a single benchmark grid row
#'
#' This is the lowest-level entry point for interactive exploration. It runs the
#' exact same data generation, fitting, and metric computation code as the main
#' benchmark (via `run_one_setting()`), while returning the fitted objects and
#' simulated truth for plotting/diagnostics.
#'
#' @param row One-row tibble/list of DGP settings. If `seed` is missing, it is
#'   created from `seed_base + 1000*dgp_id + rep_id`.
#' @param methods Optional methods subset. When provided, it is intersected with
#'   `get_applicable_methods()` for the DGP.
#' @param alpha Significance level for CIs.
#' @param rep_id Replicate id used if `row` has no `rep_id` (default 1).
#' @param seed Optional seed override.
#' @param seed_base Base seed used when generating a seed (default 2026).
#' @returns `pffr_explorer` object.
explore_row <- function(
  row,
  methods = NULL,
  alpha = 0.10,
  rep_id = 1L,
  seed = NULL,
  seed_base = 2026L
) {
  if (inherits(row, "data.frame")) {
    if (nrow(row) != 1) stop("`row` must be a single row.", call. = FALSE)
    row <- as_tibble(row)
  } else {
    stop("`row` must be a one-row tibble/data.frame.", call. = FALSE)
  }

  if (!"dgp_id" %in% names(row)) row$dgp_id <- 1L
  if (!"rep_id" %in% names(row)) row$rep_id <- as.integer(rep_id)

  if (!"seed" %in% names(row) || is.null(row$seed) || is.na(row$seed)) {
    row$seed <- if (!is.null(seed)) {
      seed
    } else {
      seed_base + 1000L * row$dgp_id + row$rep_id
    }
  } else if (!is.null(seed)) {
    row$seed <- seed
  }

  if (!"terms" %in% names(row)) {
    row$terms <- list(c("ff", "linear", "smooth", "concurrent"))
  }

  out <- run_one_setting(row, alpha = alpha, methods = methods)

  if (is.null(out$sim)) {
    dat <- NULL
    truth <- NULL
    s_grid <- NULL
    t_grid <- NULL
    formula <- NULL
  } else {
    dat <- out$sim$data
    truth <- out$sim$truth
    s_grid <- out$sim$s_grid
    t_grid <- out$sim$t_grid
    formula <- build_pffr_formula(
      out$sim,
      s_grid,
      k_smooth = 12,
      k_ff = c(12, 12)
    )
  }

  structure(
    list(
      data = dat,
      truth = truth,
      fits = out$fits,
      timings = out$timings,
      results = out$results,
      s_grid = s_grid,
      t_grid = t_grid,
      terms = row$terms[[1]],
      formula = formula,
      error = out$error,
      settings = list(
        dgp_id = row$dgp_id,
        rep_id = row$rep_id,
        seed = row$seed,
        n = row$n,
        nxgrid = row$nxgrid,
        nygrid = row$nygrid,
        snr = row$snr,
        wiggliness = row$wiggliness,
        corr_type = row$corr_type,
        corr_param = row$corr_param,
        hetero_type = row$hetero_type,
        hetero_param = row$hetero_param,
        error_dist = row$error_dist,
        alpha = alpha,
        methods = methods %||% names(out$fits)
      )
    ),
    class = "pffr_explorer"
  )
}

#' Explore a single benchmark setting
#'
#' Runs a single replicate of a benchmark setting, fitting multiple methods and
#' returning all data needed for visualization and diagnostics.
#'
#' @param snr Signal-to-noise ratio (default: 15)
#' @param wiggliness Wiggliness parameter for random truth generation (default: 0.1)
#' @param corr_type Correlation type: "iid" or "ar1" (default: "iid")
#' @param corr_param Correlation parameter (rho for ar1; period for fourier_pos).
#' @param hetero_type Heteroskedasticity: "none", "linear", "bump" (default: "none")
#' @param hetero_param Heteroskedasticity amplitude parameter.
#' @param error_dist Error distribution: "gaussian" or "t6" (default: "gaussian")
#' @param n Number of observations (default: 100)
#' @param nxgrid Grid size for functional covariate (default: 35)
#' @param nygrid Grid size for response (default: 45)
#' @param terms Character vector of terms to include (default: all)
#' @param methods Character vector of methods to fit. When NULL, uses
#'   `get_applicable_methods()` for the setting.
#' @param seed Random seed (default: NULL for random)
#' @param alpha Significance level for CIs (default: 0.10 for 90% CI).
#' @returns List with data, truth, fits, grids, standardized metrics, and settings
explore_setting <- function(
  snr = 15,
  wiggliness = 0.1,
  corr_type = c("iid", "ar1", "fourier_pos"),
  corr_param = 0.3,
  hetero_type = c("none", "linear", "bump"),
  hetero_param = 1,
  error_dist = c("gaussian", "t6"),
  n = 100,
  nxgrid = 35,
  nygrid = 45,
  terms = c("ff", "linear", "smooth", "concurrent"),
  methods = NULL,
  seed = NULL,
  alpha = 0.10
) {
  corr_type <- match.arg(corr_type)
  hetero_type <- match.arg(hetero_type)
  error_dist <- match.arg(error_dist)

  if (is.null(seed)) seed <- sample.int(1e6, 1)

  cat("Generating data with seed =", seed, "\n")

  applicable <- get_applicable_methods(corr_type, hetero_type, error_dist)
  if (is.null(methods)) {
    methods <- applicable
  } else {
    requested <- unique(methods)
    methods <- intersect(requested, applicable)
    dropped <- setdiff(requested, methods)
    if (length(dropped)) {
      message(
        "Dropping non-applicable method(s): ",
        paste(dropped, collapse = ", ")
      )
    }
  }

  row <- tibble(
    dgp_id = 1L,
    rep_id = 1L,
    seed = seed,
    n = as.integer(n),
    nxgrid = as.integer(nxgrid),
    nygrid = as.integer(nygrid),
    snr = snr,
    wiggliness = wiggliness,
    corr_type = corr_type,
    corr_param = if (corr_type == "iid") NA_real_ else corr_param,
    hetero_type = hetero_type,
    hetero_param = if (hetero_type == "none") NA_real_ else hetero_param,
    error_dist = error_dist,
    terms = list(terms)
  )
  exp <- explore_row(row, methods = methods, alpha = alpha, seed = seed)
  if (is.null(exp$error)) {
    cat("Done. Fitted methods:", paste(names(exp$fits), collapse = ", "), "\n")
  }
  exp
}

#' Explore by `dgp_id` from `make_dgp_settings()`
#'
#' Convenience wrapper for benchmarking-aligned interactive checks.
#'
#' @param dgp_id Integer id from `make_dgp_settings(<size>)`.
#' @param size One of "tiny", "small", "full".
#' @param rep_id Replicate id (used to generate the seed if not overridden).
#' @param seed_base Base seed used to generate `seed = seed_base + 1000*dgp_id + rep_id`.
#' @param methods Optional methods subset.
#' @param alpha Significance level for CIs.
#' @returns `pffr_explorer` object.
explore_dgp_id <- function(
  dgp_id,
  size = c("tiny", "small", "full"),
  rep_id = 1L,
  seed_base = 2026L,
  methods = NULL,
  alpha = 0.10
) {
  size <- match.arg(size)
  dgp <- make_dgp_settings(size)
  row <- dgp |>
    dplyr::filter(.data$dgp_id == .env$dgp_id)
  if (nrow(row) != 1) {
    stop(
      "Could not find unique dgp_id=",
      dgp_id,
      " in settings.",
      call. = FALSE
    )
  }
  row$rep_id <- as.integer(rep_id)
  row$seed <- seed_base + 1000L * row$dgp_id + row$rep_id
  explore_row(row, methods = methods, alpha = alpha)
}

#' Print method for pffr_explorer
print.pffr_explorer <- function(x, ...) {
  cat("pffr Explorer Object\n")
  cat("====================\n")
  if (!is.null(x$error)) {
    cat("ERROR:", x$error, "\n\n")
  }
  cat("Settings:\n")
  cat(
    "  n =",
    x$settings$n,
    ", nxgrid =",
    x$settings$nxgrid,
    ", nygrid =",
    x$settings$nygrid,
    "\n"
  )
  cat("  SNR =", x$settings$snr, ", wiggliness =", x$settings$wiggliness, "\n")
  cat(
    "  Errors:",
    x$settings$corr_type,
    "/",
    x$settings$hetero_type,
    "/",
    x$settings$error_dist,
    "\n"
  )
  if (!is.null(x$settings$corr_param) && is.finite(x$settings$corr_param)) {
    cat("  corr_param:", x$settings$corr_param, "\n")
  }
  if (!is.null(x$settings$hetero_param) && is.finite(x$settings$hetero_param)) {
    cat("  hetero_param:", x$settings$hetero_param, "\n")
  }
  cat("  Seed:", x$settings$seed, "\n")
  cat("Terms:", paste(x$terms, collapse = ", "), "\n")
  cat("Fitted methods:", paste(names(x$fits), collapse = ", "), "\n")
  invisible(x)
}

# Data Visualization Functions ------------------------------------------------

#' Plot response Y(t) as spaghetti or heatmap
#'
#' @param exp Explorer object from explore_setting()
#' @param type "spaghetti" or "heatmap"
#' @param n_curves Number of curves to show in spaghetti plot (default: 30)
#' @returns ggplot object
plot_response <- function(
  exp,
  type = c("spaghetti", "heatmap"),
  n_curves = 30
) {
  type <- match.arg(type)

  Y_mat <- as.matrix(exp$data$Y)
  t_grid <- exp$t_grid
  n <- nrow(Y_mat)

  if (type == "spaghetti") {
    # Sample curves if too many
    idx <- if (n > n_curves) sample(n, n_curves) else seq_len(n)

    plot_df <- expand.grid(t = t_grid, id = idx) |>
      mutate(Y = as.vector(t(Y_mat[idx, ])))

    ggplot(plot_df, aes(x = t, y = Y, group = id)) +
      geom_line(alpha = 0.4, color = "steelblue") +
      labs(
        title = "Response Y(t)",
        subtitle = sprintf("n = %d curves shown (of %d total)", length(idx), n),
        x = "t",
        y = "Y(t)"
      )
  } else {
    # Heatmap
    plot_df <- expand.grid(t = t_grid, obs = seq_len(n)) |>
      mutate(Y = as.vector(t(Y_mat)))

    ggplot(plot_df, aes(x = t, y = obs, fill = Y)) +
      geom_tile() +
      scale_fill_viridis_c() +
      labs(
        title = "Response Y(t) - Heatmap",
        subtitle = sprintf("n = %d observations", n),
        x = "t",
        y = "Observation",
        fill = "Y(t)"
      )
  }
}

#' Plot functional covariate X1(s)
#'
#' @param exp Explorer object
#' @param n_curves Number of curves to show (default: 30)
#' @returns ggplot object
plot_functional_cov <- function(exp, n_curves = 30) {
  X_mat <- as.matrix(exp$data$X1)
  s_grid <- exp$s_grid
  n <- nrow(X_mat)

  idx <- if (n > n_curves) sample(n, n_curves) else seq_len(n)

  plot_df <- expand.grid(s = s_grid, id = idx) |>
    mutate(X = as.vector(t(X_mat[idx, ])))

  ggplot(plot_df, aes(x = s, y = X, group = id)) +
    geom_line(alpha = 0.4, color = "darkorange") +
    labs(
      title = "Functional Covariate X1(s)",
      subtitle = sprintf("n = %d curves shown", length(idx)),
      x = "s",
      y = "X1(s)"
    )
}

#' Plot scalar covariates
#'
#' @param exp Explorer object
#' @returns ggplot object (patchwork)
plot_scalar_covs <- function(exp) {
  dat <- exp$data

  plots <- list()

  if (!is.null(dat$zlin)) {
    plots$zlin <- ggplot(data.frame(z = dat$zlin), aes(x = z)) +
      geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
      labs(title = "zlin (linear covariate)", x = "zlin", y = "Count")
  }

  if (!is.null(dat$zsmoo)) {
    plots$zsmoo <- ggplot(data.frame(z = dat$zsmoo), aes(x = z)) +
      geom_histogram(bins = 30, fill = "darkorange", alpha = 0.7) +
      labs(title = "zsmoo (smooth covariate)", x = "zsmoo", y = "Count")
  }

  if (length(plots) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No scalar covariates")
    )
  }

  wrap_plots(plots, ncol = 2)
}

# Term Coefficient Extraction -------------------------------------------------

get_fit_and_sandwich <- function(exp, method) {
  sandwich_type <- switch(
    method,
    pffr_sandwich = "cluster",
    pffr_hc = "hc",
    FALSE
  )
  fit <- if (method %in% c("pffr_sandwich", "pffr_hc")) {
    exp$fits$pffr
  } else {
    exp$fits[[method]]
  }
  list(fit = fit, sandwich_type = sandwich_type)
}

extract_term_df <- function(exp, term_type, method, alpha = NULL) {
  alpha <- alpha %||% exp$settings$alpha %||% 0.10
  ms <- get_fit_and_sandwich(exp, method)
  if (is.null(ms$fit)) return(NULL)

  extract_term_ci_df(
    ms$fit,
    exp$truth,
    term_type,
    alpha = alpha,
    use_sandwich = ms$sandwich_type,
    s_grid = exp$s_grid,
    t_grid = exp$t_grid,
    data = exp$data
  )
}

# Term Visualization Functions ------------------------------------------------

#' Plot a single term: estimate + CI vs truth
#'
#' @param exp Explorer object
#' @param term Term type: "ff", "linear", "smooth", "concurrent", "intercept"
#' @param method Fitting method (default: "pffr")
#' @returns ggplot object
plot_term <- function(
  exp,
  term = c("ff", "linear", "smooth", "concurrent", "intercept"),
  method = "pffr"
) {
  term <- match.arg(term)

  ms <- get_fit_and_sandwich(exp, method)
  fit <- ms$fit

  if (is.null(fit)) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = paste("Method not fitted:", method)
        )
    )
  }

  coef_df <- extract_term_df(exp, term, method)
  if (is.null(coef_df)) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = paste("Term not found:", term)
        )
    )
  }

  coverage <- mean(coef_df$covered, na.rm = TRUE)

  # Choose plot type based on dimensionality
  if (term == "ff") {
    plot_ff_panels(coef_df, term, method, coverage, exp$settings)
  } else if (term == "smooth") {
    plot_2d_slices(coef_df, term, method, coverage)
  } else {
    plot_1d_ribbon(coef_df, term, method, coverage)
  }
}

#' Plot 1D term with ribbon
plot_1d_ribbon <- function(coef_df, term, method, coverage) {
  ggplot(coef_df, aes(x = x)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "steelblue",
      alpha = 0.3
    ) +
    geom_line(aes(y = estimate), color = "steelblue", linewidth = 1) +
    geom_line(
      aes(y = truth),
      color = "red",
      linewidth = 1,
      linetype = "dashed"
    ) +
    labs(
      title = sprintf("%s term: %s", term, method),
      subtitle = sprintf(
        "Coverage: %.1f%% | Blue = estimate (90%% CI), Red = truth",
        100 * coverage
      ),
      x = "t",
      y = "Coefficient"
    )
}

#' Plot 2D term (smooth) as slices
plot_2d_slices <- function(coef_df, term, method, coverage, n_slices = 6) {
  y_vals <- unique(coef_df$y)
  slice_idx <- round(seq(
    1,
    length(y_vals),
    length.out = min(n_slices, length(y_vals))
  ))

  plot_df <- coef_df |>
    dplyr::filter(.data$y %in% y_vals[slice_idx]) |>
    dplyr::mutate(t_slice = factor(sprintf("t = %.2f", .data$y)))

  ggplot(plot_df, aes(x = x)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "steelblue",
      alpha = 0.3
    ) +
    geom_line(aes(y = estimate), color = "steelblue", linewidth = 0.8) +
    geom_line(
      aes(y = truth),
      color = "red",
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    facet_wrap(~t_slice, scales = "free_y") +
    labs(
      title = sprintf("%s term: %s", term, method),
      subtitle = sprintf(
        "Coverage: %.1f%% | Blue = estimate, Red = truth",
        100 * coverage
      ),
      x = "z",
      y = "f(z, t)"
    )
}

#' Plot ff term as 4-panel heatmaps
plot_ff_panels <- function(coef_df, term, method, coverage, settings) {
  coef_df$error <- coef_df$estimate - coef_df$truth

  # Use identical color scale for estimate and truth
  beta_range <- range(c(coef_df$estimate, coef_df$truth), na.rm = TRUE)

  p1 <- ggplot(coef_df, aes(x = x, y = y, fill = estimate)) +
    geom_tile() +
    scale_fill_viridis_c(limits = beta_range) +
    labs(title = "Fitted", x = "s", y = "t", fill = expression(beta))

  p2 <- ggplot(coef_df, aes(x = x, y = y, fill = truth)) +
    geom_tile() +
    scale_fill_viridis_c(limits = beta_range) +
    labs(title = "Truth", x = "s", y = "t", fill = expression(beta))

  p3 <- ggplot(coef_df, aes(x = x, y = y, fill = error)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(title = "Error", x = "s", y = "t", fill = "Est - Truth")

  p4 <- ggplot(coef_df, aes(x = x, y = y, fill = covered)) +
    geom_tile() +
    scale_fill_manual(values = c("TRUE" = "green3", "FALSE" = "red3")) +
    labs(
      title = sprintf("Coverage: %.1f%%", 100 * coverage),
      x = "s",
      y = "t",
      fill = "Covered"
    )

  (p1 + p2) /
    (p3 + p4) +
    plot_annotation(
      title = sprintf("ff(X1): %s", method),
      subtitle = sprintf(
        "n=%d, SNR=%d, wiggliness=%.2g, seed=%d",
        settings$n,
        settings$snr,
        settings$wiggliness,
        settings$seed
      )
    )
}

#' Plot all terms in a panel layout
#'
#' @param exp Explorer object
#' @param method Fitting method (default: "pffr")
#' @returns ggplot object (patchwork)
plot_all_terms <- function(exp, method = "pffr") {
  term_plots <- list()

  for (term in exp$terms) {
    p <- plot_term(exp, term, method)
    term_plots[[term]] <- p
  }

  # Add intercept if present
  int_plot <- tryCatch(
    plot_term(exp, "intercept", method),
    error = function(e) NULL
  )
  if (!is.null(int_plot) && !inherits(int_plot$layers[[1]], "GeomBlank")) {
    term_plots[["intercept"]] <- int_plot
  }

  if (length(term_plots) == 0) {
    return(
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No terms to plot")
    )
  }

  # Layout depends on number of terms
  wrap_plots(term_plots, ncol = 2) +
    plot_annotation(
      title = sprintf("All Terms: %s", method),
      subtitle = sprintf(
        "SNR=%d, wiggliness=%.2g, %s/%s/%s",
        exp$settings$snr,
        exp$settings$wiggliness,
        exp$settings$corr_type,
        exp$settings$hetero_type,
        exp$settings$error_dist
      )
    )
}

# Diagnostic Summary Functions ------------------------------------------------

#' Summarize term contributions to eta
#'
#' Shows how much variance each term contributes to the linear predictor.
#'
#' @param exp Explorer object
#' @returns Tibble with term variance contributions
summarize_term_contributions <- function(exp) {
  truth <- exp$truth
  eta <- truth$eta
  var_eta <- var(as.vector(eta))

  contributions <- tibble()

  # Prefer truth$etaTerms since these are the *actual* contributions used to
  # construct eta (including centering/identifiability constraints).
  if (!is.null(truth$etaTerms) && length(truth$etaTerms) > 0) {
    term_labels <- names(truth$etaTerms)
    term_labels <- term_labels[term_labels != "epsilon"]

    contributions <- purrr::map_dfr(term_labels, function(lbl) {
      eff <- truth$etaTerms[[lbl]]
      tibble(
        term = lbl,
        var_effect = var(as.vector(eff)),
        mean_abs_effect = mean(abs(eff))
      )
    })
  } else {
    # Fallback for older objects: recompute effects from beta + covariates.
    if ("ff" %in% exp$terms && !is.null(truth$beta[["ff(X1)"]])) {
      X1 <- as.matrix(exp$data$X1)
      beta_ff <- truth$beta[["ff(X1)"]]
      s_grid <- exp$s_grid
      ds <- diff(s_grid[1:2])
      ff_effect <- X1 %*% beta_ff * ds
      contributions <- bind_rows(
        contributions,
        tibble(
          term = "ff(X1)",
          var_effect = var(as.vector(ff_effect)),
          mean_abs_effect = mean(abs(ff_effect))
        )
      )
    }

    if ("linear" %in% exp$terms && !is.null(truth$beta$zlin)) {
      zlin <- exp$data$zlin
      beta_lin <- truth$beta$zlin
      lin_effect <- outer(zlin, beta_lin)
      contributions <- bind_rows(
        contributions,
        tibble(
          term = "zlin",
          var_effect = var(as.vector(lin_effect)),
          mean_abs_effect = mean(abs(lin_effect))
        )
      )
    }

    if ("smooth" %in% exp$terms && !is.null(truth$beta[["s(zsmoo)"]])) {
      zsmoo <- exp$data$zsmoo
      f_zt <- truth$beta[["s(zsmoo)"]]
      if (is.function(f_zt)) {
        smoo_effect <- f_zt(zsmoo, exp$t_grid)
      } else if (is.matrix(f_zt)) {
        smoo_effect <- f_zt
      }
      contributions <- bind_rows(
        contributions,
        tibble(
          term = "s(zsmoo)",
          var_effect = var(as.vector(smoo_effect)),
          mean_abs_effect = mean(abs(smoo_effect))
        )
      )
    }

    if ("concurrent" %in% exp$terms && !is.null(truth$beta$Xconc)) {
      Xconc <- as.matrix(exp$data$Xconc)
      beta_conc <- truth$beta$Xconc
      conc_effect <- sweep(Xconc, 2, beta_conc, "*")
      contributions <- bind_rows(
        contributions,
        tibble(
          term = "Xconc",
          var_effect = var(as.vector(conc_effect)),
          mean_abs_effect = mean(abs(conc_effect))
        )
      )
    }

    if (!is.null(truth$beta$intercept)) {
      int_effect <- matrix(
        truth$beta$intercept,
        nrow = exp$settings$n,
        ncol = exp$settings$nygrid,
        byrow = TRUE
      )
      contributions <- bind_rows(
        contributions,
        tibble(
          term = "intercept",
          var_effect = var(as.vector(int_effect)),
          mean_abs_effect = mean(abs(int_effect))
        )
      )
    }
  }

  if (nrow(contributions) == 0) {
    return(tibble(
      term = character(),
      var_effect = numeric(),
      pct_var_eta = numeric(),
      flag = character()
    ))
  }

  contributions |>
    mutate(
      pct_var_eta = 100 * var_effect / var_eta,
      flag = case_when(
        pct_var_eta < 5 ~ "WEAK (<5%)",
        pct_var_eta > 50 ~ "DOMINANT (>50%)",
        TRUE ~ ""
      )
    ) |>
    arrange(desc(pct_var_eta))
}

#' Summarize coverage and width for all terms
#'
#' @param exp Explorer object
#' @param method Fitting method (default: "pffr")
#' @returns Tibble with coverage metrics
summarize_coverage <- function(exp, method = "pffr") {
  if (is.null(exp$results) || !nrow(exp$results)) {
    return(tibble(term_type = character(), coverage = numeric()))
  }

  res <- exp$results |>
    dplyr::filter(.data$converged, .data$method == .env$method) |>
    dplyr::select(
      term_type,
      coverage,
      coverage_high_var,
      coverage_low_var,
      mean_width,
      rmse,
      bias,
      mean_se,
      z_mean,
      z_sd,
      z_kurtosis
    )

  res |>
    dplyr::mutate(
      coverage_flag = case_when(
        coverage < 0.80 ~ "LOW",
        coverage > 0.95 ~ "HIGH",
        TRUE ~ ""
      )
    )
}

#' Compare two methods side-by-side
#'
#' @param exp Explorer object
#' @param method1 First method (default: "pffr")
#' @param method2 Second method (default: "pffr_sandwich")
#' @returns Tibble with comparison metrics
compare_methods <- function(exp, method1 = "pffr", method2 = "pffr_sandwich") {
  cov1 <- summarize_coverage(exp, method1)
  cov2 <- summarize_coverage(exp, method2)

  if (nrow(cov1) == 0 || nrow(cov2) == 0) {
    return(tibble(term_type = character()))
  }

  inner_join(
    cov1 |>
      select(
        term_type,
        coverage1 = coverage,
        width1 = mean_width,
        se1 = mean_se
      ),
    cov2 |>
      select(
        term_type,
        coverage2 = coverage,
        width2 = mean_width,
        se2 = mean_se
      ),
    by = "term_type"
  ) |>
    mutate(
      width_ratio = width2 / width1,
      se_ratio = se2 / se1,
      coverage_diff = coverage2 - coverage1
    )
}

# mgcv Diagnostic Functions ---------------------------------------------------

#' Run gam.check() diagnostics
#'
#' @param exp Explorer object
#' @param method Fitting method (default: "pffr")
#' @param plot Produce diagnostic plots? (default: TRUE)
#' @returns Invisibly returns the fit; prints diagnostics
run_gam_check <- function(exp, method = "pffr", plot = TRUE) {
  fit <- if (method == "pffr_sandwich") exp$fits$pffr else exp$fits[[method]]

  if (is.null(fit)) {
    message("Method not fitted: ", method)
    return(invisible(NULL))
  }

  cat("=== gam.check() for", method, "===\n\n")

  # Run gam.check (produces plots and prints output)
  # Wrap in tryCatch as mgcv::gam.check can fail with some pffr models
  tryCatch(
    mgcv::gam.check(fit, type = "deviance"),
    error = function(e) {
      message("Note: gam.check() plotting failed: ", e$message)
      message("Running without plots...")
      if (plot) {
        # Try just the text output part
        k_check <- tryCatch(mgcv::k.check(fit), error = function(e) NULL)
        if (!is.null(k_check)) {
          cat("\n=== Basis dimension check ===\n")
          print(k_check)
        }
      }
    }
  )

  cat("\n=== Model Summary ===\n")
  cat("Deviance explained:", round(100 * summary(fit)$dev.expl, 1), "%\n")
  cat("Scale estimate:", round(fit$sig2, 4), "\n")
  cat("GCV score:", round(fit$gcv.ubre, 4), "\n")

  invisible(fit)
}

#' Check concurvity between smooth terms
#'
#' @param exp Explorer object
#' @param method Fitting method (default: "pffr")
#' @param threshold Threshold for flagging high concurvity (default: 0.8)
#' @returns Tibble with concurvity metrics
check_concurvity <- function(exp, method = "pffr", threshold = 0.8) {
  fit <- if (method == "pffr_sandwich") exp$fits$pffr else exp$fits[[method]]

  if (is.null(fit)) {
    message("Method not fitted: ", method)
    return(NULL)
  }

  cat("=== Concurvity Analysis for", method, "===\n\n")

  # Worst case concurvity
  worst <- mgcv::concurvity(fit, full = TRUE)
  cat("Worst-case concurvity:\n")
  print(round(worst, 3))

  # Pairwise concurvity
  cat("\nPairwise concurvity matrix:\n")
  pairwise <- mgcv::concurvity(fit, full = FALSE)

  # Extract the estimate matrix
  if (is.list(pairwise) && "estimate" %in% names(pairwise)) {
    pw_mat <- pairwise$estimate
    print(round(pw_mat, 3))

    # Flag high concurvity pairs
    high_conc <- which(pw_mat > threshold & pw_mat < 1, arr.ind = TRUE)
    if (nrow(high_conc) > 0) {
      cat("\n*** HIGH CONCURVITY PAIRS (>", threshold, "):\n")
      for (i in seq_len(nrow(high_conc))) {
        row_nm <- rownames(pw_mat)[high_conc[i, 1]]
        col_nm <- colnames(pw_mat)[high_conc[i, 2]]
        val <- pw_mat[high_conc[i, 1], high_conc[i, 2]]
        cat(sprintf("  %s ~ %s: %.3f\n", row_nm, col_nm, val))
      }
    }
  }

  # Return worst case as tibble
  tibble(
    smooth = names(worst["worst", ]),
    concurvity_worst = as.numeric(worst["worst", ]),
    flag = ifelse(concurvity_worst > threshold, "HIGH", "")
  )
}

#' Check collinearity of term effects
#'
#' Computes correlation matrix of term contributions to eta.
#'
#' @param exp Explorer object
#' @returns Correlation matrix
check_collinearity <- function(exp) {
  truth <- exp$truth

  effects <- list()

  # Compute effect vectors for each term
  if ("ff" %in% exp$terms && !is.null(truth$beta[["ff(X1)"]])) {
    X1 <- as.matrix(exp$data$X1)
    beta_ff <- truth$beta[["ff(X1)"]]
    ds <- diff(exp$s_grid[1:2])
    effects$ff <- as.vector(X1 %*% beta_ff * ds)
  }

  if ("linear" %in% exp$terms && !is.null(truth$beta$zlin)) {
    effects$linear <- as.vector(outer(exp$data$zlin, truth$beta$zlin))
  }

  if ("smooth" %in% exp$terms && !is.null(truth$beta[["s(zsmoo)"]])) {
    f_zt <- truth$beta[["s(zsmoo)"]]
    if (is.function(f_zt)) {
      effects$smooth <- as.vector(f_zt(exp$data$zsmoo, exp$t_grid))
    }
  }

  if ("concurrent" %in% exp$terms && !is.null(truth$beta$Xconc)) {
    Xconc <- as.matrix(exp$data$Xconc)
    effects$concurrent <- as.vector(sweep(Xconc, 2, truth$beta$Xconc, "*"))
  }

  if (length(effects) < 2) {
    cat("Need at least 2 terms to compute collinearity\n")
    return(NULL)
  }

  # Build matrix and compute correlation
  effect_mat <- do.call(cbind, effects)
  cor_mat <- cor(effect_mat)

  cat("=== Effect Collinearity Matrix ===\n")
  print(round(cor_mat, 3))

  # Flag high correlations
  high_cor <- which(abs(cor_mat) > 0.7 & cor_mat < 1, arr.ind = TRUE)
  if (nrow(high_cor) > 0) {
    cat("\n*** HIGH CORRELATIONS (|r| > 0.7):\n")
    for (i in seq_len(nrow(high_cor))) {
      if (high_cor[i, 1] < high_cor[i, 2]) {
        row_nm <- rownames(cor_mat)[high_cor[i, 1]]
        col_nm <- colnames(cor_mat)[high_cor[i, 2]]
        val <- cor_mat[high_cor[i, 1], high_cor[i, 2]]
        cat(sprintf("  %s ~ %s: %.3f\n", row_nm, col_nm, val))
      }
    }
  }

  invisible(cor_mat)
}

# Example Usage (commented out for sourcing) ----------------------------------

if (FALSE) {
  # === Interactive Workflow Example ===
  dgps <- make_dgp_settings(size = "full")
  # Run single setting
  exp <- explore_row(dgps[1, ]) # simple iid, low wiggle
  exp <- explore_row(dgps[7, ]) # simple iid, high wiggle
  exp <- explore_row(dgps[11, ]) # ar1 .9, high wiggle
  exp <- explore_row(dgps[60, ]) # het large,high wiggle
  print(exp)

  # What do the data look like?
  plot_response(exp)
  plot_response(exp, type = "heatmap")
  plot_functional_cov(exp)
  plot_scalar_covs(exp)

  # How much does each term contribute?
  summarize_term_contributions(exp)
  # -> Check if ff is weak relative to others

  # Visualize all term fits
  plot_all_terms(exp, method = "pffr")
  plot_all_terms(exp, method = "pffr_ar")
  plot_all_terms(exp, method = "pffr_sandwich")

  # Drill into problematic term
  plot_term(exp, term = "ff", method = "pffr")
  plot_term(exp, term = "linear", method = "pffr")

  # Coverage summary
  summarize_coverage(exp, method = "pffr")
  summarize_coverage(exp, method = "pffr_sandwich")

  # Compare pffr vs sandwich
  compare_methods(exp, "pffr", "pffr_sandwich")

  # mgcv diagnostics
  run_gam_check(exp, method = "pffr")
  check_concurvity(exp, method = "pffr")

  # Check collinearity of effects
  check_collinearity(exp)

  # === Test with different settings ===

  # High wiggliness (harder to fit)
  exp_wiggly <- explore_setting(snr = 15, wiggliness = 10, seed = 42)
  summarize_coverage(exp_wiggly, "pffr")

  # AR1 errors
  exp_ar <- explore_setting(snr = 15, corr_type = "ar1", seed = 42)
  compare_methods(exp_ar, "pffr", "pffr_ar")

  # Heteroskedastic errors
  exp_het <- explore_setting(snr = 15, hetero_type = "linear", seed = 42)
}
