#' ---
#' title: "pffr CI Benchmark: Gaussian (Round 2) + Non-Gaussian (Study 1) + Grid Refinement (Study 2)"
#' author: "Benchmark Analysis"
#' date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float:
#'       collapsed: false
#'     toc_depth: 3
#'     code_folding: hide
#'     fig_width: 10
#'     fig_height: 6
#'     theme: flatly
#' ---

#+ setup, include = FALSE
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 10,
  fig.height = 6
)
# knitr defaults to the Rmd's directory (ci-benchmark/); reset to repo root
# so that source() calls in sourced scripts (e.g., confint-benchmark.R) work.
if (basename(getwd()) == "ci-benchmark") {
  knitr::opts_knit$set(root.dir = normalizePath(".."))
}

#+ load-packages
library(tidyverse)
library(gt)
library(patchwork)

theme_set(theme_minimal(base_size = 11))

# Path helper: works from repo root (Rscript) or ci-benchmark/ (knitr WD)
ci_path <- function(path) {
  if (file.exists(path) || dir.exists(path)) return(path)
  # Strip ci-benchmark/ prefix and try again (knitr sets WD to ci-benchmark/)
  alt <- sub("^ci-benchmark/", "", path)
  if (file.exists(alt) || dir.exists(alt)) return(alt)
  path # return original (will fail with informative error downstream)
}

# Constants -------------------------------------------------------------------

NOMINAL_COVERAGE <- 0.90
COLORS_METHOD <- c(
  pffr = "#1b9e77",
  pffr_sandwich = "#d95f02",
  pffr_hc = "#e6ab02",
  pffr_gaulss = "#7570b3",
  pffr_ar = "#66a61e"
)

#'
#' # Theory and Study Design
#'
#' This report combines the Gaussian Round 2 benchmark with two follow-up
#' studies (non-Gaussian families and grid refinement). The main Gaussian
#' benchmark evaluates pointwise confidence intervals for functional effects in
#' a function-on-function additive model.
#'
#' Current report components:
#'
#' - **Gaussian Round 2 (main benchmark):** production Gaussian / heavy-tailed
#'   error runs loaded from the benchmark results directory, with
#'   heteroskedasticity and multiple within-curve correlation structures.
#' - **Study 1 (non-Gaussian sandwich coverage):** Poisson and Binomial DGPs for
#'   `ff(X1) + z_lin`, with IID / AR1(0.9) / `fourier_pos(0.3)`,
#'   `n = 200, 400`, and methods `default`, `hc`, `cluster`.
#' - **Study 2 (grid refinement):** Gaussian design with wiggliness = 5,
#'   `n = 20, 40, 80`, `snr = 3, 25`, IID / AR1(0.9) /
#'   `fourier_pos(0.3)`, `nygrid = 40, 80, 120` at fixed
#'   `nxgrid = 60`, targeting 120 reps per DGP \u00d7 grid cell (with
#'   completeness checked below).
#'
#' For the main Gaussian benchmark, data are generated from
#' $Y_i(t) = \eta_i(t) + \epsilon_i(t)$ where $\eta_i(t)$ combines:
#' a functional intercept, a function-on-function term
#' $\int X_i(s)\beta(s,t)ds$, a varying-coefficient linear term, a smooth
#' scalar-varying term, and a concurrent term.
#'
#' DGP difficulty in the main Gaussian Round 2 benchmark is controlled by:
#'
#' - sample size (`n`),
#' - signal-to-noise ratio (`snr`),
#' - truth wiggliness (`wiggliness`),
#' - error distribution (`gaussian`, `t6`),
#' - residual covariance structure: IID, AR(1), or periodic non-monotone
#'   (`fourier_pos`),
#' - heteroskedasticity pattern over $t$.
#'
#' Coverage is assessed pointwise on evaluation grids for each term type.
#' Study 1 and Study 2 use study-specific summaries and diagnostics defined in
#' their sections below.
#'
#' CI estimators compared in this report:
#'
#' - **Gaussian Round 2:** `pffr`, `pffr_hc`, `pffr_sandwich`, `pffr_ar`,
#'   `pffr_gaulss`.
#' - **Study 1:** `default`, `hc`, `cluster`.
#' - **Study 2:** `default`, `cluster`, with optional `hc` / `cl2` in newer
#'   result extracts.
#'
#' Gaussian Round 2 estimator labels:
#'
#' - `pffr`: default mgcv/refund Bayesian covariance-based pointwise CIs.
#' - `pffr_hc`: observation-level HC sandwich covariance (heteroskedasticity
#'   robust, not cluster/correlation robust).
#' - `pffr_sandwich`: cluster-robust sandwich covariance (clusters = curves),
#'   targeting within-curve dependence and heteroskedasticity.
#' - `pffr_ar`: `bam(..., rho=...)` AR(1) working-correlation fit, with
#'   CIs from the fitted AR model covariance.
#' - `pffr_gaulss`: Gaussian location-scale fit, modeling variance over $t$,
#'   with pointwise CIs from the fitted gaulss covariance.
#'
# Helper functions ------------------------------------------------------------

#' Format coverage with MC SE
fmt_coverage <- function(cov, se) {
  sprintf("%.1f%% (\u00b1%.1f)", cov * 100, 1.96 * se * 100)
}

#' Summarize coverage by grouping variables
#'
#' Note: MC SE is sd(coverage)/sqrt(n), NOT sqrt(p(1-p)/n).
#' Each coverage observation is already a proportion over many grid points,
#' not a single Bernoulli trial.
summarize_coverage <- function(data, ...) {
  data |>
    group_by(...) |>
    summarize(
      mean_coverage = mean(coverage, na.rm = TRUE),
      mc_se = sd(coverage, na.rm = TRUE) / sqrt(n()),
      pct_nominal = mean(coverage >= NOMINAL_COVERAGE, na.rm = TRUE) * 100,
      mean_width = mean(mean_width, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      n_obs = n(),
      .groups = "drop"
    )
}

#' Create effect plot comparing levels of a factor
plot_parameter_effect <- function(
  data,
  param,
  param_label,
  facet_by = "term_type"
) {
  p <- ggplot(
    data,
    aes(
      x = .data[[param]],
      y = mean_coverage,
      color = method,
      group = method
    )
  ) +
    geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
    geom_errorbar(
      aes(
        ymin = mean_coverage - 1.96 * mc_se,
        ymax = mean_coverage + 1.96 * mc_se
      ),
      width = 0.15,
      linewidth = 0.5,
      position = position_dodge(width = 0.3)
    ) +
    geom_hline(
      yintercept = NOMINAL_COVERAGE,
      linetype = "dashed",
      color = "red",
      linewidth = 0.7
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0.5, 1)) +
    scale_color_manual(values = COLORS_METHOD) +
    labs(x = param_label, y = "Coverage", color = "Method") +
    theme(legend.position = "bottom")

  if (!is.null(facet_by)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)), ncol = 3)
  }
  p
}

#' Create effect table for a parameter
make_effect_table <- function(data, param, param_label) {
  data |>
    mutate(coverage_fmt = fmt_coverage(mean_coverage, mc_se)) |>
    select(method, term_type, !!sym(param), coverage_fmt) |>
    pivot_wider(names_from = all_of(param), values_from = coverage_fmt) |>
    gt() |>
    tab_header(
      title = paste("Coverage by", param_label),
      subtitle = "Mean % (\u00b195% MC interval)"
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    )
}

#'
#' # Data Loading
#'
#' Load results from individual RDS files.

#+ load-data
# Resolve results directory:
# - If REFUND_BENCHMARK_RESULTS_DIR is set, use it.
# - Otherwise prefer production results over pilot outputs.
resolve_results_dir <- function() {
  env_dir <- Sys.getenv("REFUND_BENCHMARK_RESULTS_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(env_dir)
  }

  candidates <- c(
    "ci-benchmark/results",
    "results",
    "ci-benchmark/pilot-results",
    "pilot-results"
  )

  has_files <- vapply(
    candidates,
    function(dir) {
      dir.exists(dir) &&
        length(list.files(dir, pattern = "^dgp.*\\.rds$")) > 0
    },
    logical(1)
  )

  if (any(has_files)) {
    candidates[which(has_files)[1]]
  } else {
    candidates[1]
  }
}

results_dir <- resolve_results_dir()
files <- list.files(results_dir, pattern = "^dgp.*\\.rds$", full.names = TRUE)

d_raw <- map_dfr(files, readRDS)

cat("Results directory:", results_dir, "\n")
cat("Loaded", length(files), "result files\n")
cat("Total observations:", nrow(d_raw), "\n")
cat("Unique DGPs:", n_distinct(d_raw$dgp_id), "\n")
cat(
  "Replications per DGP (range):",
  paste(
    range(table(paste(d_raw$dgp_id, d_raw$term_type, d_raw$method))),
    collapse = "-"
  ),
  "\n"
)

# Guardrail: fail fast if dgp_id definitions are mixed (stale files from
# different benchmark designs in the same output directory).
allow_mixed <- identical(Sys.getenv("REFUND_ALLOW_MIXED_DGP", "0"), "1")
dgp_defs <- d_raw |>
  select(
    dgp_id,
    n,
    nxgrid,
    nygrid,
    snr,
    wiggliness,
    corr_type,
    corr_param,
    hetero_type,
    hetero_param,
    error_dist
  ) |>
  distinct()

mixed_dgp <- dgp_defs |>
  count(dgp_id, name = "n_definitions") |>
  filter(n_definitions > 1)

if (nrow(mixed_dgp) > 0) {
  msg <- paste0(
    "Detected mixed DGP definitions for ",
    nrow(mixed_dgp),
    " dgp_id(s): ",
    paste(head(mixed_dgp$dgp_id, 10), collapse = ", "),
    if (nrow(mixed_dgp) > 10) ", ..." else "",
    ". This usually means stale results from prior designs are mixed in the same ",
    "directory. Clean/recompute results or set REFUND_ALLOW_MIXED_DGP=1 to ",
    "override this check."
  )
  if (allow_mixed) {
    warning(msg, call. = FALSE)
  } else {
    stop(msg, call. = FALSE)
  }
}

#'
#' ## Parameter Grid Summary

#+ param-grid
params <- d_raw |>
  select(
    dgp_id,
    n,
    snr,
    wiggliness,
    corr_type,
    corr_param,
    hetero_type,
    hetero_param,
    error_dist
  ) |>
  distinct()

tribble(
  ~Parameter,
  ~Values,
  "n",
  paste(sort(unique(d_raw$n)), collapse = ", "),
  "SNR",
  paste(sort(unique(d_raw$snr)), collapse = ", "),
  "Wiggliness",
  paste(sort(unique(d_raw$wiggliness)), collapse = ", "),
  "Correlation type",
  paste(sort(unique(d_raw$corr_type)), collapse = ", "),
  "Correlation param (rho)",
  paste(
    sort(unique(d_raw$corr_param[!is.na(d_raw$corr_param)])),
    collapse = ", "
  ),
  "Heteroskedasticity",
  paste(sort(unique(d_raw$hetero_type)), collapse = ", "),
  "Hetero param",
  paste(
    sort(unique(d_raw$hetero_param[!is.na(d_raw$hetero_param)])),
    collapse = ", "
  ),
  "Error distribution",
  paste(sort(unique(d_raw$error_dist)), collapse = ", "),
  "Methods",
  paste(sort(unique(d_raw$method)), collapse = ", "),
  "Term types",
  paste(sort(unique(d_raw$term_type)), collapse = ", ")
) |>
  gt() |>
  tab_header(title = "Benchmark Design (Round 2)")

#'
#' ## Data Preparation

#+ prepare-data
# Derive rho from corr_param (rho=0 maps to iid, rho>0 to ar1)
d <- d_raw |>
  mutate(
    rho = case_when(
      corr_type == "iid" ~ 0,
      corr_type == "ar1" ~ corr_param,
      corr_type == "fourier_pos" ~ NA_real_,
      TRUE ~ corr_param
    ),
    method = factor(method, levels = names(COLORS_METHOD)),
    term_type = factor(
      term_type,
      levels = c(
        "intercept",
        "E(Y)",
        "linear",
        "concurrent",
        "smooth",
        "ff"
      )
    ),
    n_f = factor(n),
    snr_f = factor(snr, levels = c(5, 25), labels = c("Low (5)", "High (25)")),
    wiggliness_f = factor(
      wiggliness,
      levels = c(0.01, 5),
      labels = c("Smooth (0.01)", "Wiggly (5)")
    ),
    rho_f = factor(
      ifelse(corr_type == "fourier_pos", "fourier", as.character(rho)),
      levels = c("0", "0.3", "0.9", "fourier")
    ),
    corr_f = factor(
      corr_type,
      levels = c("iid", "ar1", "fourier_pos"),
      labels = c("IID", "AR(1)", "Fourier+")
    ),
    hetero_f = factor(
      ifelse(hetero_type == "none", "None", paste0("Bump(", hetero_param, ")")),
      levels = c("None", "Bump(1)", "Bump(3)")
    ),
    error_f = factor(
      error_dist,
      levels = c("gaussian", "t6"),
      labels = c("Gaussian", "t(6)")
    )
  )

#'
#' # Overview: Global Performance
#'
#' Before diving into parameter effects, here's the big picture.

#+ global-summary
global_summary <- d |>
  group_by(method, term_type) |>
  summarize(
    mean_coverage = mean(coverage, na.rm = TRUE),
    mc_se = sd(coverage, na.rm = TRUE) / sqrt(n()),
    pct_nominal = mean(coverage >= NOMINAL_COVERAGE, na.rm = TRUE) * 100,
    mean_width = mean(mean_width, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    .groups = "drop"
  )

#+ global-heatmap, fig.height = 5
ggplot(global_summary, aes(x = term_type, y = method, fill = mean_coverage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.0f%%", mean_coverage * 100)),
    color = "white",
    fontface = "bold",
    size = 4
  ) +
  scale_fill_gradient2(
    low = "#d73027",
    mid = "#fee08b",
    high = "#1a9850",
    midpoint = 0.85,
    limits = c(0.6, 1),
    labels = scales::percent
  ) +
  labs(
    title = "Overall Coverage by Method and Term Type",
    subtitle = "Averaged across ALL DGP conditions (includes misspecification)",
    x = "Term Type",
    y = "Method",
    fill = "Coverage"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#' **Interpretation**: This shows *marginal* performance averaging over all
#' conditions. Methods that only apply to subsets (e.g., pffr_ar on AR1 +
#' gaussian, pffr_gaulss on hetero + gaussian) face different DGP mixes.
#'
#' # Z-Score Diagnostics {.tabset .tabset-fade .tabset-pills}
#'
#' Z-scores = (estimate - truth) / SE. Under correct specification,
#' z ~ N(0,1): mean ~ 0, sd ~ 1, excess kurtosis ~ 0.
#'
#' ## Z-Score Summary

#+ zscore-summary
zscore_summary <- d |>
  group_by(method, term_type) |>
  summarize(
    mean_z_mean = mean(z_mean, na.rm = TRUE),
    mean_z_sd = mean(z_sd, na.rm = TRUE),
    mean_z_kurtosis = mean(z_kurtosis, na.rm = TRUE),
    mean_z2 = mean(z2_mean, na.rm = TRUE),
    .groups = "drop"
  )

zscore_summary |>
  gt() |>
  tab_header(
    title = "Z-Score Diagnostics by Method and Term",
    subtitle = "Under N(0,1): mean=0, sd=1, kurtosis=0, E[z^2]=1"
  ) |>
  fmt_number(
    columns = c(mean_z_mean, mean_z_sd, mean_z_kurtosis, mean_z2),
    decimals = 3
  )

#'
#' ## Z-Score SD Heatmap
#'
#' z_sd > 1 means SEs are too small (undercoverage); z_sd < 1 means SEs are
#' too large (overcoverage / conservative).

#+ zscore-sd-heatmap, fig.height = 5
ggplot(zscore_summary, aes(x = term_type, y = method, fill = mean_z_sd)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", mean_z_sd)),
    fontface = "bold",
    size = 3.5
  ) +
  scale_fill_gradient2(
    low = "#3182bd",
    mid = "white",
    high = "#e6550d",
    midpoint = 1,
    limits = c(0.5, 2)
  ) +
  labs(
    title = "Z-Score SD by Method and Term",
    subtitle = "Blue < 1 (conservative), White = 1 (ideal), Orange > 1 (anti-conservative)",
    x = "Term Type",
    y = "Method",
    fill = "z_sd"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#'
#' ## Z-Score SD by Correlation Strength

#+ zscore-by-rho, fig.height = 7
zscore_by_rho <- d |>
  filter(!is.na(rho_f)) |>
  group_by(method, term_type, rho_f) |>
  summarize(
    mean_z_sd = mean(z_sd, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  zscore_by_rho,
  aes(x = rho_f, y = mean_z_sd, color = method, group = method)
) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  facet_wrap(~term_type, ncol = 3) +
  scale_color_manual(values = COLORS_METHOD) +
  labs(
    title = "Z-Score SD by Correlation Strength",
    subtitle = "z_sd > 1 = SEs too small (undercoverage)",
    x = "Correlation (rho)",
    y = "Mean z_sd",
    color = "Method"
  ) +
  theme(legend.position = "bottom")

#'
#' # SE Calibration {.tabset .tabset-fade .tabset-pills}
#'
#' SE calibration via z-score SD: z_sd = sd((est - truth) / se) at each grid
#' point, averaged over the grid within each rep. Should be ~1 for
#' well-calibrated SEs. z_sd > 1 means SEs are too small (anti-conservative);
#' z_sd < 1 means SEs are too large (conservative).
#'
#' ## Z-SD Heatmap

#+ se-cal-compute
se_cal_summary <- d |>
  group_by(method, term_type) |>
  summarize(
    mean_z_sd = mean(z_sd, na.rm = TRUE),
    median_z_sd = median(z_sd, na.rm = TRUE),
    .groups = "drop"
  )

#+ se-cal-heatmap, fig.height = 5
ggplot(se_cal_summary, aes(x = term_type, y = method, fill = mean_z_sd)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", mean_z_sd)),
    fontface = "bold",
    size = 3.5
  ) +
  scale_fill_gradient2(
    low = "#3182bd",
    mid = "white",
    high = "#e6550d",
    midpoint = 1,
    limits = c(0.5, 2)
  ) +
  labs(
    title = "SE Calibration (z-score SD) by Method and Term",
    subtitle = "Blue < 1 (conservative), White = 1 (ideal), Orange > 1 (anti-conservative)",
    x = "Term Type",
    y = "Method",
    fill = "z_sd"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#'
#' ## Z-SD by Correlation

#+ se-cal-by-corr, fig.height = 7
z_sd_by_corr <- d |>
  mutate(
    rho_f = case_when(
      corr_type == "iid" ~ "0",
      corr_type == "fourier_pos" ~ "fourier",
      TRUE ~ as.character(corr_param)
    ),
    rho_f = factor(rho_f, levels = c("0", "0.3", "0.9", "fourier"))
  ) |>
  filter(!is.na(rho_f)) |>
  group_by(method, term_type, rho_f) |>
  summarize(mean_z_sd = mean(z_sd, na.rm = TRUE), .groups = "drop")

ggplot(
  z_sd_by_corr,
  aes(
    x = rho_f,
    y = mean_z_sd,
    color = method,
    group = method
  )
) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  facet_wrap(~term_type, ncol = 3) +
  scale_color_manual(values = COLORS_METHOD) +
  labs(
    title = "SE Calibration (z-score SD) by Correlation Strength",
    subtitle = "z_sd > 1: SEs too small (anti-conservative). z_sd < 1: SEs too large.",
    x = "Correlation",
    y = "z_sd",
    color = "Method"
  ) +
  theme(legend.position = "bottom")

#'
#' # Region-Specific Coverage {.tabset .tabset-fade .tabset-pills}
#'
#' Under heteroskedasticity, coverage may differ between high- and
#' low-variance regions. The cluster-robust sandwich should equalize this.
#'
#' ## Region Coverage Comparison

#+ region-coverage
region_data <- d |>
  filter(
    hetero_type != "none",
    !is.na(coverage_high_var),
    !is.na(coverage_low_var)
  ) |>
  select(
    method,
    term_type,
    dgp_id,
    rep_id,
    hetero_f,
    coverage_high_var,
    coverage_low_var
  ) |>
  pivot_longer(
    cols = c(coverage_high_var, coverage_low_var),
    names_to = "region",
    values_to = "coverage"
  ) |>
  mutate(
    region = ifelse(
      region == "coverage_high_var",
      "High variance",
      "Low variance"
    )
  )

if (nrow(region_data) > 0) {
  region_summary <- region_data |>
    group_by(method, term_type, region) |>
    summarize(
      mean_cov = mean(coverage, na.rm = TRUE),
      mc_se = sd(coverage, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  ggplot(
    region_summary,
    aes(
      x = term_type,
      y = mean_cov,
      fill = region
    )
  ) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
    geom_errorbar(
      aes(
        ymin = mean_cov - 1.96 * mc_se,
        ymax = mean_cov + 1.96 * mc_se
      ),
      position = position_dodge(0.8),
      width = 0.2
    ) +
    geom_hline(
      yintercept = NOMINAL_COVERAGE,
      linetype = "dashed",
      color = "red"
    ) +
    facet_wrap(~method, ncol = 3) +
    scale_y_continuous(labels = scales::percent, limits = c(0.5, 1)) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Coverage by Variance Region (Heteroskedastic DGPs)",
      subtitle = "Robust methods should equalize coverage across regions",
      x = "Term Type",
      y = "Coverage",
      fill = "Region"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
} else {
  cat("No region-specific coverage data available.\n")
}

#'
#' ## Region Coverage Gap

#+ region-gap
if (nrow(region_data) > 0) {
  region_gap <- d |>
    filter(
      hetero_type != "none",
      !is.na(coverage_high_var),
      !is.na(coverage_low_var)
    ) |>
    mutate(coverage_gap = coverage_low_var - coverage_high_var) |>
    group_by(method, term_type, hetero_f) |>
    summarize(
      mean_gap = mean(coverage_gap, na.rm = TRUE) * 100,
      se_gap = sd(coverage_gap, na.rm = TRUE) / sqrt(n()) * 100,
      .groups = "drop"
    )

  ggplot(region_gap, aes(x = term_type, y = mean_gap, fill = method)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = mean_gap - 1.96 * se_gap, ymax = mean_gap + 1.96 * se_gap),
      position = position_dodge(0.8),
      width = 0.2
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~hetero_f) +
    scale_fill_manual(values = COLORS_METHOD) +
    labs(
      title = "Coverage Gap: Low Variance - High Variance Region",
      subtitle = "Positive = better coverage in low-variance region. 0 = uniform.",
      x = "Term Type",
      y = "Coverage gap (pp)",
      fill = "Method"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
}

#'
#' # DGP Parameter Effects {.tabset .tabset-fade .tabset-pills}
#'
#' ## Correlation Strength (rho) {.tabset}
#'
#' The key design variable: rho in {0, 0.3, 0.9} for AR(1).
#'
#' ### Coverage by rho

#+ rho-effect-plot
rho_summary <- d |>
  filter(corr_type %in% c("iid", "ar1")) |>
  mutate(rho_val = ifelse(corr_type == "iid", 0, corr_param)) |>
  summarize_coverage(method, term_type, rho_val) |>
  mutate(rho_val = factor(rho_val))

plot_parameter_effect(rho_summary, "rho_val", "Correlation (rho)") +
  labs(
    title = "Effect of AR(1) Correlation on Coverage",
    subtitle = "Cluster sandwich should maintain coverage; default and HC should degrade"
  )

#'
#' ### rho Effect Table

#+ rho-effect-table
make_effect_table(rho_summary, "rho_val", "Correlation (rho)")

#'
#' ## Signal-to-Noise Ratio (SNR) {.tabset}

#+ snr-effect-plot
snr_summary <- d |> summarize_coverage(method, term_type, snr_f)

plot_parameter_effect(snr_summary, "snr_f", "Signal-to-Noise Ratio") +
  labs(
    title = "Effect of SNR on Coverage",
    subtitle = "Higher SNR generally improves coverage"
  )

#'
#' ### SNR Effect Table

#+ snr-effect-table
make_effect_table(snr_summary, "snr_f", "SNR")

#'
#' ## Wiggliness {.tabset}

#+ wig-effect-plot
wig_summary <- d |> summarize_coverage(method, term_type, wiggliness_f)

plot_parameter_effect(wig_summary, "wiggliness_f", "Wiggliness") +
  labs(
    title = "Effect of Wiggliness on Coverage",
    subtitle = "Low wiggliness can cause smoothing bias"
  )

#'
#' ### Wiggliness Effect Table

#+ wig-effect-table
make_effect_table(wig_summary, "wiggliness_f", "Wiggliness")

#'
#' ## Heteroskedasticity {.tabset}

#+ hetero-effect-plot
hetero_summary <- d |> summarize_coverage(method, term_type, hetero_f)

plot_parameter_effect(hetero_summary, "hetero_f", "Heteroskedasticity") +
  labs(
    title = "Effect of Heteroskedasticity on Coverage",
    subtitle = "Sandwich and GAULSS methods should be robust"
  )

#'
#' ### Heteroskedasticity Effect Table

#+ hetero-effect-table
make_effect_table(hetero_summary, "hetero_f", "Heteroskedasticity")

#'
#' ### Robust Methods Under Heteroskedasticity

#+ hetero-robust-compare, fig.height = 5
d |>
  filter(method %in% c("pffr", "pffr_sandwich", "pffr_hc", "pffr_gaulss")) |>
  summarize_coverage(method, hetero_f) |>
  ggplot(aes(x = hetero_f, y = mean_coverage, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  geom_errorbar(
    aes(
      ymin = mean_coverage - 1.96 * mc_se,
      ymax = mean_coverage + 1.96 * mc_se
    ),
    position = position_dodge(0.8),
    width = 0.2
  ) +
  geom_hline(
    yintercept = NOMINAL_COVERAGE,
    linetype = "dashed",
    color = "red"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0.7, 1)) +
  scale_fill_manual(values = COLORS_METHOD) +
  labs(
    title = "Methods Under Heteroskedasticity",
    subtitle = "pffr_sandwich, pffr_hc, and pffr_gaulss vs standard pffr",
    x = "Heteroskedasticity Type",
    y = "Coverage"
  )

#'
#' ## Error Distribution {.tabset}

#+ error-effect-plot
error_summary <- d |> summarize_coverage(method, term_type, error_f)

plot_parameter_effect(error_summary, "error_f", "Error Distribution") +
  labs(
    title = "Effect of Error Distribution on Coverage",
    subtitle = "t(6) has heavier tails than Gaussian"
  )

#'
#' ### Error Distribution Effect Table

#+ error-effect-table
make_effect_table(error_summary, "error_f", "Error Distribution")

#'
#' ## Sample Size {.tabset}

#+ n-effect-plot
n_summary <- d |> summarize_coverage(method, term_type, n_f)

plot_parameter_effect(n_summary, "n_f", "Sample Size (n)") +
  labs(
    title = "Effect of Sample Size on Coverage",
    subtitle = "Larger n should improve estimation but may expose bias"
  )

#'
#' # Method Comparison {.tabset .tabset-fade .tabset-pills}
#'
#' ## Overall Performance

#+ method-overall
method_summary <- d |>
  group_by(method) |>
  summarize(
    mean_coverage = mean(coverage, na.rm = TRUE),
    mc_se = sd(coverage, na.rm = TRUE) / sqrt(n()),
    pct_nominal = mean(coverage >= NOMINAL_COVERAGE, na.rm = TRUE) * 100,
    mean_width = mean(mean_width, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_time = mean(fit_time_total[term_type == "ff"], na.rm = TRUE),
    n_dgps = n_distinct(dgp_id),
    .groups = "drop"
  )

method_summary |>
  gt() |>
  tab_header(
    title = "Method Performance Summary",
    subtitle = "Averaged across all applicable DGPs"
  ) |>
  fmt_percent(columns = mean_coverage, decimals = 1) |>
  fmt_number(
    columns = c(pct_nominal, mean_width, mean_rmse, mean_time),
    decimals = 2
  ) |>
  cols_label(
    method = "Method",
    mean_coverage = "Coverage",
    mc_se = "MC SE",
    pct_nominal = "% \u2265 90%",
    mean_width = "Width",
    mean_rmse = "RMSE",
    mean_time = "Time (s)",
    n_dgps = "# DGPs"
  )

#'
#' ## Performance by Term Type

#+ method-by-term, fig.height = 6
ggplot(global_summary, aes(x = method, y = mean_coverage, fill = term_type)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  geom_errorbar(
    aes(
      ymin = mean_coverage - 1.96 * mc_se,
      ymax = mean_coverage + 1.96 * mc_se
    ),
    position = position_dodge(0.8),
    width = 0.2
  ) +
  geom_hline(
    yintercept = NOMINAL_COVERAGE,
    linetype = "dashed",
    color = "red"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0.6, 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Coverage by Method and Term Type",
    x = "Method",
    y = "Coverage",
    fill = "Term Type"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#'
#' ## Well-Specified Performance {.tabset}

#+ well-spec-def
well_spec_conditions <- tribble(
  ~method,
  ~description,
  "pffr",
  "IID (rho=0), homoskedastic, Gaussian",
  "pffr_sandwich",
  "Any correlation, any hetero, Gaussian (cluster-robust)",
  "pffr_hc",
  "IID (rho=0), any hetero, Gaussian (HC robust)",
  "pffr_ar",
  "AR(1), homoskedastic, Gaussian (estimated rho)",
  "pffr_gaulss",
  "IID (rho=0), heteroskedastic, Gaussian (modeled variance)"
)

well_spec_conditions |>
  gt() |>
  tab_header(title = "Well-Specified Conditions by Method")

#+ tag-well-spec
d <- d |>
  mutate(
    well_specified = case_when(
      # pffr: only well-specified under iid, homoskedastic, gaussian
      method == "pffr" &
        corr_type == "iid" &
        hetero_type == "none" &
        error_dist == "gaussian" ~
        TRUE,
      # sandwich: well-specified under all conditions (cluster-robust)
      method == "pffr_sandwich" & error_dist == "gaussian" ~ TRUE,
      # HC: only under iid (handles hetero but not correlation)
      method == "pffr_hc" & corr_type == "iid" & error_dist == "gaussian" ~
        TRUE,
      # AR: well-specified under AR(1), homoskedastic, gaussian
      method == "pffr_ar" &
        corr_type %in% c("ar1", "fourier_pos") &
        hetero_type == "none" &
        error_dist == "gaussian" ~
        TRUE,
      # gaulss: iid + heteroskedastic + gaussian
      method == "pffr_gaulss" &
        corr_type == "iid" &
        hetero_type != "none" &
        error_dist == "gaussian" ~
        TRUE,
      TRUE ~ FALSE
    )
  )

#'
#' ### Well-Specified Results

#+ well-spec-results
well_spec_summary <- d |>
  filter(well_specified) |>
  group_by(method, term_type) |>
  summarize(
    mean_coverage = mean(coverage, na.rm = TRUE),
    mc_se = sd(coverage, na.rm = TRUE) / sqrt(n()),
    pct_nominal = mean(coverage >= NOMINAL_COVERAGE, na.rm = TRUE) * 100,
    mean_z_sd = mean(z_sd, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

ggplot(
  well_spec_summary,
  aes(
    x = term_type,
    y = mean_coverage,
    fill = method,
    color = method
  )
) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  geom_errorbar(
    aes(
      ymin = mean_coverage - 1.96 * mc_se,
      ymax = mean_coverage + 1.96 * mc_se
    ),
    position = position_dodge(0.8),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = NOMINAL_COVERAGE,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0.6, 1)) +
  scale_fill_manual(values = COLORS_METHOD) +
  scale_color_manual(values = COLORS_METHOD) +
  labs(
    title = "Coverage Under Correct Specification",
    subtitle = "Each method on DGPs matching its assumptions",
    x = "Term Type",
    y = "Coverage"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#'
#' ### Well-Specified Table

#+ well-spec-table
well_spec_summary |>
  mutate(coverage_fmt = fmt_coverage(mean_coverage, mc_se)) |>
  select(method, term_type, coverage_fmt, pct_nominal, mean_z_sd) |>
  pivot_wider(
    names_from = term_type,
    values_from = c(coverage_fmt, pct_nominal, mean_z_sd)
  ) |>
  gt() |>
  tab_header(
    title = "Well-Specified Coverage by Method",
    subtitle = "Coverage % (\u00b195% MC interval), z_sd (should be ~1)"
  )

#'
#' ## Misspecification Heatmap

#+ misspec-heatmap, fig.height = 10, fig.width = 12
heatmap_data <- d |>
  group_by(method, term_type, corr_f, hetero_f, error_f) |>
  summarize(mean_coverage = mean(coverage, na.rm = TRUE), .groups = "drop") |>
  mutate(condition = paste(corr_f, hetero_f, error_f, sep = " | "))

ggplot(heatmap_data, aes(x = method, y = condition, fill = mean_coverage)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.0f", mean_coverage * 100)),
    size = 2.2,
    color = "black"
  ) +
  facet_wrap(~term_type, ncol = 3) +
  scale_fill_gradient2(
    low = "#d73027",
    mid = "#ffffbf",
    high = "#1a9850",
    midpoint = 0.85,
    limits = c(0.5, 1),
    labels = scales::percent,
    na.value = "gray90"
  ) +
  labs(
    title = "Coverage Heatmap: All Methods \u00d7 All Conditions",
    subtitle = "Green \u2265 90%, Yellow ~ 85%, Red < 80%",
    x = "Method",
    y = "DGP Condition",
    fill = "Coverage"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold")
  )

#'
#' # pffr vs pffr_sandwich vs pffr_hc {.tabset .tabset-fade .tabset-pills}
#'
#' Direct comparison of the three methods that share the same fitted model
#' but differ only in SE extraction.
#'
#' ## Paired Coverage Difference

#+ trio-comparison
d_trio <- d |>
  filter(method %in% c("pffr", "pffr_sandwich", "pffr_hc")) |>
  select(
    dgp_id,
    rep_id,
    term_type,
    method,
    coverage,
    mean_width,
    rho_f,
    hetero_f,
    corr_f,
    error_f
  ) |>
  pivot_wider(
    names_from = method,
    values_from = c(coverage, mean_width),
    names_sep = "."
  ) |>
  filter(
    !is.na(coverage.pffr),
    !is.na(coverage.pffr_sandwich),
    !is.na(coverage.pffr_hc)
  ) |>
  mutate(
    diff_sandwich = (coverage.pffr_sandwich - coverage.pffr) * 100,
    diff_hc = (coverage.pffr_hc - coverage.pffr) * 100
  )

#+ trio-by-rho, fig.height = 7
d_trio_long <- d_trio |>
  pivot_longer(
    cols = c(diff_sandwich, diff_hc),
    names_to = "comparison",
    values_to = "diff_pp"
  ) |>
  mutate(
    comparison = ifelse(
      comparison == "diff_sandwich",
      "Sandwich - Default",
      "HC - Default"
    )
  )

ggplot(d_trio_long, aes(x = rho_f, y = diff_pp, fill = comparison)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~term_type, ncol = 3) +
  scale_fill_manual(
    values = c(
      "Sandwich - Default" = "#d95f02",
      "HC - Default" = "#e6ab02"
    )
  ) +
  labs(
    title = "Coverage Gain Over Default pffr by Correlation",
    subtitle = "Positive = better than default. HC should NOT help under correlation.",
    x = "Correlation (rho)",
    y = "Coverage difference (pp)",
    fill = NULL
  ) +
  theme(legend.position = "bottom")

#'
#' ## Summary Table

#+ trio-summary-table
d_trio |>
  group_by(rho_f, hetero_f, term_type) |>
  summarize(
    cov_default = mean(coverage.pffr, na.rm = TRUE),
    cov_sandwich = mean(coverage.pffr_sandwich, na.rm = TRUE),
    cov_hc = mean(coverage.pffr_hc, na.rm = TRUE),
    gain_sandwich = mean(diff_sandwich, na.rm = TRUE),
    gain_hc = mean(diff_hc, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) |>
  filter(term_type %in% c("ff", "smooth", "linear", "concurrent")) |>
  gt() |>
  tab_header(
    title = "Coverage: Default vs Sandwich vs HC",
    subtitle = "Gain = coverage improvement in percentage points"
  ) |>
  fmt_percent(columns = c(cov_default, cov_sandwich, cov_hc), decimals = 1) |>
  fmt_number(columns = c(gain_sandwich, gain_hc), decimals = 1)

#'
#' # CI Width Analysis {.tabset .tabset-fade .tabset-pills}
#'
#' ## Width by Method

#+ width-method
width_summary <- d |>
  group_by(method, term_type) |>
  summarize(
    mean_width = mean(mean_width, na.rm = TRUE),
    sd_width = sd(mean_width, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(width_summary, aes(x = term_type, y = mean_width, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = COLORS_METHOD) +
  labs(
    title = "Mean CI Width by Method and Term",
    x = "Term Type",
    y = "Mean CI Width",
    fill = "Method"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#'
#' ## Width vs Coverage Tradeoff

#+ width-coverage-tradeoff, fig.height = 7
tradeoff_data <- d |>
  group_by(method, term_type) |>
  summarize(
    coverage = mean(coverage, na.rm = TRUE),
    width = mean(mean_width, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  tradeoff_data,
  aes(
    x = width,
    y = coverage,
    color = method,
    shape = term_type
  )
) +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_hline(
    yintercept = NOMINAL_COVERAGE,
    linetype = "dashed",
    color = "red"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0.65, 1)) +
  scale_color_manual(values = COLORS_METHOD) +
  labs(
    title = "Coverage vs CI Width Tradeoff",
    subtitle = "Upper-left is ideal: high coverage, narrow CIs",
    x = "Mean CI Width",
    y = "Coverage",
    color = "Method",
    shape = "Term"
  ) +
  theme(legend.position = "right")

#'
#' # Timing Analysis {.tabset}
#'
#' ## Computation Time by Method

#+ timing-summary
timing_data <- d |>
  filter(term_type == "ff") |>
  group_by(method) |>
  summarize(
    mean_time = mean(fit_time_total, na.rm = TRUE),
    median_time = median(fit_time_total, na.rm = TRUE),
    sd_time = sd(fit_time_total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    time_ratio = mean_time / mean_time[method == "pffr"]
  )

timing_data |>
  gt() |>
  tab_header(
    title = "Computation Time by Method",
    subtitle = "Seconds per model fit"
  ) |>
  fmt_number(
    columns = c(mean_time, median_time, sd_time, time_ratio),
    decimals = 1
  )

#'
#' # Summary and Recommendations
#'
#' ## Key Findings

#+ key-findings-table
findings <- d |>
  group_by(method) |>
  summarize(
    overall_cov = mean(coverage, na.rm = TRUE),
    well_spec_cov = mean(coverage[well_specified], na.rm = TRUE),
    misspec_cov = mean(coverage[!well_specified], na.rm = TRUE),
    mean_z_sd = mean(z_sd, na.rm = TRUE),
    mean_time = mean(fit_time_total[term_type == "ff"], na.rm = TRUE),
    .groups = "drop"
  )

findings |>
  gt() |>
  tab_header(title = "Method Summary") |>
  fmt_percent(
    columns = c(overall_cov, well_spec_cov, misspec_cov),
    decimals = 1
  ) |>
  fmt_number(columns = c(mean_z_sd, mean_time), decimals = 2) |>
  cols_label(
    overall_cov = "Overall",
    well_spec_cov = "Well-Spec",
    misspec_cov = "Misspec",
    mean_z_sd = "z_sd",
    mean_time = "Time (s)"
  )

#'
#' ## Recommendations
#'
#' Based on this benchmark:
#'
#' 1. **Default `pffr`**: Adequate under IID, homoskedastic, Gaussian conditions.
#'
#' 2. **`pffr_sandwich` (cluster-robust)**: Recommended as default -- maintains
#'    valid coverage under within-curve autocorrelation AND heteroskedasticity,
#'    with minimal cost under IID.
#'
#' 3. **`pffr_hc` (HC sandwich)**: Does NOT help under correlation (by design).
#'    Only useful for pure heteroskedasticity without correlation.
#'
#' 4. **`pffr_ar` (estimated AR1)**: Useful when correlation is known to be AR(1).
#'
#' 5. **`pffr_gaulss` (modeled variance)**: Specialized for heteroskedasticity
#'    modeling when the variance pattern matters.
#'
#'
#' ---
#'
#' # Study 1: Non-Gaussian Sandwich Coverage
#'
#' This section extends the Gaussian Round 2 benchmark above to **non-Gaussian**
#' families: Poisson (log link) and Binomial (logit link). The DGP includes
#' `ff(X1) + z_lin` with three within-curve correlation structures
#' (IID, AR1(0.9), fourier_pos(0.3)) and two sample sizes (n = 200, 400).
#'
#' Three CI methods are compared:
#'
#' - **default**: standard mgcv Bayesian CIs (no sandwich correction)
#' - **hc**: heteroscedasticity-consistent (HC) sandwich
#' - **cluster**: cluster-robust sandwich (clusters = curves)
#'
#' **Key design notes:**
#'
#' - Poisson DGP: additive eps on linear predictor, log link is collapsible
#'   (conditional beta = marginal beta). Intercept corrected for sigma_eps^2/2.
#' - Binomial DGP: Gaussian copula preserves marginal beta (logit is NOT
#'   collapsible). P(Y=1|X) = p exactly; within-curve correlation via copula.
#' - Grid: 30 x 40 (s x t). Basis: k = c(12, 12) for ff, k = 12 for linear.
#' - 150 replications per DGP cell (12 cells = 2 families x 3 corr x 2 n).
#'

#+ s1-colors
COLORS_S1 <- c(
  default = "#1b9e77",
  hc = "#e6ab02",
  cluster = "#d95f02"
)

# Redefine here in case Study 1 section is rendered standalone
if (!exists("PRACTICAL_THRESHOLD")) PRACTICAL_THRESHOLD <- 0.05

#'
#' ## Study 1: Data Loading
#'

#+ s1-load-data
s1_dir <- ci_path("ci-benchmark/study1-nongaussian")
s1_files <- list.files(
  s1_dir,
  pattern = "^dgp\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)
s1_combined <- file.path(s1_dir, "results_combined.rds")

s1 <- if (file.exists(s1_combined)) {
  readRDS(s1_combined)
} else if (length(s1_files) > 0) {
  bind_rows(lapply(s1_files, readRDS))
} else {
  stop("No Study 1 results found in ", s1_dir)
}

s1 <- s1 |>
  mutate(
    method = factor(method, levels = c("default", "hc", "cluster")),
    term_type = factor(term_type, levels = c("E(Y)", "ff", "linear")),
    family_f = factor(
      family,
      levels = c("poisson", "binomial"),
      labels = c("Poisson", "Binomial")
    ),
    corr_f = factor(
      corr_type,
      levels = c("iid", "ar1", "fourier_pos"),
      labels = c("IID", "AR1(0.9)", "Fourier+(0.3)")
    ),
    n_f = factor(n)
  )

cat("Study 1 loaded:", nrow(s1), "rows\n")
cat("Families:", paste(levels(s1$family_f), collapse = ", "), "\n")
cat("Correlations:", paste(levels(s1$corr_f), collapse = ", "), "\n")
cat("Methods:", paste(levels(s1$method), collapse = ", "), "\n")
cat("Convergence:", sum(s1$converged), "/", nrow(s1), "\n")

#'
#' ## Study 1: Coverage Summary
#'

#+ s1-summary
s1_summary <- s1 |>
  filter(!is.na(coverage)) |>
  group_by(family_f, corr_f, n_f, method, term_type) |>
  summarize(
    coverage = mean(coverage, na.rm = TRUE),
    mc_se = sd(coverage, na.rm = TRUE) / sqrt(n()),
    z_sd = mean(z_sd, na.rm = TRUE),
    width = mean(mean_width, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE),
    bias = mean(bias, na.rm = TRUE),
    fit_time = mean(fit_time, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )

#+ s1-summary-table
s1_summary |>
  select(family_f, corr_f, n_f, method, term_type, coverage, mc_se, z_sd) |>
  mutate(
    coverage_fmt = fmt_coverage(coverage, mc_se),
    z_sd = sprintf("%.2f", z_sd)
  ) |>
  select(family_f, corr_f, n_f, method, term_type, coverage_fmt, z_sd) |>
  pivot_wider(
    names_from = c(method),
    values_from = c(coverage_fmt, z_sd),
    names_glue = "{method}_{.value}"
  ) |>
  gt() |>
  tab_header(
    title = "Study 1: Coverage by Family × Correlation × Method",
    subtitle = "Coverage % (±95% MC interval) and z_sd (target = 1)"
  ) |>
  tab_spanner(label = "default", columns = starts_with("default_")) |>
  tab_spanner(label = "HC", columns = starts_with("hc_")) |>
  tab_spanner(label = "cluster", columns = starts_with("cluster_"))

#'
#' ## Study 1: Coverage Heatmap
#'
#' White = nominal (90%). Red = undercoverage. Blue = overcoverage.
#'

#+ s1-heatmap, fig.width = 12, fig.height = 8
s1_summary |>
  ggplot(aes(
    x = method,
    y = interaction(corr_f, paste0("n=", n_f)),
    fill = coverage
  )) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", coverage * 100)), size = 3) +
  scale_fill_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#4575b4",
    midpoint = NOMINAL_COVERAGE,
    limits = c(0.4, 1),
    labels = scales::percent
  ) +
  facet_grid(family_f ~ term_type) +
  labs(
    title = "Study 1: Coverage Heatmap",
    subtitle = "Non-Gaussian families: Poisson (log) and Binomial (logit)",
    x = "Method",
    y = "Correlation × Sample Size",
    fill = "Coverage"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#'
#' ## Study 1: Coverage Improvement (cluster − default)
#'
#' Paired comparison: how much does cluster-robust improve over default CIs?
#' Positive = cluster helps. Dotted lines at ±5pp practical significance.
#'

#+ s1-paired-diff, fig.width = 12, fig.height = 7
s1_paired <- s1_summary |>
  filter(method %in% c("default", "cluster")) |>
  select(family_f, corr_f, n_f, method, term_type, coverage) |>
  pivot_wider(names_from = method, values_from = coverage) |>
  mutate(diff = cluster - default)

ggplot(
  s1_paired,
  aes(
    x = interaction(corr_f, paste0("n=", n_f)),
    y = diff,
    fill = family_f
  )
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(
    yintercept = c(-PRACTICAL_THRESHOLD, PRACTICAL_THRESHOLD),
    linetype = 3,
    color = "gray50"
  ) +
  facet_wrap(~term_type) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Study 1: Coverage Improvement (cluster − default)",
    subtitle = "Positive = cluster helps. Dotted = ±5pp threshold.",
    x = "Correlation × n",
    y = "Coverage Difference",
    fill = "Family"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#'
#' ## Study 1: Three-Way Method Comparison
#'
#' Direct comparison: default vs HC vs cluster, broken out by correlation.
#' Under IID, all three should be near 90%. Under AR1/Fourier, only cluster
#' should reach nominal.
#'

#+ s1-method-facet, fig.width = 12, fig.height = 8
s1_summary |>
  ggplot(aes(
    x = corr_f,
    y = coverage,
    color = method,
    group = method
  )) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_line(linewidth = 0.6, position = position_dodge(width = 0.4)) +
  geom_errorbar(
    aes(
      ymin = coverage - 1.96 * mc_se,
      ymax = coverage + 1.96 * mc_se
    ),
    width = 0.2,
    linewidth = 0.5,
    position = position_dodge(width = 0.4)
  ) +
  geom_hline(
    yintercept = NOMINAL_COVERAGE,
    linetype = "dashed",
    color = "red"
  ) +
  facet_grid(family_f ~ term_type + n_f) +
  scale_y_continuous(labels = scales::percent, limits = c(0.4, 1)) +
  scale_color_manual(values = COLORS_S1) +
  labs(
    title = "Study 1: Coverage by Correlation Structure",
    subtitle = "Error bars = 95% MC intervals. Red dashed = 90% nominal.",
    x = "Correlation",
    y = "Coverage",
    color = "Method"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

#'
#' ## Study 1: SE Calibration (z_sd)
#'
#' z_sd = SD of (estimate − truth) / SE, averaged over grid points within each
#' rep. Target = 1. Values > 1 mean SEs are too small (anti-conservative);
#' values < 1 mean SEs are too large (conservative).
#'

#+ s1-zsd, fig.width = 12, fig.height = 8
s1_summary |>
  ggplot(aes(x = method, y = z_sd, color = corr_f, shape = term_type)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = c(0.85, 1.15), linetype = 3, color = "gray50") +
  facet_grid(family_f ~ n_f) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Study 1: SE Calibration (z_sd, target = 1)",
    subtitle = "Dotted band = ±0.15 tolerance. Deviations indicate SE miscalibration.",
    x = "Method",
    y = "z_sd",
    color = "Correlation",
    shape = "Term"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

#'
#' ## Study 1: HC vs Cluster Under Correlation
#'
#' HC sandwich corrects for heteroscedasticity but NOT within-curve dependence.
#' Under AR1/Fourier, HC should not reach nominal — only cluster should.
#' This plot shows the coverage gain of each method over default.
#'

#+ s1-hc-vs-cluster, fig.width = 12, fig.height = 7
s1_trio <- s1 |>
  select(dgp_id, rep_id, term_type, method, coverage, family_f, corr_f, n_f) |>
  pivot_wider(
    names_from = method,
    values_from = coverage,
    names_prefix = "cov_"
  ) |>
  mutate(
    gain_hc = (cov_hc - cov_default) * 100,
    gain_cluster = (cov_cluster - cov_default) * 100
  )

s1_trio_long <- s1_trio |>
  pivot_longer(
    cols = c(gain_hc, gain_cluster),
    names_to = "comparison",
    values_to = "gain_pp"
  ) |>
  mutate(
    comparison = factor(
      comparison,
      levels = c("gain_cluster", "gain_hc"),
      labels = c("Cluster − Default", "HC − Default")
    )
  )

ggplot(s1_trio_long, aes(x = corr_f, y = gain_pp, fill = comparison)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(family_f ~ term_type) +
  scale_fill_manual(
    values = c(
      "Cluster − Default" = "#d95f02",
      "HC − Default" = "#e6ab02"
    )
  ) +
  labs(
    title = "Study 1: Coverage Gain Over Default by Correlation",
    subtitle = "HC should NOT help under correlation. Only cluster addresses within-curve dependence.",
    x = "Correlation",
    y = "Coverage Gain (pp)",
    fill = NULL
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

#'
#' ## Study 1: Poisson IID — The Overdispersion Story
#'
#' Under IID errors, default pffr CIs should be near nominal for Binomial
#' but undercover for Poisson. Why? The Poisson DGP adds latent eps on the
#' log-linear predictor, creating overdispersion that mgcv's default SEs
#' ignore. Both HC and cluster fix this since they are model-free variance
#' estimators.
#'

#+ s1-iid-family-comparison, fig.width = 10, fig.height = 6
s1_iid <- s1_summary |>
  filter(corr_f == "IID")

ggplot(s1_iid, aes(x = method, y = coverage, fill = family_f)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(
      ymin = coverage - 1.96 * mc_se,
      ymax = coverage + 1.96 * mc_se
    ),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  geom_hline(
    yintercept = NOMINAL_COVERAGE,
    linetype = "dashed",
    color = "red"
  ) +
  facet_grid(term_type ~ n_f) +
  scale_y_continuous(labels = scales::percent, limits = c(0.6, 1)) +
  labs(
    title = "Study 1: IID Coverage — Poisson vs Binomial",
    subtitle = "IID errors only. Compare Poisson vs Binomial default coverage.",
    x = "Method",
    y = "Coverage",
    fill = "Family"
  ) +
  theme(legend.position = "bottom")

#'
#' ## Study 1: Empirical DGP Verification
#'
#' Before interpreting coverage, verify that the DGPs produce the intended
#' properties: overdispersion for Poisson, and within-curve correlation
#' matching the AR1/Fourier specification.
#'

#+ s1-dgp-verify-setup
# Source Study 1 DGP functions (guarded main blocks won't execute)
if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(refund)
}
source(ci_path("ci-benchmark/benchmark-utils.R"))
source(ci_path("ci-benchmark/confint-benchmark.R"))
source(ci_path("ci-benchmark/sim-study-nongaussian-sandwich.R"))
# Re-assert dplyr::filter after sourcing (confint-benchmark.R can mask it)
filter <- dplyr::filter

#+ s1-dgp-overdispersion, fig.width = 10, fig.height = 5
# Poisson IID: check realized overdispersion
set.seed(42)
dgp_check_results <- map_dfr(1:20, function(rep) {
  # Poisson IID
  sim_pois <- simulate_poisson_ff_linear(
    n = 200,
    corr_type = "iid",
    seed = 42000 + rep
  )
  Y_pois <- sim_pois$data$Y
  # Per-curve mean and variance (across t)
  pois_disp <- tibble(
    curve_mean = rowMeans(Y_pois),
    curve_var = apply(Y_pois, 1, var),
    family = "Poisson IID",
    rep = rep
  )

  # Binomial IID
  sim_binom <- simulate_binomial_ff_linear(
    n = 200,
    corr_type = "iid",
    seed = 42000 + rep
  )
  Y_binom <- sim_binom$data$Y
  p_hat <- rowMeans(Y_binom)
  binom_disp <- tibble(
    curve_mean = p_hat,
    curve_var = apply(Y_binom, 1, var),
    family = "Binomial IID",
    rep = rep
  )

  bind_rows(pois_disp, binom_disp)
})

# Poisson: plot var vs mean, overlay Var=mu line
p_disp_pois <- dgp_check_results |>
  filter(family == "Poisson IID") |>
  ggplot(aes(x = curve_mean, y = curve_var)) +
  geom_point(alpha = 0.05, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Poisson IID: Realized Var vs Mean (per curve)",
    subtitle = "Red = Var(Y) = E(Y) (Poisson assumption). Points above = overdispersion.",
    x = "Curve Mean",
    y = "Curve Variance"
  )

# Compute summary dispersion ratio
pois_disp_ratio <- dgp_check_results |>
  filter(family == "Poisson IID") |>
  summarize(
    mean_dispersion = mean(curve_var / curve_mean, na.rm = TRUE),
    median_dispersion = median(curve_var / curve_mean, na.rm = TRUE),
    q25 = quantile(curve_var / curve_mean, 0.25, na.rm = TRUE),
    q75 = quantile(curve_var / curve_mean, 0.75, na.rm = TRUE)
  )

cat(sprintf(
  "Poisson IID dispersion ratio (Var/Mean): median=%.2f, IQR=[%.2f, %.2f]\n",
  pois_disp_ratio$median_dispersion,
  pois_disp_ratio$q25,
  pois_disp_ratio$q75
))
cat("(Poisson assumption: ratio = 1. Values > 1 confirm overdispersion.)\n\n")

print(p_disp_pois)

#+ s1-dgp-correlation, fig.width = 10, fig.height = 6
# Within-curve correlation: sample residuals from a few datasets
set.seed(43)
corr_check <- map_dfr(c("iid", "ar1", "fourier_pos"), function(ct) {
  cp <- switch(ct, ar1 = 0.9, fourier_pos = 0.3, NULL)
  map_dfr(1:5, function(rep) {
    sim <- simulate_poisson_ff_linear(
      n = 200,
      corr_type = ct,
      corr_param = cp,
      seed = 43000 + rep
    )
    Y <- sim$data$Y
    # Compute residual correlation: subtract row means (crude detrending)
    resid <- Y - rowMeans(Y)
    # Average lag-1 correlation across curves
    lag1_cors <- sapply(seq_len(nrow(resid)), function(i) {
      cor(resid[i, -ncol(resid)], resid[i, -1])
    })

    tibble(
      corr_type = ct,
      rep = rep,
      mean_lag1 = mean(lag1_cors, na.rm = TRUE),
      sd_lag1 = sd(lag1_cors, na.rm = TRUE)
    )
  })
}) |>
  mutate(
    corr_label = factor(
      corr_type,
      levels = c("iid", "ar1", "fourier_pos"),
      labels = c("IID", "AR1(0.9)", "Fourier+(0.3)")
    )
  )

corr_summary <- corr_check |>
  group_by(corr_label) |>
  summarize(
    mean_lag1 = mean(mean_lag1),
    sd_lag1 = mean(sd_lag1),
    .groups = "drop"
  )

cat(
  "Realized lag-1 within-curve correlation (Poisson, averaged over curves):\n"
)
print(corr_summary)
cat("\n(IID should be ~0. AR1/Fourier realized lag-1 < nominal rho=0.9/0.3\n")
cat("because Poisson count variance dilutes latent correlation.)\n")

#+ s1-dgp-correlation-plot, fig.width = 8, fig.height = 4
ggplot(corr_check, aes(x = corr_label, y = mean_lag1)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_lag1 - 1.96 * sd_lag1 / sqrt(200),
      ymax = mean_lag1 + 1.96 * sd_lag1 / sqrt(200)
    ),
    width = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.9, linetype = "dotted", color = "red") +
  labs(
    title = "Study 1: Realized Lag-1 Within-Curve Correlation",
    subtitle = "5 replications each. Error bars = ±1.96 SE across curves. Dotted = nominal AR1.",
    x = "DGP Correlation Type",
    y = "Mean Lag-1 Correlation"
  )

#'
#' ## Study 1: CI Width
#'
#' Cluster-robust CIs are wider (as expected — they account for within-curve
#' dependence). The width premium is the cost of robustness.
#'

#+ s1-width, fig.width = 12, fig.height = 7
s1_summary |>
  ggplot(aes(x = method, y = width, fill = corr_f)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(family_f ~ term_type, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Study 1: Mean CI Width by Method and Correlation",
    subtitle = "Cluster CIs are wider — the cost of accounting for within-curve dependence.",
    x = "Method",
    y = "Mean CI Width",
    fill = "Correlation"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#'
#' ## Study 1: Width vs Coverage Tradeoff
#'

#+ s1-width-coverage, fig.width = 10, fig.height = 7
s1_summary |>
  ggplot(aes(x = width, y = coverage, color = method, shape = corr_f)) +
  geom_point(size = 3.5) +
  geom_hline(
    yintercept = NOMINAL_COVERAGE,
    linetype = "dashed",
    color = "red"
  ) +
  facet_grid(family_f ~ term_type, scales = "free_x") +
  scale_y_continuous(labels = scales::percent, limits = c(0.4, 1)) +
  scale_color_manual(values = COLORS_S1) +
  labs(
    title = "Study 1: Coverage vs CI Width Tradeoff",
    subtitle = "Upper-left ideal (high coverage, narrow CIs).",
    x = "Mean CI Width",
    y = "Coverage",
    color = "Method",
    shape = "Correlation"
  )

#'
#' ## Study 1: Timing
#'

#+ s1-timing
s1_timing <- s1 |>
  filter(term_type == "ff") |>
  group_by(family_f, corr_f, n_f, method) |>
  summarize(
    mean_time = mean(fit_time, na.rm = TRUE),
    median_time = median(fit_time, na.rm = TRUE),
    .groups = "drop"
  )

#' Note: fit_time measures total pffr + CI extraction. The three methods share
#' the same fitted model — timing differences come from SE computation only.
#'

s1_timing |>
  pivot_wider(
    names_from = method,
    values_from = c(mean_time, median_time),
    names_glue = "{method}_{.value}"
  ) |>
  gt() |>
  tab_header(
    title = "Study 1: Computation Time (seconds)",
    subtitle = "Per model fit — methods share the same pffr fit"
  )

#'
#' ## Study 1: Sanity Checks
#'

#+ s1-sanity
cat("=== Study 1 Sanity Checks ===\n\n")

# Check 1: Under AR1, cluster >> default by >10pp
s1_ar1 <- s1_summary |>
  filter(
    corr_f == "AR1(0.9)",
    term_type %in% c("ff", "linear"),
    method %in% c("default", "cluster")
  ) |>
  group_by(method) |>
  summarize(cov = mean(coverage), .groups = "drop") |>
  pivot_wider(names_from = method, values_from = cov)

diff_ar1 <- s1_ar1$cluster - s1_ar1$default
cat(sprintf(
  "AR1: cluster − default = %.1fpp [%s]\n",
  diff_ar1 * 100,
  if (diff_ar1 > 0.10) "PASS" else "FAIL"
))

# Check 2: Under AR1, cluster >> hc
s1_ar1_hc <- s1_summary |>
  filter(
    corr_f == "AR1(0.9)",
    term_type %in% c("ff", "linear"),
    method %in% c("hc", "cluster")
  ) |>
  group_by(method) |>
  summarize(cov = mean(coverage), .groups = "drop") |>
  pivot_wider(names_from = method, values_from = cov)

diff_ar1_hc <- s1_ar1_hc$cluster - s1_ar1_hc$hc
cat(sprintf(
  "AR1: cluster − hc = %.1fpp (cluster should beat HC)\n",
  diff_ar1_hc * 100
))

# Check 3: Poisson IID: hc and cluster should be near 90%
s1_pois_iid <- s1_summary |>
  filter(
    family_f == "Poisson",
    corr_f == "IID",
    method == "cluster",
    term_type == "linear"
  )
cat(sprintf(
  "Poisson IID cluster (linear): %.1f%% [%s]\n",
  mean(s1_pois_iid$coverage) * 100,
  if (abs(mean(s1_pois_iid$coverage) - 0.90) < 0.05) "PASS" else "CHECK"
))

# Check 4: Binomial IID default should be adequate
s1_binom_iid <- s1_summary |>
  filter(
    family_f == "Binomial",
    corr_f == "IID",
    method == "default",
    term_type == "linear"
  )
cat(sprintf(
  "Binomial IID default (linear): %.1f%% [%s]\n",
  mean(s1_binom_iid$coverage) * 100,
  if (abs(mean(s1_binom_iid$coverage) - 0.90) < 0.05) "PASS" else "CHECK"
))

# Check 5: z_sd for cluster under AR1/Fourier should be near 1
s1_zsd_cluster <- s1_summary |>
  filter(
    method == "cluster",
    corr_f %in% c("AR1(0.9)", "Fourier+(0.3)"),
    term_type %in% c("ff", "linear")
  )
cat(sprintf(
  "Cluster z_sd under correlation: mean=%.2f, range=[%.2f, %.2f] [%s]\n",
  mean(s1_zsd_cluster$z_sd),
  min(s1_zsd_cluster$z_sd),
  max(s1_zsd_cluster$z_sd),
  if (all(abs(s1_zsd_cluster$z_sd - 1) < 0.15)) "PASS" else "CHECK"
))

#'
#' ## Study 1: Key Findings
#'
#' ### Cluster sandwich is essential under within-curve correlation
#'
#' Under AR1(0.9) and Fourier+(0.3), default CIs collapse to 44–64% coverage
#' (nominal 90%). Cluster-robust restores coverage to 84–89%:
#'
#' - **Poisson AR1**: default 44% → cluster 87% (+43pp for ff term)
#' - **Binomial AR1**: default 55% → cluster 88% (+33pp for ff term)
#' - **Poisson Fourier**: default 57% → cluster 87% (+29pp for ff term)
#' - **Binomial Fourier**: default 64% → cluster 84% (+21pp for ff term)
#'
#' Coverage gains are comparable for both sample sizes (n=200, 400),
#' suggesting the pattern is driven by the structural mismatch between
#' the assumed and actual covariance rather than by sample size.
#'
#' ### HC sandwich does not fix within-curve correlation
#'
#' HC corrects for heteroscedasticity but NOT within-curve dependence. Under
#' AR1/Fourier, HC barely improves over default. For Binomial AR1, HC is
#' actually *worse* than default (ff: 52% vs 56%; linear: 59% vs 60%). This
#' occurs because HC ignores within-curve correlation in the meat matrix
#' X'eₑ'X, which for high-dimensional functional bases can anti-conservatively
#' bias SE estimates relative to even the naive default. The coverage gain of
#' cluster over HC under AR1 is ~29pp on average.
#'
#' ### Poisson default undercovers even under IID
#'
#' Poisson default coverage under IID is only 69–79% (both terms). The
#' latent eps on the log-linear predictor creates overdispersion that
#' mgcv's default SEs ignore. HC (89–90%) and cluster (87–89%) both fix
#' this since they are model-free variance estimators. This parallels the
#' Gaussian Round 2 finding that the sandwich is not just for correlation.
#'
#' ### Binomial IID: cluster anti-conservative for ff term
#'
#' Binomial ff under IID: cluster gives 76–77% coverage, worse than default
#' (85%). The z_sd for cluster is 1.37 (> 1), indicating SEs are too *small*
#' — this is anti-conservative, not overcorrection. The cluster-robust
#' sandwich with n=200–400 clusters and the high-dimensional ff basis
#' (k = 12 × 12 tensor) may suffer from finite-sample downward bias in the
#' variance estimate, a known issue when the number of parameters is large
#' relative to the number of clusters.
#' The linear term is unaffected (cluster ≈ 89–91% under IID, z_sd ≈ 0.9).
#'
#' ### Cluster coverage gap under correlation (84–89% vs 90%)
#'
#' Even under correlation where cluster is the clear winner, coverage falls
#' 1–6pp below nominal. Possible contributors:
#'
#' 1. **Finite-sample sandwich bias**: n=200 clusters may be insufficient
#'    for asymptotic sandwich theory, especially for the ff term with
#'    effectively ~144 tensor basis functions.
#' 2. **Grid discretization**: The coarse 30 × 40 grid is known to introduce
#'    intercept estimation bias (see Study 2 grid refinement). This may also
#'    affect SE calibration for functional effects.
#' 3. **Penalization interaction**: Sandwich SEs assume a fixed (non-random)
#'    design; the penalty selection introduces additional variability that
#'    the sandwich does not account for.
#'
#' The gap narrows from Poisson (86%) to Binomial (89%), suggesting
#' family-specific residual patterns also matter. Study 2 investigates
#' grid effects; bias-corrected sandwich variants (CR2/CR3) are left for
#' future work.
#'
#' ### Consistency with Gaussian benchmark
#'
#' The non-Gaussian results reinforce the Gaussian Round 2 findings:
#'
#' 1. Cluster-robust sandwich is the right default under within-curve correlation
#' 2. HC is insufficient for dependent functional data
#' 3. The cost of cluster under IID is mild for low-dimensional terms (linear)
#'    but noticeable for high-dimensional ones (ff)
#' 4. Coverage patterns are stable across n=200 and n=400
#'
#' ### Suggested follow-ups
#'
#' 1. **CR2/CR3 sandwich**: Finite-sample bias corrections may close the
#'    remaining coverage gap under correlation and fix the Binomial IID ff issue.
#' 2. **Basis sensitivity**: Rerun one DGP cell with k = 8, 12, 16 for ff to
#'    test whether the anti-conservative bias scales with basis dimension.
#' 3. **Empirical DGP verification**: Compute realized overdispersion ratios
#'    for Poisson and realized within-curve correlations to confirm the
#'    theoretical DGP properties.
#' 4. **Larger n**: Add n = 800 to test whether the cluster coverage gap and
#'    Binomial IID ff issue diminish with more clusters.
#'
#'
#' ---
#'
#' # Study 2: Grid Refinement and Sandwich Covariance Quality
#'
#' This study quantifies whether finer response grids (nygrid) improve
#' cluster-robust sandwich CI calibration, while holding the covariate grid
#' fixed (nxgrid = 60) to isolate the CL2-relevant scaling dimension.
#'
#' Target design (factorial):
#'
#' - **Gaussian family** (fixed), wiggliness = 5
#' - **3 sample sizes**: n = 20, 40, 80
#' - **2 SNR levels**: 3, 25
#' - **3 correlation structures**: IID, AR1(0.9), fourier_pos(0.3)
#' - **Full model**: ff + linear + smooth + concurrent
#' - **Main grid factor**: nygrid in {40, 80, 120} with **nxgrid fixed at 60**
#' - **Methods**: default (no sandwich), cluster-robust sandwich; newer runs may
#'   also include HC and CL2 variants
#' - **120 reps** per DGP x grid cell (planned production target)
#'
#' Note: the downloaded archive is checked for completeness below, and this
#' section reports preliminary results if only a subset of cells completed.
#'

#+ s2-colors
COLORS_S2 <- c(
  default = "#1b9e77",
  cluster = "#d95f02",
  cl2 = "#7570b3",
  hc = "#e6ab02"
)
S2_EXPECTED_N <- c(20L, 40L, 80L)
S2_EXPECTED_SNR <- c(3, 25)
S2_EXPECTED_CORR <- c("iid", "ar1", "fourier_pos")
S2_FIXED_NXGRID <- 60L
S2_EXPECTED_NYGRID <- c(40L, 80L, 120L)
S2_EXPECTED_GRIDS <- paste0("y", S2_EXPECTED_NYGRID)
S2_TARGET_REPS <- 120L

normalize_s2_grid_label <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x == "coarse" ~ "y40",
    x == "medium" ~ "y80",
    x == "fine" ~ "y120",
    TRUE ~ x
  )
}

s2_extract_grid_info <- function(grid_labels) {
  grid_labels <- as.character(grid_labels)
  x_side <- stringr::str_match(grid_labels, "^x(\\d+)_y(\\d+)$")
  x_pair <- stringr::str_match(grid_labels, "^(\\d+)x(\\d+)$")
  y_only <- stringr::str_match(grid_labels, "^y(\\d+)$")

  tibble(grid_label = grid_labels) |>
    mutate(
      nxgrid = dplyr::case_when(
        !is.na(x_side[, 1]) ~ as.integer(x_side[, 2]),
        !is.na(x_pair[, 1]) ~ as.integer(x_pair[, 2]),
        !is.na(y_only[, 1]) ~ S2_FIXED_NXGRID,
        TRUE ~ NA_integer_
      ),
      nygrid = dplyr::case_when(
        !is.na(x_side[, 1]) ~ as.integer(x_side[, 3]),
        !is.na(x_pair[, 1]) ~ as.integer(x_pair[, 3]),
        !is.na(y_only[, 1]) ~ as.integer(y_only[, 2]),
        TRUE ~ NA_integer_
      )
    )
}

format_s2_grid_label <- function(x) {
  x <- as.character(x)
  info <- s2_extract_grid_info(x)
  dplyr::case_when(
    grepl("^y\\d+$", x) ~ paste0("ny = ", info$nygrid),
    grepl("^x\\d+_y\\d+$", x) ~
      paste0("nx = ", info$nxgrid, ", ny = ", info$nygrid),
    grepl("^\\d+x\\d+$", x) ~ gsub("x", "\u00d7", x, fixed = TRUE),
    TRUE ~ x
  )
}

order_s2_grid_levels <- function(grid_levels) {
  s2_extract_grid_info(unique(grid_levels)) |>
    dplyr::arrange(nxgrid, nygrid) |>
    dplyr::pull(grid_label)
}

s2_safe_first <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  x[[1]]
}

s2_fmt_pct <- function(x) {
  if (is.na(x)) return("NA")
  sprintf("%.1f%%", 100 * x)
}

s2_fmt_pp <- function(x) {
  if (is.na(x)) return("NA")
  sprintf("%+.1fpp", 100 * x)
}

#'
#' ## Study 2: Data Loading
#'

#+ s2-load-data
s2_dir <- ci_path("ci-benchmark/study2-grid-refinement")
s2_combined <- file.path(s2_dir, "main_results_combined.rds")
s2_cov_file <- file.path(s2_dir, "cov_quality.rds")
s2_main_dir <- file.path(s2_dir, "main")

s2_new_files <- list.files(
  s2_main_dir,
  pattern = "^dgp\\d+_n\\d+_y\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)
s2_legacy_factorial_files <- list.files(
  s2_main_dir,
  pattern = "^dgp\\d+_n\\d+_grid\\d+x\\d+_rep\\d+\\.rds$",
  full.names = TRUE
)
s2_old_files <- list.files(
  s2_main_dir,
  pattern = "^dgp\\d+_grid\\w+_rep\\d+\\.rds$",
  full.names = TRUE
)

if (
  length(s2_new_files) > 0 &&
    (length(s2_legacy_factorial_files) > 0 || length(s2_old_files) > 0)
) {
  stop(
    "Detected both redesigned Study 2 files (_y###_) and legacy Study 2 files ",
    "in ",
    s2_main_dir,
    ". Remove legacy Study 2 outputs before analysis to ",
    "avoid mixing incompatible designs."
  )
}

s2 <- if (length(s2_new_files) > 0) {
  bind_rows(lapply(s2_new_files, function(f) {
    obj <- readRDS(f)
    if (is.list(obj) && "metrics" %in% names(obj)) obj$metrics else obj
  }))
} else if (length(s2_legacy_factorial_files) > 0) {
  bind_rows(lapply(s2_legacy_factorial_files, function(f) {
    obj <- readRDS(f)
    if (is.list(obj) && "metrics" %in% names(obj)) obj$metrics else obj
  }))
} else if (length(s2_old_files) > 0) {
  bind_rows(lapply(s2_old_files, function(f) {
    obj <- readRDS(f)
    if (is.list(obj) && "metrics" %in% names(obj)) obj$metrics else obj
  }))
} else if (file.exists(s2_combined)) {
  readRDS(s2_combined)
} else {
  stop("No Study 2 results found in ", s2_dir)
}

s2 <- s2 |>
  mutate(grid_label = normalize_s2_grid_label(grid_label))

if (!all(c("nxgrid", "nygrid") %in% names(s2))) {
  s2$nxgrid <- NA_integer_
  s2$nygrid <- NA_integer_
}
grid_dims_from_label <- s2_extract_grid_info(s2$grid_label)
s2$nxgrid <- ifelse(is.na(s2$nxgrid), grid_dims_from_label$nxgrid, s2$nxgrid)
s2$nygrid <- ifelse(is.na(s2$nygrid), grid_dims_from_label$nygrid, s2$nygrid)

s2_grid_levels <- order_s2_grid_levels(s2$grid_label)
s2_grid_labels <- setNames(
  format_s2_grid_label(s2_grid_levels),
  s2_grid_levels
)
s2_grid_extrema <- s2_extract_grid_info(s2_grid_levels) |>
  mutate(n_points = nxgrid * nygrid) |>
  arrange(n_points, nxgrid, nygrid)
s2_coarse_grid <- s2_grid_extrema$grid_label[[1]]
s2_fine_grid <- s2_grid_extrema$grid_label[[nrow(s2_grid_extrema)]]
s2_coarse_display <- unname(s2_grid_labels[s2_coarse_grid])
s2_fine_display <- unname(s2_grid_labels[s2_fine_grid])

s2 <- s2 |>
  mutate(
    method = factor(method, levels = c("default", "cluster", "cl2", "hc")),
    n_f = factor(n, levels = sort(unique(c(S2_EXPECTED_N, n)))),
    term_type = factor(
      term_type,
      levels = c("intercept", "E(Y)", "linear", "concurrent", "smooth", "ff")
    ),
    grid_f = factor(
      grid_label,
      levels = s2_grid_levels,
      labels = unname(s2_grid_labels[s2_grid_levels])
    ),
    corr_f = factor(
      corr_type,
      levels = c("iid", "ar1", "fourier_pos"),
      labels = c("IID", "AR1(0.9)", "Fourier+(0.3)")
    )
  )

s2_obs_reps <- s2 |>
  distinct(n, snr, corr_type, grid_label, rep_id) |>
  count(n, snr, corr_type, grid_label, name = "n_rep")
s2_obs_cells <- s2_obs_reps |>
  select(n, snr, corr_type, grid_label)
s2_expected_cells <- tidyr::crossing(
  n = S2_EXPECTED_N,
  snr = S2_EXPECTED_SNR,
  corr_type = S2_EXPECTED_CORR,
  grid_label = S2_EXPECTED_GRIDS
)
s2_missing_cells <- anti_join(
  s2_expected_cells,
  s2_obs_cells,
  by = c("n", "snr", "corr_type", "grid_label")
)
s2_extra_cells <- anti_join(
  s2_obs_cells,
  s2_expected_cells,
  by = c("n", "snr", "corr_type", "grid_label")
)
s2_complete <- nrow(s2_missing_cells) == 0 &&
  nrow(s2_extra_cells) == 0 &&
  nrow(s2_obs_reps) > 0 &&
  min(s2_obs_reps$n_rep) >= S2_TARGET_REPS

cat("Study 2 loaded:", nrow(s2), "rows\n")
cat("Grid levels:", paste(levels(s2$grid_f), collapse = ", "), "\n")
cat("Fixed nxgrid (expected main design):", S2_FIXED_NXGRID, "\n")
cat("Methods:", paste(levels(s2$method), collapse = ", "), "\n")
cat("Correlations:", paste(levels(s2$corr_f), collapse = ", "), "\n")
cat(
  "Reps per cell:",
  paste(
    range(table(paste(
      s2$dgp_id,
      s2$grid_label,
      s2$method,
      s2$term_type
    ))),
    collapse = "-"
  ),
  "\n"
)
cat(
  "Expected Study 2 cells:",
  nrow(s2_expected_cells),
  "(target reps/cell:",
  S2_TARGET_REPS,
  ")\n"
)
cat("Observed Study 2 cells:", nrow(s2_obs_cells), "\n")
cat("Missing expected cells:", nrow(s2_missing_cells), "\n")
if (nrow(s2_obs_reps) > 0) {
  cat(
    "Observed reps per cell (range):",
    paste(range(s2_obs_reps$n_rep), collapse = "-"),
    "\n"
  )
}
cat(
  "Study 2 completeness:",
  if (s2_complete) "COMPLETE" else "INCOMPLETE",
  "\n"
)
if (nrow(s2_missing_cells) > 0) {
  cat("First missing cells:\n")
  print(head(s2_missing_cells, 10))
}

# Load covariance quality metrics if available
s2_cov <- if (file.exists(s2_cov_file)) readRDS(s2_cov_file) else NULL
if (!is.null(s2_cov))
  cat("Covariance quality metrics loaded:", nrow(s2_cov), "rows\n")

#'
#' ## Study 2: Coverage Summary
#'

#+ s2-summary
s2_summary_by_n <- s2 |>
  dplyr::filter(!is.na(coverage)) |>
  group_by(n_f, corr_f, grid_f, method, term_type) |>
  summarize(
    mc_se = sd(coverage, na.rm = TRUE) / sqrt(n()),
    coverage = mean(coverage, na.rm = TRUE),
    z_sd = mean(z_sd, na.rm = TRUE),
    z_mean_sd = sd(z_mean, na.rm = TRUE),
    z_mean = mean(z_mean, na.rm = TRUE),
    width = mean(mean_width, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE),
    bias = mean(bias, na.rm = TRUE),
    fit_time = mean(fit_time, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )

s2_summary <- s2_summary_by_n |>
  group_by(corr_f, grid_f, method, term_type) |>
  summarize(
    mc_se = sd(coverage, na.rm = TRUE) / sqrt(n()),
    coverage = mean(coverage, na.rm = TRUE),
    z_sd = mean(z_sd, na.rm = TRUE),
    z_mean_sd = mean(z_mean_sd, na.rm = TRUE),
    z_mean = mean(z_mean, na.rm = TRUE),
    width = mean(width, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE),
    bias = mean(bias, na.rm = TRUE),
    fit_time = mean(fit_time, na.rm = TRUE),
    n_reps = sum(n_reps, na.rm = TRUE),
    .groups = "drop"
  )

#+ s2-summary-table
s2_summary |>
  dplyr::filter(term_type %in% c("intercept", "ff", "linear", "E(Y)")) |>
  select(corr_f, grid_f, method, term_type, coverage, mc_se, z_sd) |>
  mutate(
    coverage_fmt = fmt_coverage(coverage, mc_se),
    z_sd = sprintf("%.2f", z_sd)
  ) |>
  select(corr_f, grid_f, method, term_type, coverage_fmt, z_sd) |>
  pivot_wider(
    names_from = method,
    values_from = c(coverage_fmt, z_sd),
    names_glue = "{method}_{.value}"
  ) |>
  gt() |>
  tab_header(
    title = "Study 2: Coverage by Grid Level \u00d7 Correlation \u00d7 Method",
    subtitle = "Coverage % (\u00b195% MC interval) and z_sd (target = 1), averaged over n"
  ) |>
  tab_spanner(label = "default", columns = starts_with("default_")) |>
  tab_spanner(label = "cluster", columns = starts_with("cluster_"))

#'
#' ## Study 2: Coverage Heatmap
#'
#' Does finer grid resolution improve coverage? White = nominal (90%).
#'

#+ s2-heatmap, fig.width = 12, fig.height = 10
s2_summary |>
  dplyr::filter(
    term_type %in%
      c("intercept", "ff", "linear", "smooth", "concurrent", "E(Y)")
  ) |>
  ggplot(aes(
    x = method,
    y = interaction(corr_f, grid_f, sep = " | "),
    fill = coverage
  )) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", coverage * 100)), size = 2.8) +
  scale_fill_gradient2(
    low = "#d73027",
    mid = "white",
    high = "#4575b4",
    midpoint = NOMINAL_COVERAGE,
    limits = c(0.5, 1),
    labels = scales::percent
  ) +
  facet_wrap(~term_type, ncol = 3) +
  labs(
    title = "Study 2: Coverage Heatmap by Grid Level",
    subtitle = "Paired design (same data subsampled), averaged over n; nxgrid fixed at 60 in redesigned Study 2.",
    x = "Method",
    y = "Correlation | nygrid",
    fill = "Coverage"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#'
#' ## Study 2: Coverage by Grid Level (Line Plot)
#'
#' The key visual: does coverage improve monotonically with grid resolution?
#'

#+ s2-coverage-by-grid, fig.width = 12, fig.height = 8
s2_summary_by_n |>
  dplyr::filter(
    term_type %in%
      c("intercept", "ff", "linear", "smooth", "concurrent", "E(Y)")
  ) |>
  ggplot(aes(
    x = grid_f,
    y = coverage,
    color = method,
    group = method
  )) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_errorbar(
    aes(
      ymin = coverage - 1.96 * mc_se,
      ymax = coverage + 1.96 * mc_se
    ),
    width = 0.2,
    linewidth = 0.5,
    position = position_dodge(width = 0.3)
  ) +
  geom_hline(
    yintercept = NOMINAL_COVERAGE,
    linetype = "dashed",
    color = "red"
  ) +
  facet_grid(n_f + corr_f ~ term_type) +
  scale_y_continuous(labels = scales::percent, limits = c(0.3, 1)) +
  scale_color_manual(values = COLORS_S2) +
  labs(
    title = "Study 2: Coverage by nygrid (n shown explicitly)",
    subtitle = "Red dashed = 90% nominal. Error bars = 95% MC intervals. nxgrid fixed at 60.",
    x = "Response Grid (nygrid)",
    y = "Coverage",
    color = "Method"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

#'
#' ## Study 2: Coverage Improvement (cluster - default) by Grid
#'
#' Paired comparison: how much does cluster-robust improve over default,
#' and does the improvement change with grid resolution?
#'

#+ s2-paired-diff, fig.width = 12, fig.height = 7
build_s2_paired_diff_summary <- function(s2_data, num_method, den_method) {
  s2_data |>
    dplyr::filter(
      method %in% c(den_method, num_method),
      term_type %in%
        c("intercept", "ff", "linear", "smooth", "concurrent", "E(Y)")
    ) |>
    group_by(dgp_id, rep_id, n_f, term_type, method, corr_f, grid_f) |>
    summarize(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(
      names_from = method,
      values_from = coverage,
      values_fn = mean
    ) |>
    mutate(
      diff = .data[[num_method]] - .data[[den_method]]
    ) |>
    dplyr::filter(!is.na(diff)) |>
    group_by(n_f, corr_f, grid_f, term_type) |>
    summarize(
      mean_diff = mean(diff, na.rm = TRUE),
      se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

plot_s2_paired_diff <- function(summary_df, title, subtitle) {
  ggplot(
    summary_df,
    aes(x = grid_f, y = mean_diff, fill = corr_f)
  ) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(
      aes(
        ymin = mean_diff - 1.96 * se_diff,
        ymax = mean_diff + 1.96 * se_diff
      ),
      position = position_dodge(width = 0.7),
      width = 0.2
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(n_f ~ term_type) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Response Grid (nygrid)",
      y = "Coverage Difference",
      fill = "Correlation"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

s2_paired_summary <- build_s2_paired_diff_summary(
  s2_data = s2,
  num_method = "cluster",
  den_method = "default"
)

plot_s2_paired_diff(
  s2_paired_summary,
  title = "Study 2: Coverage Gain (cluster \u2212 default) by Grid Level",
  subtitle = "Positive = cluster helps. n shown explicitly; nxgrid fixed at 60."
)

#'
#' ## Study 2: Coverage Improvement (cl2 - cluster) by Grid
#'
#' Paired comparison: incremental gain (or loss) from CL2 relative to cluster.
#'

#+ s2-paired-diff-cl2, fig.width = 12, fig.height = 7
if (sum(as.character(s2$method) == "cl2", na.rm = TRUE) == 0) {
  cat("No CL2 results available in Study 2 extract.\n")
} else {
  s2_paired_summary_cl2 <- build_s2_paired_diff_summary(
    s2_data = s2,
    num_method = "cl2",
    den_method = "cluster"
  )

  if (nrow(s2_paired_summary_cl2) == 0) {
    cat("No paired CL2 vs cluster comparisons available in Study 2 extract.\n")
  } else {
    plot_s2_paired_diff(
      s2_paired_summary_cl2,
      title = "Study 2: Coverage Gain (CL2 \u2212 cluster) by Grid Level",
      subtitle = "Positive = CL2 improves coverage over cluster. n shown explicitly; nxgrid fixed at 60."
    )
  }
}

#'
#' ## Study 2: SE Calibration (z_sd) by Grid
#'
#' z_sd = SD of (estimate - truth) / SE averaged over grid points.
#' Target = 1. The key question: does finer grid bring z_sd closer to 1?
#'

#+ s2-zsd-by-grid, fig.width = 12, fig.height = 8
s2_summary_by_n |>
  dplyr::filter(
    term_type %in%
      c("intercept", "ff", "linear", "smooth", "concurrent", "E(Y)")
  ) |>
  ggplot(aes(
    x = grid_f,
    y = z_sd,
    color = method,
    group = method
  )) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_hline(
    yintercept = c(0.85, 1.15),
    linetype = "dotted",
    color = "gray70"
  ) +
  facet_grid(n_f + corr_f ~ term_type) +
  scale_color_manual(values = COLORS_S2) +
  labs(
    title = "Study 2: SE Calibration (z_sd) by nygrid (n shown explicitly)",
    subtitle = "Dotted band = \u00b10.15 tolerance. z_sd > 1 = SEs too small. nxgrid fixed at 60.",
    x = "Response Grid (nygrid)",
    y = "z_sd",
    color = "Method"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

#'
#' ## Study 2: Intercept Diagnostics by Grid
#'
#' On the smallest response grid (nygrid = 40), intercept coverage can drop
#' sharply despite
#' rep-averaged z_mean near 0 and z_sd near 1. The paradox resolves when
#' we look at between-rep variance: sd(z_mean) is large on coarse grids,
#' meaning individual reps have globally-shifted intercepts that cancel
#' in the average but cause low coverage per rep.
#'

#+ s2-intercept-bias, fig.width = 10, fig.height = 9
s2_intercept <- s2_summary_by_n |>
  dplyr::filter(term_type == "intercept")

p_zmean_sd <- ggplot(
  s2_intercept,
  aes(x = grid_f, y = z_mean_sd, color = method, group = method)
) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_grid(n_f ~ corr_f) +
  scale_color_manual(values = COLORS_S2) +
  labs(
    title = "Intercept: Between-Rep Variance of z_mean",
    subtitle = "sd(z_mean) across reps. Large values = global intercept shifts driving undercoverage.",
    x = "Response Grid (nygrid)",
    y = "sd(z_mean) across reps",
    color = "Method"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

p_zsd <- ggplot(
  s2_intercept,
  aes(x = grid_f, y = z_sd, color = method, group = method)
) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  facet_grid(n_f ~ corr_f) +
  scale_color_manual(values = COLORS_S2) +
  labs(
    title = "Intercept z_sd by Grid Resolution",
    subtitle = "Within-rep z_sd (target = 1). Conditionally well-calibrated even on coarse grids.",
    x = "Response Grid (nygrid)",
    y = "z_sd",
    color = "Method"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

p_cov <- ggplot(
  s2_intercept,
  aes(x = grid_f, y = coverage, color = method, group = method)
) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_errorbar(
    aes(
      ymin = coverage - 1.96 * mc_se,
      ymax = coverage + 1.96 * mc_se
    ),
    width = 0.2,
    linewidth = 0.5,
    position = position_dodge(width = 0.3)
  ) +
  geom_hline(
    yintercept = NOMINAL_COVERAGE,
    linetype = "dashed",
    color = "red"
  ) +
  facet_grid(n_f ~ corr_f) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_manual(values = COLORS_S2) +
  labs(
    title = "Intercept Coverage by Grid Resolution",
    subtitle = "nygrid effect by n (nxgrid fixed at 60).",
    x = "Response Grid (nygrid)",
    y = "Coverage",
    color = "Method"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

p_zmean_sd / p_zsd / p_cov

#'
#' ## Study 2: Covariance Quality Metrics
#'
#' Diagonal ratio = mean(model SE^2) / MC Var(est - truth). Values < 1 mean
#' model SEs underestimate true variability. The cluster sandwich should have
#' ratio closer to 1 than default under correlation.
#'

#+ s2-cov-quality, fig.width = 12, fig.height = 8
if (!is.null(s2_cov)) {
  s2_cov_plot <- s2_cov |>
    mutate(
      grid_label = normalize_s2_grid_label(grid_label),
      n_f = factor(n, levels = sort(unique(c(S2_EXPECTED_N, n)))),
      grid_f = factor(
        grid_label,
        levels = s2_grid_levels,
        labels = unname(s2_grid_labels[s2_grid_levels])
      ),
      corr_f = factor(
        corr_type,
        levels = c("iid", "ar1", "fourier_pos"),
        labels = c("IID", "AR1(0.9)", "Fourier+(0.3)")
      ),
      term_type = factor(
        term_type,
        levels = c("intercept", "linear", "concurrent", "smooth", "ff")
      )
    )

  ggplot(
    s2_cov_plot,
    aes(x = grid_f, y = diag_ratio_median, color = method, group = method)
  ) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
    geom_errorbar(
      aes(ymin = diag_ratio_iqr_low, ymax = diag_ratio_iqr_high),
      width = 0.2,
      linewidth = 0.5,
      position = position_dodge(width = 0.3)
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    facet_grid(n_f + corr_f ~ term_type) +
    scale_color_manual(values = COLORS_S2) +
    labs(
      title = "Study 2: Diagonal Variance Ratio by nygrid (n shown explicitly)",
      subtitle = "Median (IQR) of mean(SE\u00b2) / Var(est \u2212 truth). Target = 1. nxgrid fixed at 60.",
      x = "Response Grid (nygrid)",
      y = "Diagonal Ratio",
      color = "Method"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
} else {
  cat("No covariance quality metrics available.\n")
}

#'
#' ## Study 2: CI Width by Grid
#'

#+ s2-width, fig.width = 12, fig.height = 8
s2_summary_by_n |>
  dplyr::filter(term_type %in% c("ff", "linear", "smooth", "concurrent")) |>
  ggplot(aes(x = grid_f, y = width, color = method, group = method)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  facet_grid(n_f + corr_f ~ term_type, scales = "free_y") +
  scale_color_manual(values = COLORS_S2) +
  labs(
    title = "Study 2: Mean CI Width by nygrid (n shown explicitly)",
    subtitle = "Cluster CIs are wider under correlation; nxgrid fixed at 60.",
    x = "Response Grid (nygrid)",
    y = "Mean CI Width",
    color = "Method"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

#'
#' ## Study 2: CL2 Inference Overhead (Targeted Timing Benchmark)
#'
#' Separate timing benchmark for method-specific `coef()` inference cost on a
#' simplified IID Gaussian model (`ff + linear`, low wiggliness, high SNR). This
#' isolates the incremental CL2 cost better than raw `fit_time`.
#'

#+ s2-cl2-inference-overhead, fig.width = 12, fig.height = 9
s2_cl2_timing_dir <- file.path("ci-benchmark", "study2-cl2-timing")
s2_cl2_timing_raw_file <- file.path(s2_cl2_timing_dir, "raw_timings.rds")
s2_cl2_timing_overhead_csv <- file.path(
  s2_cl2_timing_dir,
  "summary_overheads.csv"
)
s2_cl2_nx_side_dir <- file.path("ci-benchmark", "study2-nx-sidecheck")
s2_cl2_nx_side_raw_file <- file.path(s2_cl2_nx_side_dir, "raw_timings.rds")
s2_cl2_nx_side_overhead_csv <- file.path(
  s2_cl2_nx_side_dir,
  "summary_overheads.csv"
)

summarize_cl2_timing_pairs <- function(raw_timings) {
  raw_timings |>
    dplyr::filter(ok %in% TRUE | is.na(ok)) |>
    dplyr::filter(method %in% c("default", "cluster", "hc", "cl2")) |>
    mutate(grid_label = normalize_s2_grid_label(grid_label)) |>
    select(
      rep_id,
      seed,
      n,
      snr,
      wiggliness,
      corr_type,
      grid_label,
      nxgrid,
      nygrid,
      method,
      fit_time,
      infer_time,
      total_time
    ) |>
    pivot_wider(
      names_from = method,
      values_from = c(infer_time, total_time),
      names_sep = "__",
      values_fn = mean
    ) |>
    mutate(
      cl2_minus_cluster = infer_time__cl2 - infer_time__cluster,
      cluster_minus_default = infer_time__cluster - infer_time__default,
      cl2_over_cluster_ratio = infer_time__cl2 / infer_time__cluster,
      cl2_overhead_share_total = cl2_minus_cluster / total_time__cluster
    ) |>
    dplyr::filter(
      is.finite(cl2_minus_cluster),
      is.finite(cl2_over_cluster_ratio),
      is.finite(cl2_overhead_share_total)
    )
}

summarize_cl2_timing_overhead <- function(pair_df) {
  pair_df |>
    group_by(n, grid_label, nxgrid, nygrid) |>
    summarize(
      n_rep = n(),
      mean_cl2_minus_cluster = mean(cl2_minus_cluster, na.rm = TRUE),
      mc_se_cl2_minus_cluster = sd(cl2_minus_cluster, na.rm = TRUE) /
        sqrt(n_rep),
      mean_cl2_over_cluster_ratio = mean(cl2_over_cluster_ratio, na.rm = TRUE),
      mc_se_cl2_over_cluster_ratio = sd(cl2_over_cluster_ratio, na.rm = TRUE) /
        sqrt(n_rep),
      mean_cl2_overhead_share_total = mean(
        cl2_overhead_share_total,
        na.rm = TRUE
      ),
      mc_se_cl2_overhead_share_total = sd(
        cl2_overhead_share_total,
        na.rm = TRUE
      ) /
        sqrt(n_rep),
      .groups = "drop"
    )
}

if (file.exists(s2_cl2_timing_raw_file)) {
  s2_cl2_pairs <- summarize_cl2_timing_pairs(readRDS(s2_cl2_timing_raw_file))

  if (nrow(s2_cl2_pairs) == 0) {
    cat("No valid paired CL2 timing results available.\n")
  } else {
    s2_cl2_grid_levels <- order_s2_grid_levels(unique(s2_cl2_pairs$grid_label))
    s2_cl2_grid_labels <- setNames(
      format_s2_grid_label(s2_cl2_grid_levels),
      s2_cl2_grid_levels
    )

    s2_cl2_overhead <- summarize_cl2_timing_overhead(s2_cl2_pairs) |>
      mutate(
        n_f = factor(n, levels = sort(unique(n))),
        grid_f = factor(
          grid_label,
          levels = s2_cl2_grid_levels,
          labels = unname(s2_cl2_grid_labels[s2_cl2_grid_levels])
        )
      ) |>
      arrange(grid_f)

    p_delta <- ggplot(
      s2_cl2_overhead,
      aes(x = grid_f, y = mean_cl2_minus_cluster)
    ) +
      geom_col(fill = "#c44e52", width = 0.65) +
      geom_errorbar(
        aes(
          ymin = mean_cl2_minus_cluster - 1.96 * mc_se_cl2_minus_cluster,
          ymax = mean_cl2_minus_cluster + 1.96 * mc_se_cl2_minus_cluster
        ),
        width = 0.2
      ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~n_f, nrow = 1) +
      labs(
        title = "CL2 Incremental Inference Cost (CL2 - cluster)",
        subtitle = "Paired per-fit `coef()` inference time difference; main timing benchmark (nxgrid fixed at 60)",
        x = "Response Grid (nygrid)",
        y = "Seconds"
      ) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    p_ratio <- ggplot(
      s2_cl2_overhead,
      aes(x = grid_f, y = mean_cl2_over_cluster_ratio)
    ) +
      geom_col(fill = "#4c72b0", width = 0.65) +
      geom_errorbar(
        aes(
          ymin = mean_cl2_over_cluster_ratio -
            1.96 * mc_se_cl2_over_cluster_ratio,
          ymax = mean_cl2_over_cluster_ratio +
            1.96 * mc_se_cl2_over_cluster_ratio
        ),
        width = 0.2
      ) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      facet_wrap(~n_f, nrow = 1) +
      labs(
        title = "Relative Inference Cost (CL2 / cluster)",
        subtitle = "Paired per-fit ratio on inference stage only; nxgrid fixed at 60",
        x = "Response Grid (nygrid)",
        y = "Ratio"
      ) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    p_share <- ggplot(
      s2_cl2_overhead,
      aes(x = grid_f, y = mean_cl2_overhead_share_total)
    ) +
      geom_col(fill = "#55a868", width = 0.65) +
      geom_errorbar(
        aes(
          ymin = mean_cl2_overhead_share_total -
            1.96 * mc_se_cl2_overhead_share_total,
          ymax = mean_cl2_overhead_share_total +
            1.96 * mc_se_cl2_overhead_share_total
        ),
        width = 0.2
      ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~n_f, nrow = 1) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        title = "CL2 Overhead as Share of Total Runtime",
        subtitle = "Share of total time under cluster baseline (`fit + cluster inference`)",
        x = "Response Grid (nygrid)",
        y = "Share of total runtime"
      ) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    p_delta / p_ratio / p_share

    s2_cl2_overhead |>
      select(
        n_f,
        grid_f,
        n_rep,
        mean_cl2_minus_cluster,
        mean_cl2_over_cluster_ratio,
        mean_cl2_overhead_share_total
      ) |>
      gt() |>
      tab_header(
        title = "Study 2: CL2 Inference Overhead by Grid (Targeted Benchmark)",
        subtitle = "Simplified IID Gaussian `ff + linear`; values are paired CL2 vs cluster summaries"
      ) |>
      cols_label(
        n_f = "n",
        grid_f = "Grid",
        n_rep = "Reps",
        mean_cl2_minus_cluster = "CL2 - cluster (s)",
        mean_cl2_over_cluster_ratio = "CL2 / cluster",
        mean_cl2_overhead_share_total = "CL2 overhead / total"
      ) |>
      fmt_number(
        columns = c(mean_cl2_minus_cluster, mean_cl2_over_cluster_ratio),
        decimals = 2
      ) |>
      fmt_percent(columns = mean_cl2_overhead_share_total, decimals = 1)
  }
} else if (file.exists(s2_cl2_timing_overhead_csv)) {
  s2_cl2_overhead_csv <- read.csv(
    s2_cl2_timing_overhead_csv,
    stringsAsFactors = FALSE
  )
  if (!("n" %in% names(s2_cl2_overhead_csv))) {
    s2_cl2_overhead_csv$n <- NA_integer_
  }
  s2_cl2_overhead_csv <- s2_cl2_overhead_csv |>
    mutate(
      grid_label = normalize_s2_grid_label(grid_label),
      n_f = factor(n, levels = sort(unique(n)))
    )
  s2_cl2_grid_levels <- order_s2_grid_levels(unique(
    s2_cl2_overhead_csv$grid_label
  ))
  s2_cl2_grid_labels <- setNames(
    format_s2_grid_label(s2_cl2_grid_levels),
    s2_cl2_grid_levels
  )

  s2_cl2_overhead_csv <- s2_cl2_overhead_csv |>
    mutate(
      grid_f = factor(
        grid_label,
        levels = s2_cl2_grid_levels,
        labels = unname(s2_cl2_grid_labels[s2_cl2_grid_levels])
      )
    ) |>
    arrange(grid_f)

  s2_cl2_overhead_csv |>
    select(
      n_f,
      grid_f,
      n_rep,
      mean_cl2_minus_cluster,
      mean_cl2_over_cluster_ratio
    ) |>
    gt() |>
    tab_header(
      title = "Study 2: CL2 Inference Overhead by Grid (CSV Summary)",
      subtitle = "Raw timing file not found; showing pre-aggregated summaries"
    ) |>
    cols_label(
      n_f = "n",
      grid_f = "Grid",
      n_rep = "Reps",
      mean_cl2_minus_cluster = "CL2 - cluster (s)",
      mean_cl2_over_cluster_ratio = "CL2 / cluster"
    ) |>
    fmt_number(
      columns = c(mean_cl2_minus_cluster, mean_cl2_over_cluster_ratio),
      decimals = 2
    )
} else {
  cat(
    "No targeted CL2 timing benchmark results found in ci-benchmark/study2-cl2-timing.\n"
  )
}

#'
#' ## Study 2: nxgrid Side-Check for CL2 Timing (Optional)
#'
#' Tiny side-check to verify nxgrid is secondary for CL2 timing, run at fixed
#' n and nygrid. Not part of the main Study 2 design.
#'

#+ s2-cl2-inference-overhead-nx-sidecheck, fig.width = 11, fig.height = 6
if (file.exists(s2_cl2_nx_side_raw_file)) {
  s2_cl2_nx_pairs <- summarize_cl2_timing_pairs(readRDS(
    s2_cl2_nx_side_raw_file
  ))

  if (nrow(s2_cl2_nx_pairs) == 0) {
    cat("No valid nxgrid side-check timing results available.\n")
  } else {
    s2_cl2_nx_overhead <- summarize_cl2_timing_overhead(s2_cl2_nx_pairs) |>
      mutate(
        n_f = factor(n, levels = sort(unique(n))),
        nx_f = factor(nxgrid, levels = sort(unique(nxgrid))),
        ny_f = factor(nygrid, levels = sort(unique(nygrid)))
      )

    p_nx_side <- ggplot(
      s2_cl2_nx_overhead,
      aes(x = nx_f, y = mean_cl2_minus_cluster)
    ) +
      geom_col(fill = "#8172b2", width = 0.65) +
      geom_errorbar(
        aes(
          ymin = mean_cl2_minus_cluster - 1.96 * mc_se_cl2_minus_cluster,
          ymax = mean_cl2_minus_cluster + 1.96 * mc_se_cl2_minus_cluster
        ),
        width = 0.2
      ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(n_f ~ ny_f, labeller = label_both) +
      labs(
        title = "CL2 Timing nxgrid Side-Check (Optional)",
        subtitle = "Paired CL2 - cluster inference time at fixed nygrid (not part of main Study 2)",
        x = "Covariate Grid (nxgrid)",
        y = "CL2 - cluster (seconds)"
      )

    p_nx_side

    s2_cl2_nx_overhead |>
      select(
        n_f,
        nygrid,
        nxgrid,
        n_rep,
        mean_cl2_minus_cluster,
        mean_cl2_over_cluster_ratio,
        mean_cl2_overhead_share_total
      ) |>
      gt() |>
      tab_header(
        title = "Study 2: nxgrid Side-Check (CL2 Timing Overhead)",
        subtitle = "Separate side-check; use for sensitivity only, not main Study 2 conclusions"
      ) |>
      cols_label(
        n_f = "n",
        nygrid = "nygrid",
        nxgrid = "nxgrid",
        n_rep = "Reps",
        mean_cl2_minus_cluster = "CL2 - cluster (s)",
        mean_cl2_over_cluster_ratio = "CL2 / cluster",
        mean_cl2_overhead_share_total = "CL2 overhead / total"
      ) |>
      fmt_number(
        columns = c(mean_cl2_minus_cluster, mean_cl2_over_cluster_ratio),
        decimals = 2
      ) |>
      fmt_percent(columns = mean_cl2_overhead_share_total, decimals = 1)
  }
} else if (file.exists(s2_cl2_nx_side_overhead_csv)) {
  s2_cl2_nx_overhead_csv <- read.csv(
    s2_cl2_nx_side_overhead_csv,
    stringsAsFactors = FALSE
  )
  if (!("n" %in% names(s2_cl2_nx_overhead_csv))) {
    s2_cl2_nx_overhead_csv$n <- NA_integer_
  }
  s2_cl2_nx_overhead_csv <- s2_cl2_nx_overhead_csv |>
    mutate(
      n_f = factor(n, levels = sort(unique(n)))
    ) |>
    arrange(n, nygrid, nxgrid)

  s2_cl2_nx_overhead_csv |>
    select(
      n_f,
      nygrid,
      nxgrid,
      n_rep,
      mean_cl2_minus_cluster,
      mean_cl2_over_cluster_ratio
    ) |>
    gt() |>
    tab_header(
      title = "Study 2: nxgrid Side-Check (CSV Summary)",
      subtitle = "Raw timing file not found; showing pre-aggregated summaries"
    ) |>
    cols_label(
      n_f = "n",
      nygrid = "nygrid",
      nxgrid = "nxgrid",
      n_rep = "Reps",
      mean_cl2_minus_cluster = "CL2 - cluster (s)",
      mean_cl2_over_cluster_ratio = "CL2 / cluster"
    ) |>
    fmt_number(
      columns = c(mean_cl2_minus_cluster, mean_cl2_over_cluster_ratio),
      decimals = 2
    )
} else {
  cat(
    "No nxgrid side-check timing results found in ci-benchmark/study2-nx-sidecheck.\n"
  )
}

#'
#' ## Study 2: Sanity Checks
#'

#+ s2-sanity
cat("=== Study 2 Sanity / Completeness Checks ===\n\n")

cat(
  "Expected cells:",
  nrow(s2_expected_cells),
  "| Observed cells:",
  nrow(s2_obs_cells),
  "| Missing cells:",
  nrow(s2_missing_cells),
  "\n"
)
if (nrow(s2_obs_reps) > 0) {
  cat(
    "Observed reps per cell (min-max):",
    paste(range(s2_obs_reps$n_rep), collapse = "-"),
    "\n"
  )
}
cat("Completeness flag:", if (s2_complete) "COMPLETE" else "INCOMPLETE", "\n")

if (nrow(s2_missing_cells) > 0) {
  cat("\nFirst missing expected cells:\n")
  print(head(s2_missing_cells, 12), n = 12)
}

cat("\nObserved cells with rep counts:\n")
print(
  s2_obs_reps |>
    arrange(n, snr, corr_type, grid_label),
  n = min(30, nrow(s2_obs_reps))
)

s2_fine_diff <- s2_summary |>
  filter(
    term_type %in% c("ff", "linear"),
    grid_f == s2_fine_display
  ) |>
  select(corr_f, term_type, method, coverage) |>
  tidyr::pivot_wider(names_from = method, values_from = coverage) |>
  mutate(cluster_minus_default = cluster - default)

if (nrow(s2_fine_diff) > 0) {
  cat("\nFine-grid coverage gain (cluster - default):\n")
  print(s2_fine_diff, n = nrow(s2_fine_diff))
} else {
  cat(
    "\nFine-grid cluster-default comparison not available in current extract.\n"
  )
}

#'
#' ## Study 2: Key Findings
#'
#' These findings are preliminary unless the completeness check above reports
#' **COMPLETE**. The current summary is computed from all available cells and
#' should be interpreted as a pilot snapshot when cells/reps are missing.

#+ s2-key-findings
n_cells_obs <- nrow(s2_obs_cells)
n_cells_exp <- nrow(s2_expected_cells)
rep_range_text <- if (nrow(s2_obs_reps) > 0) {
  paste(range(s2_obs_reps$n_rep), collapse = "-")
} else {
  "none"
}

cat("1) Data completeness\n")
cat(
  sprintf(
    "   Observed %d / %d expected cells, reps per observed cell: %s (target: %d).\n",
    n_cells_obs,
    n_cells_exp,
    rep_range_text,
    S2_TARGET_REPS
  )
)

cluster_default_fine <- s2_summary |>
  filter(
    term_type %in% c("ff", "linear"),
    grid_f == s2_fine_display
  ) |>
  select(term_type, method, coverage) |>
  group_by(term_type, method) |>
  summarize(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = method, values_from = coverage) |>
  mutate(diff = cluster - default)

cat("2) Fine-grid cluster vs default (ff, linear)\n")
if (nrow(cluster_default_fine) == 0) {
  cat("   Not available in the current extract.\n")
} else {
  apply(cluster_default_fine, 1, function(row) {
    cat(
      sprintf(
        "   %s: default=%s, cluster=%s, cluster-default=%s\n",
        row[["term_type"]],
        s2_fmt_pct(as.numeric(row[["default"]])),
        s2_fmt_pct(as.numeric(row[["cluster"]])),
        s2_fmt_pp(as.numeric(row[["diff"]]))
      )
    )
  })
}

cluster_intercept_grid <- s2_summary |>
  filter(method == "cluster", term_type == "intercept") |>
  group_by(grid_f) |>
  summarize(coverage = mean(coverage, na.rm = TRUE), .groups = "drop")
s2_intercept_coarse <- cluster_intercept_grid |>
  filter(grid_f == s2_coarse_display) |>
  pull(coverage) |>
  s2_safe_first()
s2_intercept_fine <- cluster_intercept_grid |>
  filter(grid_f == s2_fine_display) |>
  pull(coverage) |>
  s2_safe_first()

cat("3) Intercept coverage change from coarsest to finest observed grid\n")
cat(
  sprintf(
    "   cluster intercept: %s -> %s (%s)\n",
    s2_fmt_pct(s2_intercept_coarse),
    s2_fmt_pct(s2_intercept_fine),
    s2_fmt_pp(s2_intercept_fine - s2_intercept_coarse)
  )
)
#'

#' ---
#'
#' *Report generated on `r Sys.time()`*
