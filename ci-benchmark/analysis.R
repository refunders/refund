#' ---
#' title: "pffr CI Benchmark: Round 2 Analysis"
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

#+ load-packages
library(tidyverse)
library(gt)
library(patchwork)

theme_set(theme_minimal(base_size = 11))

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
#' # Theory
#'
#' This benchmark evaluates pointwise confidence intervals for functional effects
#' in a function-on-function additive model. Data are generated from
#' \eqn{Y_i(t) = \eta_i(t) + \epsilon_i(t)} where \eqn{\eta_i(t)} combines:
#' a functional intercept, a function-on-function term \eqn{\int X_i(s)\beta(s,t)ds},
#' a varying-coefficient linear term, a smooth scalar-varying term, and a
#' concurrent term.
#'
#' DGP difficulty is controlled by:
#' - sample size (\code{n}),
#' - signal-to-noise ratio (\code{snr}),
#' - truth wiggliness (\code{wiggliness}),
#' - error distribution (\code{gaussian}, \code{t6}),
#' - residual covariance structure: IID, AR(1), or periodic non-monotone
#'   (\code{fourier_pos}),
#' - heteroskedasticity pattern over \eqn{t}.
#'
#' Coverage is assessed pointwise on evaluation grids for each term type.
#'
#' CI estimators compared:
#' - \code{pffr}: default mgcv/refund Bayesian covariance-based pointwise CIs.
#' - \code{pffr_hc}: observation-level HC sandwich covariance (heteroskedasticity
#'   robust, not cluster/correlation robust).
#' - \code{pffr_sandwich}: cluster-robust sandwich covariance (clusters = curves),
#'   targeting within-curve dependence and heteroskedasticity.
#' - \code{pffr_ar}: \code{bam(..., rho=...)} AR(1) working-correlation fit, with
#'   CIs from the fitted AR model covariance.
#' - \code{pffr_gaulss}: Gaussian location-scale fit, modeling variance over
#'   \eqn{t}, with pointwise CIs from the fitted gaulss covariance.
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
#' ---
#'
#' *Report generated on `r Sys.time()`*
