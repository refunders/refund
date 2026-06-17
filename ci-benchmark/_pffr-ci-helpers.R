# _pffr-ci-helpers.R
# Thin data-loading and helper script sourced by pffr-ci-report.qmd.
# Extracts reusable pieces from analysis.R without executing any plots.

# Packages ------------------------------------------------------------------

library(tidyverse)
library(gt)
library(patchwork)

theme_set(theme_minimal(base_size = 11))

# Path helper ---------------------------------------------------------------

ci_path <- function(path) {
  if (file.exists(path) || dir.exists(path)) return(path)
  alt <- sub("^ci-benchmark/", "", path)
  if (file.exists(alt) || dir.exists(alt)) return(alt)
  path
}

# Constants -----------------------------------------------------------------

NOMINAL_COVERAGE <- 0.90
PRACTICAL_THRESHOLD <- 0.05

COLORS_S1 <- c(
  default = "#1b9e77",
  hc = "#e6ab02",
  cluster = "#d95f02",
  cl2 = "#7570b3"
)

COLORS_S2 <- c(
  default = "#1b9e77",
  cluster = "#d95f02",
  cl2 = "#7570b3",
  hc = "#e6ab02"
)

# Formatting helpers --------------------------------------------------------

fmt_coverage <- function(cov, se) {
  sprintf("%.1f%% (\u00b1%.1f)", cov * 100, 1.96 * se * 100)
}

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

# Study 2 grid label helpers ------------------------------------------------

S2_EXPECTED_N <- c(20L, 40L, 80L)
S2_EXPECTED_SNR <- c(3, 25)
S2_EXPECTED_CORR <- c("iid", "ar1", "fourier_pos")
S2_FIXED_NXGRID <- 60L
S2_EXPECTED_NYGRID <- c(40L, 80L, 120L)
S2_EXPECTED_GRIDS <- paste0("y", S2_EXPECTED_NYGRID)
S2_TARGET_REPS <- 50L

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

# Study 1 data loading ------------------------------------------------------

load_study1 <- function() {
  s1_dir <- ci_path("ci-benchmark/study1-nongaussian")
  s1_files <- list.files(
    s1_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  s1_combined_cl2 <- file.path(s1_dir, "results_combined_with_cl2.rds")
  s1_combined <- file.path(s1_dir, "results_combined.rds")

  s1 <- if (file.exists(s1_combined_cl2)) {
    readRDS(s1_combined_cl2)
  } else if (file.exists(s1_combined)) {
    readRDS(s1_combined)
  } else if (length(s1_files) > 0) {
    bind_rows(lapply(s1_files, readRDS))
  } else {
    stop("No non-Gaussian results found in ", s1_dir)
  }

  s1 |>
    mutate(
      method = factor(method, levels = c("default", "hc", "cluster", "cl2")),
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
}

# Study 2 data loading ------------------------------------------------------

load_study2 <- function() {
  s2_dir <- ci_path("ci-benchmark/study2-grid-refinement")
  s2_combined <- file.path(s2_dir, "main_results_combined.rds")
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
      "Detected both redesigned Study 2 files and legacy files in ",
      s2_main_dir,
      ". Remove legacy outputs before analysis."
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
  grid_dims <- s2_extract_grid_info(s2$grid_label)
  s2$nxgrid <- ifelse(is.na(s2$nxgrid), grid_dims$nxgrid, s2$nxgrid)
  s2$nygrid <- ifelse(is.na(s2$nygrid), grid_dims$nygrid, s2$nygrid)

  s2_grid_levels <- order_s2_grid_levels(s2$grid_label)
  s2_grid_labels <- setNames(
    format_s2_grid_label(s2_grid_levels),
    s2_grid_levels
  )

  s2 <- s2 |>
    mutate(
      method = factor(method, levels = c("default", "cluster", "hc", "cl2")),
      n_f = factor(n, levels = sort(unique(c(S2_EXPECTED_N, n)))),
      snr_f = factor(
        snr,
        levels = sort(unique(c(S2_EXPECTED_SNR, snr))),
        labels = paste0("SNR = ", sort(unique(c(S2_EXPECTED_SNR, snr))))
      ),
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

  attr(s2, "grid_levels") <- s2_grid_levels
  attr(s2, "grid_labels") <- s2_grid_labels
  s2
}

load_study2_cov <- function(s2_grid_levels, s2_grid_labels) {
  s2_cov_file <- ci_path(
    file.path("ci-benchmark/study2-grid-refinement", "cov_quality.rds")
  )
  if (!file.exists(s2_cov_file)) return(NULL)

  readRDS(s2_cov_file) |>
    mutate(
      grid_label = normalize_s2_grid_label(grid_label),
      n_f = factor(n, levels = sort(unique(c(S2_EXPECTED_N, n)))),
      snr_f = factor(
        snr,
        levels = sort(unique(c(S2_EXPECTED_SNR, snr))),
        labels = paste0("SNR = ", sort(unique(c(S2_EXPECTED_SNR, snr))))
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
      ),
      term_type = factor(
        term_type,
        levels = c("intercept", "linear", "concurrent", "smooth", "ff")
      ),
      method = factor(method, levels = c("default", "cluster", "hc", "cl2"))
    ) |>
    dplyr::filter(method %in% c("default", "cluster", "hc", "cl2"))
}

# Study 2 paired-diff helper ------------------------------------------------

build_s2_paired_diff_summary <- function(s2_data, num_method, den_method) {
  s2_data |>
    dplyr::filter(
      method %in% c(den_method, num_method),
      term_type %in% c("ff", "linear", "smooth", "concurrent")
    ) |>
    group_by(dgp_id, rep_id, n_f, snr_f, corr_f, grid_f, term_type, method) |>
    summarize(coverage = mean(coverage, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(
      names_from = method,
      values_from = coverage,
      values_fn = mean
    ) |>
    mutate(diff = .data[[num_method]] - .data[[den_method]]) |>
    dplyr::filter(!is.na(diff)) |>
    group_by(snr_f, corr_f, grid_f, term_type) |>
    summarize(
      mean_diff = mean(diff, na.rm = TRUE),
      se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}
