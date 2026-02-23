# ===========================================================================
# Study 2 CL2 Inference Timing (Standalone Script)
# ===========================================================================
#
# Purpose:
#   Measure method-specific post-fit inference time (especially CL2 overhead)
#   for a simpler pffr model (ff + linear), aligned with the Study 2 main design:
#   factorial in n x nygrid with nxgrid fixed. Includes an optional tiny nxgrid
#   side-check (separate output directory) to verify nxgrid is secondary.
#
# Default main timing design:
#   - Gaussian family
#   - IID errors
#   - low wiggliness (0.01)
#   - high SNR (25)
#   - simpler model: ff + linear
#   - n in {20, 40, 80}
#   - nygrid in {40, 80, 120}
#   - nxgrid fixed at 60
#
# Timed quantities:
#   - fit_time: base pffr fit (same for all inference methods on that dataset)
#   - infer_time: time for coef(..., sandwich = <method>) on the fitted model
#   - paired overheads (e.g., cl2 - cluster) computed post hoc
#
# Usage examples:
#   Rscript ci-benchmark/study2-cl2-inference-timing.R
#   Rscript ci-benchmark/study2-cl2-inference-timing.R --n-rep=5
#   Rscript ci-benchmark/study2-cl2-inference-timing.R --n-levels=20,40
#   Rscript ci-benchmark/study2-cl2-inference-timing.R --mode=nx-sidecheck
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
source("ci-benchmark/sim-study-grid-refinement.R")

# Constants -------------------------------------------------------------------

TIMING_METHODS <- c(
  default = "none",
  hc = "hc",
  cluster = "cluster",
  cl2 = "cl2"
)
TIMING_BASE_SEED <- 91001L

# Argument parsing ------------------------------------------------------------

parse_int_vec <- function(x) {
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_args <- function(args) {
  out <- list(
    mode = "main",
    n_rep = 10L,
    n_levels = c(20L, 40L, 80L),
    snr = 25,
    wiggliness = 0.01,
    nygrid_levels = STUDY2_NYGRID_LEVELS,
    nxgrid = STUDY2_FIXED_NXGRID,
    out_dir = "ci-benchmark/study2-cl2-timing",
    nx_side_n = 40L,
    nx_side_snr = 25,
    nx_side_wiggliness = 0.01,
    nx_side_nygrid = 80L,
    nx_side_nxgrid_levels = STUDY2_NXGRID_LEVELS_ALL,
    out_dir_nx_sidecheck = "ci-benchmark/study2-nx-sidecheck",
    k_ff = c(12L, 12L),
    k_yindex = 12L
  )

  if (length(args) == 0) return(out)

  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) >= 2) kv[2] else ""

    if (key == "mode" && nzchar(val)) out$mode <- val
    if (key == "n-rep" && nzchar(val)) out$n_rep <- as.integer(val)
    if (key == "n-levels" && nzchar(val)) out$n_levels <- parse_int_vec(val)
    if (key == "n" && nzchar(val)) out$n_levels <- as.integer(val)
    if (key == "snr" && nzchar(val)) out$snr <- as.numeric(val)
    if (key == "wiggliness" && nzchar(val)) out$wiggliness <- as.numeric(val)
    if (key == "nygrids" && nzchar(val)) out$nygrid_levels <- parse_int_vec(val)
    if (key == "nxgrid" && nzchar(val)) out$nxgrid <- as.integer(val)
    if (key == "out-dir" && nzchar(val)) out$out_dir <- val

    if (key == "nx-side-n" && nzchar(val)) out$nx_side_n <- as.integer(val)
    if (key == "nx-side-snr" && nzchar(val)) out$nx_side_snr <- as.numeric(val)
    if (key == "nx-side-wiggliness" && nzchar(val)) {
      out$nx_side_wiggliness <- as.numeric(val)
    }
    if (key == "nx-side-nygrid" && nzchar(val)) {
      out$nx_side_nygrid <- as.integer(val)
    }
    if (key == "nx-side-nxgrids" && nzchar(val)) {
      out$nx_side_nxgrid_levels <- parse_int_vec(val)
    }
    if (key == "out-dir-nx-sidecheck" && nzchar(val)) {
      out$out_dir_nx_sidecheck <- val
    }

    if (key == "k-ff" && nzchar(val)) {
      k_ff <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
      if (length(k_ff) == 2) out$k_ff <- k_ff
    }
    if (key == "k-yindex" && nzchar(val)) out$k_yindex <- as.integer(val)
  }

  out
}

validate_config <- function(cfg) {
  stopifnot(
    "mode must be 'main' or 'nx-sidecheck'" = cfg$mode %in%
      c("main", "nx-sidecheck"),
    "n_rep must be positive" = is.numeric(cfg$n_rep) && cfg$n_rep >= 1,
    "snr must be positive" = is.numeric(cfg$snr) && cfg$snr > 0,
    "wiggliness must be non-negative" = is.numeric(cfg$wiggliness) &&
      cfg$wiggliness >= 0,
    "nxgrid must be positive" = is.numeric(cfg$nxgrid) && cfg$nxgrid >= 1,
    "n_levels must be non-empty" = length(cfg$n_levels) >= 1,
    "nygrid_levels must be non-empty" = length(cfg$nygrid_levels) >= 1
  )

  cfg$n_rep <- as.integer(cfg$n_rep)
  cfg$n_levels <- sort(unique(as.integer(cfg$n_levels)))
  cfg$nygrid_levels <- sort(unique(as.integer(cfg$nygrid_levels)))
  cfg$nxgrid <- as.integer(cfg$nxgrid)
  cfg$nx_side_n <- as.integer(cfg$nx_side_n)
  cfg$nx_side_nygrid <- as.integer(cfg$nx_side_nygrid)
  cfg$nx_side_nxgrid_levels <- sort(unique(as.integer(
    cfg$nx_side_nxgrid_levels
  )))
  cfg$k_ff <- as.integer(cfg$k_ff)
  cfg$k_yindex <- as.integer(cfg$k_yindex)

  if (any(cfg$nygrid_levels > STUDY2_FINE_NY)) {
    stop("All nygrid levels must be <= STUDY2_FINE_NY (", STUDY2_FINE_NY, ")")
  }
  if (cfg$nxgrid > STUDY2_FINE_NX) {
    stop("nxgrid must be <= STUDY2_FINE_NX (", STUDY2_FINE_NX, ")")
  }
  if (any(cfg$nx_side_nxgrid_levels > STUDY2_FINE_NX)) {
    stop(
      "All nx-side nxgrid levels must be <= STUDY2_FINE_NX (",
      STUDY2_FINE_NX,
      ")"
    )
  }
  if (cfg$nx_side_nygrid > STUDY2_FINE_NY) {
    stop("nx-side nygrid must be <= STUDY2_FINE_NY (", STUDY2_FINE_NY, ")")
  }

  cfg
}

# Grid specs ------------------------------------------------------------------

build_main_timing_grids <- function(nygrid_levels, nxgrid) {
  grids <- lapply(as.integer(nygrid_levels), function(ny) {
    list(nxgrid = as.integer(nxgrid), nygrid = as.integer(ny))
  })
  names(grids) <- vapply(
    as.integer(nygrid_levels),
    study2_grid_label_y,
    character(1)
  )
  grids
}

build_nx_sidecheck_grids <- function(nxgrid_levels, nygrid) {
  grids <- lapply(as.integer(nxgrid_levels), function(nx) {
    list(nxgrid = as.integer(nx), nygrid = as.integer(nygrid))
  })
  names(grids) <- sprintf(
    "x%d_y%d",
    as.integer(nxgrid_levels),
    as.integer(nygrid)
  )
  grids
}

# Data generation -------------------------------------------------------------

generate_timing_sims <- function(n, snr, wiggliness, seed, grid_specs) {
  sim_fine <- generate_benchmark_data(
    n = n,
    nxgrid = STUDY2_FINE_NX,
    nygrid = STUDY2_FINE_NY,
    snr = snr,
    error_dist = "gaussian",
    corr_type = "iid",
    corr_param = NA_real_,
    hetero_type = "none",
    hetero_param = NA_real_,
    terms = c("ff", "linear"),
    wiggliness = wiggliness,
    seed = seed
  )

  sims <- list()
  for (grid_label in names(grid_specs)) {
    g <- grid_specs[[grid_label]]
    if (g$nxgrid == STUDY2_FINE_NX && g$nygrid == STUDY2_FINE_NY) {
      sims[[grid_label]] <- sim_fine
    } else {
      sims[[grid_label]] <- subsample_grid(
        sim = sim_fine,
        target_nxgrid = g$nxgrid,
        target_nygrid = g$nygrid
      )
    }
  }
  sims
}

# Timing helpers --------------------------------------------------------------

fit_base_model <- function(sim, k_ff = c(12L, 12L), k_yindex = 12L) {
  formula <- build_pffr_formula(sim, sim$s_grid, k_smooth = 12L, k_ff = k_ff)
  bs_yindex <- list(bs = "ps", k = k_yindex, m = c(2, 1))

  t0 <- Sys.time()
  fit <- pffr(
    formula,
    yind = sim$t_grid,
    data = sim$data,
    bs.yindex = bs_yindex,
    sandwich = "none"
  )
  fit_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  list(fit = fit, fit_time = fit_time)
}

time_coef_inference <- function(fit, nxgrid, nygrid, sandwich_type) {
  t0 <- Sys.time()
  res <- try(
    coef(
      fit,
      sandwich = sandwich_type,
      seWithMean = FALSE,
      n1 = nxgrid,
      n2 = nygrid
    ),
    silent = TRUE
  )
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  list(
    ok = !inherits(res, "try-error"),
    infer_time = dt,
    error_msg = if (inherits(res, "try-error")) as.character(res) else
      NA_character_
  )
}

time_one_grid <- function(sim, grid_label, grid_info, rep_id, seed, n, cfg) {
  fit_res <- fit_base_model(sim, k_ff = cfg$k_ff, k_yindex = cfg$k_yindex)
  fit <- fit_res$fit
  fit_time <- fit_res$fit_time

  rows <- vector("list", length(TIMING_METHODS))
  i <- 1L
  for (method_name in names(TIMING_METHODS)) {
    sw <- TIMING_METHODS[[method_name]]
    infer <- time_coef_inference(
      fit = fit,
      nxgrid = grid_info$nxgrid,
      nygrid = grid_info$nygrid,
      sandwich_type = sw
    )

    rows[[i]] <- tibble(
      rep_id = rep_id,
      seed = seed,
      n = as.integer(n),
      snr = cfg$snr_current,
      wiggliness = cfg$wiggliness_current,
      corr_type = "iid",
      nxgrid = grid_info$nxgrid,
      nygrid = grid_info$nygrid,
      grid_label = grid_label,
      method = method_name,
      sandwich = sw,
      fit_time = fit_time,
      infer_time = infer$infer_time,
      total_time = fit_time + infer$infer_time,
      ok = infer$ok,
      error_msg = infer$error_msg
    )
    i <- i + 1L
  }

  bind_rows(rows)
}

# Summaries -------------------------------------------------------------------

summarize_timing_by_method <- function(raw_timings) {
  raw_timings |>
    group_by(n, grid_label, nxgrid, nygrid, method) |>
    summarize(
      n_rep = n(),
      n_ok = sum(ok, na.rm = TRUE),
      mean_fit_time = mean(fit_time, na.rm = TRUE),
      mean_infer_time = mean(infer_time, na.rm = TRUE),
      median_infer_time = median(infer_time, na.rm = TRUE),
      sd_infer_time = sd(infer_time, na.rm = TRUE),
      mean_total_time = mean(total_time, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(n, nygrid, nxgrid, method)
}

summarize_paired_overheads <- function(raw_timings) {
  paired <- raw_timings |>
    dplyr::filter(ok) |>
    dplyr::select(
      rep_id,
      n,
      grid_label,
      nxgrid,
      nygrid,
      method,
      infer_time,
      total_time
    ) |>
    tidyr::pivot_wider(
      names_from = method,
      values_from = c(infer_time, total_time),
      names_sep = "__"
    )

  paired |>
    dplyr::mutate(
      cluster_minus_default = infer_time__cluster - infer_time__default,
      cl2_minus_cluster = infer_time__cl2 - infer_time__cluster,
      cl2_over_cluster_ratio = infer_time__cl2 / infer_time__cluster,
      cl2_overhead_share_total = cl2_minus_cluster / total_time__cluster,
      hc_minus_default = infer_time__hc - infer_time__default
    ) |>
    dplyr::group_by(n, grid_label, nxgrid, nygrid) |>
    dplyr::summarize(
      n_rep = n(),
      mean_cluster_minus_default = mean(cluster_minus_default, na.rm = TRUE),
      mean_cl2_minus_cluster = mean(cl2_minus_cluster, na.rm = TRUE),
      median_cl2_minus_cluster = median(cl2_minus_cluster, na.rm = TRUE),
      mean_cl2_over_cluster_ratio = mean(cl2_over_cluster_ratio, na.rm = TRUE),
      mean_cl2_overhead_share_total = mean(
        cl2_overhead_share_total,
        na.rm = TRUE
      ),
      mean_hc_minus_default = mean(hc_minus_default, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(n, nygrid, nxgrid)
}

# Main runners ----------------------------------------------------------------

run_timing_grid_benchmark <- function(
  n_levels,
  grid_specs,
  cfg,
  out_dir,
  study_label
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat(study_label, "\n")
  cat(strrep("=", nchar(study_label)), "\n", sep = "")
  cat("n_rep:", cfg$n_rep, "\n")
  cat("n levels:", paste(n_levels, collapse = ", "), "\n")
  cat("snr:", cfg$snr_current, "\n")
  cat("wiggliness:", cfg$wiggliness_current, "\n")
  cat("grids:", paste(names(grid_specs), collapse = ", "), "\n")
  cat("methods:", paste(names(TIMING_METHODS), collapse = ", "), "\n")
  cat("output:", out_dir, "\n\n")

  raw <- list()
  idx <- 1L
  n_total <- length(n_levels) * cfg$n_rep * length(grid_specs)

  for (n_val in n_levels) {
    for (rep_id in seq_len(cfg$n_rep)) {
      seed <- TIMING_BASE_SEED + 10000L * as.integer(n_val) + rep_id
      sims <- generate_timing_sims(
        n = n_val,
        snr = cfg$snr_current,
        wiggliness = cfg$wiggliness_current,
        seed = seed,
        grid_specs = grid_specs
      )

      for (grid_label in names(grid_specs)) {
        g <- grid_specs[[grid_label]]
        cat(sprintf(
          "\r[%d/%d] n=%d rep=%d grid=%s",
          idx,
          n_total,
          n_val,
          rep_id,
          grid_label
        ))
        raw[[idx]] <- tryCatch(
          time_one_grid(
            sim = sims[[grid_label]],
            grid_label = grid_label,
            grid_info = g,
            rep_id = rep_id,
            seed = seed,
            n = n_val,
            cfg = cfg
          ),
          error = function(e) {
            tibble(
              rep_id = rep_id,
              seed = seed,
              n = as.integer(n_val),
              snr = cfg$snr_current,
              wiggliness = cfg$wiggliness_current,
              corr_type = "iid",
              nxgrid = g$nxgrid,
              nygrid = g$nygrid,
              grid_label = grid_label,
              method = NA_character_,
              sandwich = NA_character_,
              fit_time = NA_real_,
              infer_time = NA_real_,
              total_time = NA_real_,
              ok = FALSE,
              error_msg = conditionMessage(e)
            )
          }
        )
        idx <- idx + 1L
      }
    }
  }
  cat("\n")

  raw_timings <- bind_rows(raw)
  summary_method <- summarize_timing_by_method(raw_timings)
  summary_overhead <- summarize_paired_overheads(raw_timings)

  saveRDS(raw_timings, file.path(out_dir, "raw_timings.rds"))
  write.csv(
    summary_method,
    file.path(out_dir, "summary_by_method.csv"),
    row.names = FALSE
  )
  write.csv(
    summary_overhead,
    file.path(out_dir, "summary_overheads.csv"),
    row.names = FALSE
  )

  cat("\n=== Mean inference time by n/grid/method (seconds) ===\n")
  print(
    summary_method |>
      select(
        n,
        grid_label,
        nxgrid,
        nygrid,
        method,
        n_rep,
        n_ok,
        mean_infer_time,
        median_infer_time
      ),
    n = nrow(summary_method)
  )

  cat("\n=== Paired overhead summaries (seconds; ratio unitless) ===\n")
  print(summary_overhead, n = nrow(summary_overhead))

  invisible(list(
    raw_timings = raw_timings,
    summary_by_method = summary_method,
    summary_overheads = summary_overhead
  ))
}

run_main_timing_study <- function(cfg) {
  cfg$snr_current <- cfg$snr
  cfg$wiggliness_current <- cfg$wiggliness
  grid_specs <- build_main_timing_grids(
    nygrid_levels = cfg$nygrid_levels,
    nxgrid = cfg$nxgrid
  )

  run_timing_grid_benchmark(
    n_levels = cfg$n_levels,
    grid_specs = grid_specs,
    cfg = cfg,
    out_dir = cfg$out_dir,
    study_label = sprintf(
      "Study 2 CL2 inference timing (n x nygrid; nxgrid fixed = %d)",
      cfg$nxgrid
    )
  )
}

run_nx_sidecheck_timing_study <- function(cfg) {
  cfg$snr_current <- cfg$nx_side_snr
  cfg$wiggliness_current <- cfg$nx_side_wiggliness
  grid_specs <- build_nx_sidecheck_grids(
    nxgrid_levels = cfg$nx_side_nxgrid_levels,
    nygrid = cfg$nx_side_nygrid
  )

  run_timing_grid_benchmark(
    n_levels = cfg$nx_side_n,
    grid_specs = grid_specs,
    cfg = cfg,
    out_dir = cfg$out_dir_nx_sidecheck,
    study_label = sprintf(
      "Study 2 CL2 timing nxgrid side-check (n=%d, nygrid=%d)",
      cfg$nx_side_n,
      cfg$nx_side_nygrid
    )
  )
}

if (sys.nframe() == 0) {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  cfg <- validate_config(cfg)

  if (identical(cfg$mode, "nx-sidecheck")) {
    run_nx_sidecheck_timing_study(cfg)
  } else {
    run_main_timing_study(cfg)
  }
}
