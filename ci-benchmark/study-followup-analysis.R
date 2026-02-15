# ===========================================================================
# Combined Analysis: Follow-up Studies 1 & 2
# ===========================================================================
#
# Reads outputs from both follow-up studies and produces:
#   - Coverage/width/z_sd tables (gt format with kable fallback)
#   - Heatmaps: method x family x corr_type for Study 1
#   - Paired coverage difference plots for both studies
#   - Per-term covariance quality comparison plots for Study 2
#   - Concise markdown summary
#
# Usage:
#   Rscript ci-benchmark/study-followup-analysis.R
# ===========================================================================

# Setup -----------------------------------------------------------------------

library(tidyverse)
library(ggplot2)

theme_set(theme_bw(base_size = 12))

NOMINAL_COVERAGE <- 0.90
PRACTICAL_THRESHOLD <- 0.05 # 5 percentage points

# Load results ----------------------------------------------------------------

load_study1 <- function(dir = "ci-benchmark/study1-nongaussian") {
  path <- file.path(dir, "results_combined.rds")
  if (!file.exists(path)) {
    # Try loading from individual files
    files <- list.files(
      dir,
      pattern = "^dgp\\d+_rep\\d+\\.rds$",
      full.names = TRUE
    )
    if (length(files) == 0) stop("No Study 1 results found in ", dir)
    results <- bind_rows(lapply(files, readRDS))
    return(results)
  }
  readRDS(path)
}

load_study2 <- function(dir = "ci-benchmark/study2-grid-refinement") {
  main_path <- file.path(dir, "main_results_combined.rds")
  cov_path <- file.path(dir, "cov_quality.rds")

  results <- if (file.exists(main_path)) {
    readRDS(main_path)
  } else {
    # Try loading from individual files
    main_dir <- file.path(dir, "main")
    files <- list.files(
      main_dir,
      pattern = "^dgp\\d+_grid\\w+_rep\\d+\\.rds$",
      full.names = TRUE
    )
    if (length(files) == 0)
      return(list(results = tibble(), cov_quality = tibble()))
    bind_rows(lapply(files, function(f) {
      obj <- readRDS(f)
      if (is.list(obj) && "metrics" %in% names(obj)) obj$metrics else obj
    }))
  }

  cov_quality <- if (file.exists(cov_path)) readRDS(cov_path) else tibble()

  list(results = results, cov_quality = cov_quality)
}

# Formatting helpers ----------------------------------------------------------

fmt_pct <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f%%"), x * 100)
}

fmt_num <- function(x, digits = 2) {
  sprintf(paste0("%.", digits, "f"), x)
}

# Table helper: gt if available, else kable
make_table <- function(df, title = NULL) {
  if (requireNamespace("gt", quietly = TRUE)) {
    tbl <- gt::gt(df)
    if (!is.null(title)) tbl <- gt::tab_header(tbl, title = title)
    return(tbl)
  }
  if (requireNamespace("knitr", quietly = TRUE)) {
    return(knitr::kable(df, caption = title, digits = 3))
  }
  print(df)
}

# ===========================================================================
# STUDY 1 ANALYSIS
# ===========================================================================

analyze_study1 <- function(results) {
  cat("\n========== STUDY 1: NON-GAUSSIAN SANDWICH ==========\n\n")

  if (nrow(results) == 0) {
    cat("No Study 1 results available.\n")
    return(invisible(NULL))
  }

  # Count failures before filtering
  n_total_reps <- nrow(results)
  n_failed <- sum(results$converged == FALSE, na.rm = TRUE)
  n_na_coverage <- sum(is.na(results$coverage))
  if (n_failed > 0 || n_na_coverage > 0) {
    cat(sprintf(
      "Note: %d total rows, %d failed fits, %d NA coverage excluded.\n",
      n_total_reps,
      n_failed,
      n_na_coverage
    ))
  }

  # Filter converged
  res <- results |> filter(!is.na(coverage))

  # --- Summary table ---
  summary1 <- res |>
    group_by(family, corr_type, n, method, term_type) |>
    summarise(
      coverage = mean(coverage, na.rm = TRUE),
      se_cov = sd(coverage, na.rm = TRUE) / sqrt(n()),
      width = mean(mean_width, na.rm = TRUE),
      z_sd = mean(z_sd, na.rm = TRUE),
      rmse = mean(rmse, na.rm = TRUE),
      bias = mean(bias, na.rm = TRUE),
      fit_time = mean(fit_time, na.rm = TRUE),
      n_reps = n(),
      .groups = "drop"
    )

  # Print compact table
  cat("Coverage by family x correlation x method (ff + linear terms):\n")
  compact <- summary1 |>
    filter(term_type %in% c("ff", "linear")) |>
    select(family, corr_type, n, method, term_type, coverage, z_sd, n_reps) |>
    mutate(coverage = fmt_pct(coverage), z_sd = fmt_num(z_sd))

  make_table(compact, "Study 1: Coverage Summary")

  # --- Coverage heatmap ---
  p1_heat <- summary1 |>
    filter(term_type %in% c("ff", "linear")) |>
    ggplot(aes(
      x = method,
      y = interaction(corr_type, paste0("n=", n)),
      fill = coverage
    )) +
    geom_tile(color = "white") +
    geom_text(aes(label = fmt_pct(coverage)), size = 3) +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = NOMINAL_COVERAGE,
      limits = c(0.5, 1)
    ) +
    facet_grid(family ~ term_type) +
    geom_hline(yintercept = 0, color = "gray80") +
    labs(
      title = "Study 1: Coverage Heatmap",
      x = "Method",
      y = "Correlation x Sample Size",
      fill = "Coverage"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # --- Paired coverage difference: cluster - default ---
  paired <- summary1 |>
    filter(
      method %in% c("default", "cluster"),
      term_type %in% c("ff", "linear")
    ) |>
    select(family, corr_type, n, method, term_type, coverage) |>
    pivot_wider(names_from = method, values_from = coverage) |>
    mutate(coverage_diff = cluster - default)

  p1_diff <- paired |>
    ggplot(aes(
      x = interaction(corr_type, paste0("n=", n)),
      y = coverage_diff,
      fill = family
    )) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_hline(
      yintercept = c(-PRACTICAL_THRESHOLD, PRACTICAL_THRESHOLD),
      linetype = 3,
      color = "gray50"
    ) +
    facet_wrap(~term_type) +
    labs(
      title = "Study 1: Coverage Improvement (cluster - default)",
      x = "Correlation x n",
      y = "Coverage Difference",
      fill = "Family"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # --- z_sd plot ---
  p1_zsd <- summary1 |>
    filter(term_type %in% c("ff", "linear")) |>
    ggplot(aes(x = method, y = z_sd, color = corr_type, shape = term_type)) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_hline(yintercept = c(0.85, 1.15), linetype = 3, color = "gray50") +
    facet_grid(family ~ paste0("n=", n)) +
    labs(
      title = "Study 1: SE Calibration (z_sd, target = 1)",
      x = "Method",
      y = "z_sd"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  list(
    summary = summary1,
    plots = list(heatmap = p1_heat, diff = p1_diff, zsd = p1_zsd),
    paired = paired
  )
}

# ===========================================================================
# STUDY 2 ANALYSIS
# ===========================================================================

analyze_study2 <- function(results, cov_quality) {
  cat("\n========== STUDY 2: GRID REFINEMENT ==========\n\n")

  if (nrow(results) == 0) {
    cat("No Study 2 results available.\n")
    return(invisible(NULL))
  }

  res <- results |> filter(!is.na(coverage))

  # --- Summary table ---
  summary2 <- res |>
    group_by(corr_type, grid_label, nxgrid, nygrid, method, term_type) |>
    summarise(
      coverage = mean(coverage, na.rm = TRUE),
      se_cov = sd(coverage, na.rm = TRUE) / sqrt(n()),
      width = mean(mean_width, na.rm = TRUE),
      z_sd = mean(z_sd, na.rm = TRUE),
      rmse = mean(rmse, na.rm = TRUE),
      fit_time = mean(fit_time, na.rm = TRUE),
      n_reps = n(),
      .groups = "drop"
    )

  cat("Coverage by grid x correlation x method:\n")
  compact2 <- summary2 |>
    filter(term_type %in% c("ff", "linear", "E(Y)")) |>
    select(corr_type, grid_label, method, term_type, coverage, z_sd, n_reps) |>
    mutate(coverage = fmt_pct(coverage), z_sd = fmt_num(z_sd))

  make_table(compact2, "Study 2: Coverage by Grid Level")

  # --- Paired coverage difference: cluster - default ---
  paired2 <- summary2 |>
    filter(
      method %in% c("default", "cluster"),
      term_type %in% c("ff", "linear", "smooth", "concurrent", "E(Y)")
    ) |>
    select(corr_type, grid_label, method, term_type, coverage) |>
    pivot_wider(names_from = method, values_from = coverage) |>
    mutate(coverage_diff = cluster - default)

  p2_diff <- paired2 |>
    ggplot(aes(
      x = interaction(corr_type, grid_label),
      y = coverage_diff,
      fill = term_type
    )) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
      title = "Study 2: Coverage Improvement (cluster - default)",
      x = "Correlation x Grid",
      y = "Coverage Difference",
      fill = "Term"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # --- Grid effect: paired comparison across grid levels ---
  grid_paired <- summary2 |>
    filter(method == "cluster", term_type %in% c("ff", "linear", "E(Y)")) |>
    select(corr_type, grid_label, term_type, coverage, z_sd) |>
    pivot_wider(
      names_from = grid_label,
      values_from = c(coverage, z_sd),
      names_sep = "_"
    )

  cat("\nGrid effect on cluster coverage (fine - coarse):\n")
  print(grid_paired)

  # --- Covariance quality plots ---
  cov_plots <- NULL
  if (nrow(cov_quality) > 0) {
    cat("\nCovariance Quality Summary:\n")
    print(
      cov_quality |>
        arrange(corr_type, grid_label, method, term_type)
    )

    p2_diag <- cov_quality |>
      ggplot(aes(
        x = interaction(grid_label, method),
        y = diag_ratio_median,
        ymin = diag_ratio_iqr_low,
        ymax = diag_ratio_iqr_high,
        color = corr_type
      )) +
      geom_pointrange(position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 1, linetype = 2) +
      facet_wrap(~term_type, scales = "free_y") +
      labs(
        title = "Study 2: Diagonal Variance Ratio (model/MC)",
        subtitle = "Median with IQR. Target = 1.",
        x = "Grid x Method",
        y = "Variance Ratio",
        color = "Correlation"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    cov_plots <- list(diag_ratio = p2_diag)
  }

  list(
    summary = summary2,
    plots = c(list(diff = p2_diff), cov_plots),
    paired = paired2,
    grid_paired = grid_paired,
    cov_quality = cov_quality
  )
}

# ===========================================================================
# COMBINED MARKDOWN SUMMARY
# ===========================================================================

write_markdown_summary <- function(
  s1_analysis,
  s2_analysis,
  outpath = "ci-benchmark/study-followup-analysis.md"
) {
  lines <- character()
  add <- function(...) lines <<- c(lines, paste0(...))

  add("# Follow-up Studies: Summary\n")
  add("*Auto-generated by study-followup-analysis.R*\n")

  # Study 1
  add("## Study 1: Non-Gaussian Sandwich with ff + linear\n")

  if (!is.null(s1_analysis)) {
    s1 <- s1_analysis$summary

    add("### Key Findings\n")

    # Check: does cluster help under correlation?
    s1_ff <- s1 |>
      filter(term_type == "ff", method %in% c("default", "cluster"))

    for (fam in unique(s1_ff$family)) {
      add(sprintf("**%s family:**\n", stringr::str_to_title(fam)))

      fam_data <- s1_ff |> filter(family == fam)
      for (ct in c("iid", "ar1", "fourier_pos")) {
        ct_data <- fam_data |> filter(corr_type == ct)
        if (nrow(ct_data) == 0) next

        default_cov <- ct_data |>
          filter(method == "default") |>
          pull(coverage) |>
          mean()
        cluster_cov <- ct_data |>
          filter(method == "cluster") |>
          pull(coverage) |>
          mean()
        diff <- cluster_cov - default_cov

        flag <- if (abs(diff) > PRACTICAL_THRESHOLD) {
          if (diff > 0) " **[cluster helps]**" else " **[cluster worse]**"
        } else {
          " (negligible difference)"
        }

        add(sprintf(
          "- %s: default=%.1f%%, cluster=%.1f%% (diff=%+.1fpp)%s",
          ct,
          default_cov * 100,
          cluster_cov * 100,
          diff * 100,
          flag
        ))
      }
      add("")
    }

    # Failure rates (use raw summary which aggregated from filtered data)
    add("### Convergence\n")
    # s1 is already the aggregated summary — report n_reps info
    total_reps <- sum(s1$n_reps, na.rm = TRUE)
    add(sprintf("- Total DGP x method x term cells: %d", nrow(s1)))
    add(sprintf("- Total replicate-level observations: %d\n", total_reps))
  } else {
    add("*No Study 1 results available.*\n")
  }

  # Study 2
  add("## Study 2: Grid Refinement\n")

  if (!is.null(s2_analysis)) {
    s2 <- s2_analysis$summary

    add("### Key Findings\n")

    # Grid effect
    if (
      !is.null(s2_analysis$grid_paired) && nrow(s2_analysis$grid_paired) > 0
    ) {
      add("**Grid effect on cluster coverage:**\n")
      gp <- s2_analysis$grid_paired
      # Find coverage columns dynamically
      cov_cols <- grep("^coverage_", names(gp), value = TRUE)
      if (length(cov_cols) >= 2) {
        add("| Correlation | Term | Coarse | Fine |")
        add("|---|---|---|---|")
        for (i in seq_len(nrow(gp))) {
          row <- gp[i, ]
          add(sprintf(
            "| %s | %s | %s | %s |",
            row$corr_type,
            row$term_type,
            fmt_pct(row[[cov_cols[1]]]),
            fmt_pct(row[[cov_cols[2]]])
          ))
        }
        add("")
      }
    }

    # Covariance quality
    if (nrow(s2_analysis$cov_quality) > 0) {
      add("### Covariance Quality (diagonal variance ratio)\n")
      add("Target: median ratio = 1.0 (model SE matches MC SD)\n")

      cq <- s2_analysis$cov_quality |>
        filter(term_type %in% c("linear", "smooth", "concurrent"))

      if (nrow(cq) > 0) {
        add("| Corr | Grid | Method | Term | Median Ratio | IQR |")
        add("|---|---|---|---|---|---|")
        for (i in seq_len(nrow(cq))) {
          r <- cq[i, ]
          add(sprintf(
            "| %s | %s | %s | %s | %.2f | [%.2f, %.2f] |",
            r$corr_type,
            r$grid_label,
            r$method,
            r$term_type,
            r$diag_ratio_median,
            r$diag_ratio_iqr_low,
            r$diag_ratio_iqr_high
          ))
        }
        add("")
      }
    }
  } else {
    add("*No Study 2 results available.*\n")
  }

  # Practical significance note
  add("## Methodology Notes\n")
  add(sprintf(
    "- Practical significance threshold: %d percentage points",
    PRACTICAL_THRESHOLD * 100
  ))
  add("- z_sd target: 1.0 (meaningful deviation if |z_sd - 1| > 0.15)")
  add("- MC SEs at 150 reps for 90%% coverage: ~2.4 percentage points")
  add(
    "- Binomial DGP uses Gaussian copula to preserve marginal regression coefficients"
  )
  add(
    "- Poisson DGP uses additive eps on linear predictor (collapsible log link)"
  )
  add(
    "- Grid comparison uses paired design (fine grid + deterministic subsampling)\n"
  )

  writeLines(lines, outpath)
  cat("Markdown summary written to", outpath, "\n")
  invisible(lines)
}

# ===========================================================================
# SANITY CHECKS
# ===========================================================================

run_sanity_checks <- function(s1_results, s2_results) {
  cat("\n========== SANITY CHECKS ==========\n\n")
  all_pass <- TRUE

  # Check 1: Under strong AR1, cluster > default by >10pp
  if (nrow(s1_results) > 0) {
    ar1_check <- s1_results |>
      filter(
        !is.na(coverage),
        corr_type == "ar1",
        term_type %in% c("ff", "linear"),
        method %in% c("default", "cluster")
      ) |>
      group_by(method) |>
      summarise(cov = mean(coverage, na.rm = TRUE), .groups = "drop") |>
      pivot_wider(names_from = method, values_from = cov)

    if (
      nrow(ar1_check) > 0 &&
        "cluster" %in% names(ar1_check) &&
        "default" %in% names(ar1_check) &&
        !is.na(ar1_check$cluster[1]) &&
        !is.na(ar1_check$default[1])
    ) {
      diff <- ar1_check$cluster[1] - ar1_check$default[1]
      pass <- diff > 0.10
      cat(sprintf(
        "Study 1 AR1: cluster - default = %.1fpp [%s]\n",
        diff * 100,
        if (pass) "PASS" else "FAIL"
      ))
      if (!pass) all_pass <- FALSE
    }
  }

  # Check 2: HC should not match cluster improvements under correlation
  if (nrow(s1_results) > 0) {
    hc_check <- s1_results |>
      filter(
        !is.na(coverage),
        corr_type == "ar1",
        term_type %in% c("ff", "linear"),
        method %in% c("hc", "cluster")
      ) |>
      group_by(method) |>
      summarise(cov = mean(coverage, na.rm = TRUE), .groups = "drop") |>
      pivot_wider(names_from = method, values_from = cov)

    if (
      nrow(hc_check) > 0 &&
        "cluster" %in% names(hc_check) &&
        "hc" %in% names(hc_check) &&
        !is.na(hc_check$cluster[1]) &&
        !is.na(hc_check$hc[1])
    ) {
      diff <- hc_check$cluster[1] - hc_check$hc[1]
      cat(sprintf(
        "Study 1 AR1: cluster - hc = %.1fpp (cluster should beat hc)\n",
        diff * 100
      ))
    }
  }

  # Check 3: IID default should have good coverage (~90%)
  if (nrow(s1_results) > 0) {
    iid_check <- s1_results |>
      filter(
        !is.na(coverage),
        corr_type == "iid",
        method == "default",
        term_type == "linear"
      ) |>
      summarise(cov = mean(coverage, na.rm = TRUE))

    if (nrow(iid_check) > 0) {
      pass <- abs(iid_check$cov - 0.90) < 0.10
      cat(sprintf(
        "Study 1 IID default coverage: %.1f%% [%s]\n",
        iid_check$cov * 100,
        if (pass) "PASS" else "FAIL"
      ))
      if (!pass) all_pass <- FALSE
    }
  }

  # Study 2 checks
  if (nrow(s2_results) > 0) {
    s2_ar1 <- s2_results |>
      filter(
        !is.na(coverage),
        corr_type == "ar1",
        method %in% c("default", "cluster"),
        term_type %in% c("ff", "linear")
      ) |>
      group_by(method) |>
      summarise(cov = mean(coverage, na.rm = TRUE), .groups = "drop") |>
      pivot_wider(names_from = method, values_from = cov)

    if (
      nrow(s2_ar1) > 0 &&
        "cluster" %in% names(s2_ar1) &&
        "default" %in% names(s2_ar1) &&
        !is.na(s2_ar1$cluster[1]) &&
        !is.na(s2_ar1$default[1])
    ) {
      diff <- s2_ar1$cluster[1] - s2_ar1$default[1]
      pass <- diff > 0.10
      cat(sprintf(
        "Study 2 AR1: cluster - default = %.1fpp [%s]\n",
        diff * 100,
        if (pass) "PASS" else "FAIL"
      ))
      if (!pass) all_pass <- FALSE
    }
  }

  cat(sprintf(
    "\nOverall sanity: %s\n",
    if (all_pass) "ALL PASS" else "SOME FAILED"
  ))
  all_pass
}

# ===========================================================================
# SAVE PLOTS
# ===========================================================================

save_analysis_plots <- function(
  s1_analysis,
  s2_analysis,
  outdir = "ci-benchmark"
) {
  if (!is.null(s1_analysis$plots)) {
    for (name in names(s1_analysis$plots)) {
      p <- s1_analysis$plots[[name]]
      if (!is.null(p)) {
        path <- file.path(outdir, paste0("study1_", name, ".pdf"))
        ggsave(path, p, width = 10, height = 7)
        cat("Saved", path, "\n")
      }
    }
  }

  if (!is.null(s2_analysis$plots)) {
    for (name in names(s2_analysis$plots)) {
      p <- s2_analysis$plots[[name]]
      if (!is.null(p)) {
        path <- file.path(outdir, paste0("study2_", name, ".pdf"))
        ggsave(path, p, width = 10, height = 7)
        cat("Saved", path, "\n")
      }
    }
  }
}

# ===========================================================================
# MAIN
# ===========================================================================

if (sys.nframe() == 0) {
  cat("Follow-up Studies Analysis\n")
  cat("==========================\n\n")

  # Load results
  s1_results <- tryCatch(load_study1(), error = function(e) {
    cat("Study 1 results not found:", e$message, "\n")
    tibble()
  })
  s2_data <- tryCatch(load_study2(), error = function(e) {
    cat("Study 2 results not found:", e$message, "\n")
    list(results = tibble(), cov_quality = tibble())
  })

  # Analyze
  s1_analysis <- if (nrow(s1_results) > 0) analyze_study1(s1_results) else NULL
  s2_analysis <- if (nrow(s2_data$results) > 0) {
    analyze_study2(s2_data$results, s2_data$cov_quality)
  } else {
    NULL
  }

  # Sanity checks
  run_sanity_checks(s1_results, s2_data$results)

  # Save outputs
  save_analysis_plots(s1_analysis, s2_analysis)
  write_markdown_summary(s1_analysis, s2_analysis)

  cat("\nAnalysis complete.\n")
}
