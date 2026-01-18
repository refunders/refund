# viz-functions.R
# Visualization functions for confint-benchmark results
# Moved from confint-benchmark.R to reduce main file size

plot_summary <- function(summary_df, nominal = 0.90) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }
  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = method, y = coverage, color = fit_family)
  ) +
    ggplot2::geom_hline(yintercept = nominal, linetype = 2, color = "gray40") +
    ggplot2::geom_hline(yintercept = 0.90, linetype = 3, color = "gray60") +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.4),
      size = 2
    ) +
    ggplot2::facet_grid(
      term_type ~ error_dist + corr_type + hetero_type + snr,
      labeller = ggplot2::label_both
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    ggplot2::labs(
      title = "CI Coverage by Term Type and DGP Settings",
      y = "Empirical Coverage",
      x = "Method"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(size = 8)
    )
}

plot_metric <- function(
  summary_df,
  metric = c("coverage", "width", "rmse"),
  nominal = 0.90
) {
  metric <- match.arg(metric)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = method, y = .data[[metric]], color = fit_family)
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.4),
      size = 2
    ) +
    ggplot2::facet_grid(
      term_type ~ error_dist + corr_type + hetero_type + snr,
      labeller = ggplot2::label_both
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(size = 8)
    )

  if (metric == "coverage") {
    p <- p +
      ggplot2::geom_hline(
        yintercept = nominal,
        linetype = 2,
        color = "gray40"
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(y = "Coverage", title = "CI Coverage")
  } else if (metric == "width") {
    p <- p + ggplot2::labs(y = "Mean CI Width", title = "CI Width")
  } else {
    p <- p + ggplot2::labs(y = "RMSE", title = "Estimation Accuracy")
  }
  p
}

plot_surface_metric <- function(
  scored,
  term_type = c("ff", "smooth"),
  metric = "covered",
  show_delta = FALSE,
  nominal = 0.90
) {
  term_type <- match.arg(term_type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  df <- scored[scored$term_type == term_type, , drop = FALSE]
  if (nrow(df) == 0) {
    message("No data for term_type: ", term_type)
    return(NULL)
  }

  agg <- aggregate(
    df[[metric]],
    by = list(x = df$x, y = df$y, method = df$method),
    FUN = mean
  )
  names(agg)[4] <- metric

  if (show_delta && metric == "covered") {
    agg$fill_val <- agg[[metric]] - nominal
    fill_label <- paste0("Coverage - ", nominal)
    fill_limits <- c(-0.5, 0.5)
    fill_colors <- c("darkred", "white", "darkblue")
  } else {
    agg$fill_val <- agg[[metric]]
    fill_label <- if (metric == "covered") "Coverage" else metric
    fill_limits <- if (metric == "covered") c(0, 1) else NULL
    fill_colors <- NULL
  }

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = x, y = y, fill = fill_val)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~method) +
    ggplot2::labs(
      title = sprintf("%s: %s by Grid Position", term_type, fill_label),
      x = if (term_type == "ff") "s" else "z",
      y = "t",
      fill = fill_label
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_fixed()

  if (!is.null(fill_colors)) {
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = fill_colors[1],
        mid = fill_colors[2],
        high = fill_colors[3],
        midpoint = 0,
        limits = fill_limits
      )
  } else if (!is.null(fill_limits)) {
    p <- p + ggplot2::scale_fill_viridis_c(limits = fill_limits)
  } else {
    p <- p + ggplot2::scale_fill_viridis_c()
  }
  p
}

plot_curve_metric <- function(
  scored,
  term_type = c("concurrent", "linear"),
  metric = "covered",
  show_ribbon = TRUE,
  nominal = 0.90
) {
  term_type <- match.arg(term_type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  df <- scored[scored$term_type == term_type, , drop = FALSE]
  if (nrow(df) == 0) {
    message("No data for term_type: ", term_type)
    return(NULL)
  }

  agg <- aggregate(
    df[[metric]],
    by = list(x = df$x, method = df$method),
    FUN = function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v)))
  )
  agg <- cbind(agg[, 1:2], as.data.frame(agg[[3]]))
  names(agg)[3:4] <- c("mean_val", "se_val")

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = x, y = mean_val, color = method))

  if (show_ribbon) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = mean_val - 1.96 * se_val,
          ymax = mean_val + 1.96 * se_val,
          fill = method
        ),
        alpha = 0.2,
        color = NA
      )
  }

  p <- p +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(
      title = sprintf("%s: %s vs t", term_type, metric),
      x = "t",
      y = if (metric == "covered") "Coverage" else metric
    ) +
    ggplot2::theme_bw()

  if (metric == "covered") {
    p <- p +
      ggplot2::geom_hline(
        yintercept = nominal,
        linetype = 2,
        color = "gray40"
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1))
  }
  p
}

plot_zscore_dist <- function(scored, facet_by = "term_type") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  scored$z_score <- (scored$est - scored$truth) / scored$se

  p <- ggplot2::ggplot(scored, ggplot2::aes(x = z_score)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 50,
      fill = "steelblue",
      alpha = 0.7
    ) +
    ggplot2::stat_function(fun = dnorm, color = "red", linewidth = 1) +
    ggplot2::facet_wrap(as.formula(paste("~", facet_by))) +
    ggplot2::labs(
      title = "Z-score Distribution vs Standard Normal",
      x = "Z-score = (estimate - truth) / SE",
      y = "Density"
    ) +
    ggplot2::theme_bw() +
    ggplot2::xlim(-5, 5)
  p
}

plot_estimate_vs_truth <- function(
  scored,
  term_type,
  method = NULL,
  rep_id = 1
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  df <- scored[scored$term_type == term_type, , drop = FALSE]
  if (!is.null(method)) df <- df[df$method == method, , drop = FALSE]
  df <- df[df$rep_id == rep_id, , drop = FALSE]

  if (nrow(df) == 0) {
    message("No data for specified filters")
    return(NULL)
  }

  is_2d <- !all(is.na(df$y))

  if (is_2d) {
    df_long <- rbind(
      data.frame(x = df$x, y = df$y, value = df$truth, type = "Truth"),
      data.frame(x = df$x, y = df$y, value = df$est, type = "Estimate")
    )
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~type) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = sprintf("%s: Estimate vs Truth", term_type)) +
      ggplot2::theme_bw()
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        fill = "steelblue",
        alpha = 0.3
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = est, color = "Estimate"),
        linewidth = 1
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = truth, color = "Truth"),
        linewidth = 1,
        linetype = 2
      ) +
      ggplot2::scale_color_manual(
        values = c("Estimate" = "steelblue", "Truth" = "red")
      ) +
      ggplot2::labs(
        title = sprintf("%s: Estimate vs Truth (rep %d)", term_type, rep_id),
        x = "t",
        y = "Coefficient",
        color = ""
      ) +
      ggplot2::theme_bw()
  }
  p
}

plot_coverage_vs_snr <- function(summary_df, nominal = 0.90) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = factor(snr),
      y = coverage,
      color = method,
      shape = fit_family
    )
  ) +
    ggplot2::geom_hline(yintercept = nominal, linetype = 2, color = "gray40") +
    ggplot2::geom_hline(yintercept = 0.90, linetype = 3, color = "gray60") +
    ggplot2::geom_point(
      size = 3,
      position = ggplot2::position_dodge(width = 0.3)
    ) +
    ggplot2::facet_wrap(~term_type) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Coverage vs SNR by Term Type",
      x = "SNR",
      y = "Coverage"
    ) +
    ggplot2::theme_bw()
}

plot_dgp_summary <- function(scored, summary_df, dgp_label = "") {
  if (
    !requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)
  ) {
    message("ggplot2 and/or patchwork not installed; skipping plots.")
    return(NULL)
  }

  p1 <- plot_metric(summary_df, "coverage") +
    ggplot2::ggtitle(paste("Coverage", dgp_label))
  p2 <- plot_zscore_dist(scored) +
    ggplot2::ggtitle("Z-score Distribution")

  patchwork::wrap_plots(p1, p2, ncol = 1)
}
