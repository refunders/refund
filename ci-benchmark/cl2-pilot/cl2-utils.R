# cl2-utils.R
#
# Shared utilities for CL2-style cluster-robust covariance experiments.
# Supports standard single-linear-predictor GLM families (e.g., gaussian,
# poisson, binomial, Gamma) and gaulss for mgcv::gam/pffr fits.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

cl2_build_cluster_id <- function(pffr_meta) {
  if (isTRUE(pffr_meta$is_sparse)) {
    cluster_id <- pffr_meta$ydata$.obs
  } else {
    cluster_id <- rep(seq_len(pffr_meta$nobs), each = pffr_meta$nyindex)
  }
  if (!is.null(pffr_meta$missing_indices)) {
    cluster_id <- cluster_id[-pffr_meta$missing_indices]
  }
  cluster_id
}

cl2_sym_inv_sqrt <- function(M, tol = 1e-8) {
  M <- 0.5 * (M + t(M))
  ee <- eigen(M, symmetric = TRUE)
  vals <- pmax(ee$values, tol)
  ee$vectors %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(ee$vectors)
}

cl2_build_working_standard <- function(b, cluster_id) {
  X <- model.matrix(b)
  y <- as.vector(b$y)
  mu <- as.vector(b$fitted.values)
  eta <- as.vector(b$linear.predictors)

  var_mu <- as.vector(b$family$variance(mu))
  mu_eta <- as.vector(b$family$mu.eta(eta))
  sig2 <- b$sig2 %||% 1
  if (!is.finite(sig2) || sig2 <= 0) sig2 <- 1

  pw <- b$prior.weights
  if (is.null(pw)) pw <- rep(1, length(y))
  pw <- as.vector(pw)

  denom <- sig2 * var_mu
  bad <- !is.finite(denom) | denom <= 0 | !is.finite(mu_eta)
  denom[bad] <- NA_real_

  # Score factorization:
  #   U = X' * [pw * mu_eta * (y - mu) / (sig2 * V(mu))]
  #     = (X * sqrt_w)' * z
  sqrt_common <- sqrt(pw / denom)
  x_scale <- mu_eta * sqrt_common
  z <- (y - mu) * sqrt_common

  x_scale[!is.finite(x_scale)] <- 0
  z[!is.finite(z)] <- 0

  Xw <- X * x_scale

  list(Xw = Xw, z = z, cluster_id = cluster_id)
}

cl2_build_working_gaulss <- function(b, cluster_id) {
  X <- model.matrix(b)
  lpi <- attr(X, "lpi")
  if (is.null(lpi) || length(lpi) < 2) {
    stop(
      "gaulss fit is missing 'lpi' structure on model matrix.",
      call. = FALSE
    )
  }

  n_obs <- length(b$y)
  if (length(cluster_id) != n_obs) {
    stop(
      "cluster_id length does not match gaulss observation count.",
      call. = FALSE
    )
  }

  mu <- b$fitted.values[seq_len(n_obs)]
  tau <- b$fitted.values[n_obs + seq_len(n_obs)]
  eta1 <- b$linear.predictors[seq_len(n_obs)]
  eta2 <- b$linear.predictors[n_obs + seq_len(n_obs)]
  r <- b$y - mu

  # Same score components used in compute_gaulss_scores() in R/pffr-core.R.
  w1 <- (tau^2 * r) * b$family$linfo[[1]]$mu.eta(eta1)
  w2 <- (1 / tau - tau * r^2) * b$family$linfo[[2]]$mu.eta(eta2)

  pw <- b$prior.weights
  if (!is.null(pw) && any(pw != 1)) {
    w1 <- pw * w1
    w2 <- pw * w2
  }

  make_block <- function(w, lp_idx) {
    abs_sqrt <- sqrt(abs(w))
    z <- sign(w) * abs_sqrt
    z[!is.finite(z)] <- 0

    Xw_block <- matrix(0, nrow = n_obs, ncol = ncol(X))
    if (!is.null(lp_idx) && length(lp_idx)) {
      Xw_block[, lp_idx] <- X[, lp_idx, drop = FALSE] * abs_sqrt
      Xw_block[!is.finite(Xw_block)] <- 0
    }
    list(Xw = Xw_block, z = z)
  }

  block1 <- make_block(w1, lpi[[1]])
  block2 <- make_block(w2, lpi[[2]])

  list(
    Xw = rbind(block1$Xw, block2$Xw),
    z = c(block1$z, block2$z),
    cluster_id = rep(cluster_id, times = 2L)
  )
}

cl2_build_working_representation <- function(b, cluster_id) {
  fam <- tolower(as.character(b$family$family %||% ""))

  if (fam == "gaulss") {
    return(cl2_build_working_gaulss(b, cluster_id))
  }

  if (!is.null(b$family$sandwich)) {
    stop(
      "CL2 pilot path is not implemented for family '",
      b$family$family,
      "' (custom family$sandwich).",
      call. = FALSE
    )
  }

  cl2_build_working_standard(b, cluster_id)
}

#' Compute CL2-style sandwich covariance for GAM (standard GLM families)
#'
#' @param b Fitted gam object (class "pffr" stripped).
#' @param cluster_id Cluster membership vector for vectorized observations.
#' @param freq If TRUE, use frequentist sandwich (B2=0); else Bayesian (Vp-Ve).
#' @param tol Eigenvalue floor for numerical stability.
#' @param leverage_cap Cap on cluster leverage eigenvalues (<1).
#' @returns Covariance matrix with attribute n_capped_clusters.
cl2_gam_sandwich_cluster <- function(
  b,
  cluster_id,
  freq = FALSE,
  tol = 1e-8,
  leverage_cap = 0.999
) {
  work <- cl2_build_working_representation(b, cluster_id)
  Xw <- work$Xw
  z <- work$z
  cluster_id_work <- work$cluster_id
  Vp <- b$Vp
  B2 <- if (freq) 0 else (b$Vp - b$Ve)

  groups <- unique(cluster_id_work)
  G <- length(groups)
  if (G < 2) {
    stop("Need at least two clusters for cluster sandwich.", call. = FALSE)
  }

  p <- ncol(Xw)
  meat <- matrix(0, nrow = p, ncol = p)
  n_capped_clusters <- 0L

  for (g in groups) {
    idx <- which(cluster_id_work == g)
    Xwg <- Xw[idx, , drop = FALSE]
    zg <- z[idx]

    Hgg <- Xwg %*% Vp %*% t(Xwg)
    Hgg <- 0.5 * (Hgg + t(Hgg))

    ee_H <- eigen(Hgg, symmetric = TRUE)
    if (any(ee_H$values > leverage_cap, na.rm = TRUE)) {
      n_capped_clusters <- n_capped_clusters + 1L
      ee_H$values <- pmin(ee_H$values, leverage_cap)
      Hgg <- ee_H$vectors %*%
        diag(ee_H$values, nrow = length(ee_H$values)) %*%
        t(ee_H$vectors)
    }

    Mg <- diag(length(idx)) - Hgg
    Ag <- cl2_sym_inv_sqrt(Mg, tol = tol)

    Ug <- crossprod(Xwg, Ag %*% zg)
    meat <- meat + Ug %*% t(Ug)
  }

  hc1 <- G / (G - 1)
  V <- hc1 * Vp %*% meat %*% Vp + B2
  V <- 0.5 * (V + t(V))
  attr(V, "n_capped_clusters") <- n_capped_clusters
  V
}

#' Apply CL2 covariance to a pffr fit
#'
#' @param fit Fitted pffr model.
#' @param tol Eigenvalue floor.
#' @param leverage_cap Leverage cap (<1).
#' @returns List with modified fit and number of capped clusters.
cl2_apply_to_fit <- function(fit, tol = 1e-8, leverage_cap = 0.999) {
  algorithm <- fit$pffr$algorithm %||% "gam"
  use_gam_slot <- as.character(algorithm) %in% c("gamm", "gamm4")

  gam_obj <- if (use_gam_slot) fit$gam else fit
  gam_obj_stripped <- gam_obj
  class(gam_obj_stripped) <- setdiff(class(gam_obj_stripped), "pffr")

  cluster_id <- cl2_build_cluster_id(gam_obj$pffr)

  Vp_cl2 <- cl2_gam_sandwich_cluster(
    gam_obj_stripped,
    cluster_id = cluster_id,
    freq = FALSE,
    tol = tol,
    leverage_cap = leverage_cap
  )
  Ve_cl2 <- cl2_gam_sandwich_cluster(
    gam_obj_stripped,
    cluster_id = cluster_id,
    freq = TRUE,
    tol = tol,
    leverage_cap = leverage_cap
  )

  gam_obj$Vp <- gam_obj$Vc <- Vp_cl2
  gam_obj$Ve <- Ve_cl2

  fit_cl2 <- fit
  if (use_gam_slot) {
    fit_cl2$gam <- gam_obj
  } else {
    fit_cl2 <- gam_obj
  }

  list(
    fit = fit_cl2,
    n_capped_clusters = attr(Vp_cl2, "n_capped_clusters") %||% 0L
  )
}

#' Compute benchmark metrics for CL2-adjusted fit
#'
#' @param fit Fitted pffr model.
#' @param sim Simulation object expected by compute_all_metrics().
#' @param alpha Significance level.
#' @param method_name Method label to assign in output.
#' @param tol Eigenvalue floor.
#' @param leverage_cap Leverage cap (<1).
#' @returns List with metrics tibble and n_capped_clusters.
cl2_compute_metrics <- function(
  fit,
  sim,
  alpha = 0.10,
  method_name = "pffr_cl2",
  tol = 1e-8,
  leverage_cap = 0.999
) {
  fit_cl2 <- cl2_apply_to_fit(fit, tol = tol, leverage_cap = leverage_cap)

  metrics <- compute_all_metrics(
    fit = fit_cl2$fit,
    sim = sim,
    method = method_name,
    use_sandwich = "none",
    alpha = alpha,
    err_struct = sim$err_struct
  )

  list(
    metrics = metrics,
    n_capped_clusters = fit_cl2$n_capped_clusters
  )
}
