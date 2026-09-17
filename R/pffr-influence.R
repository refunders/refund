# Fixed-fit cluster inference core, GPL (>= 2), as in refund.
# Does not propagate smoothing selection or remove smoothing bias.

#' Compressed fixed-fit cluster influence geometry
#'
#' Uses Z_g = Q_g T_g and R_g = I - Z_g C Z_g', C = 2B - B Z'Z B.
#' Numerical rank is measured; no universal response-basis rank bound is used.
#' @param Xw Finite working-likelihood-scaled design.
#' @param Vp Penalized bread in the same scaling.
#' @param cluster_id One nonmissing cluster label per working row.
#' @param z Optional finite working residuals.
#' @param adjustment Exact CL2, historical shortcut, or none.
#' @param leverage_cap Numerical cap strictly between zero and one.
#' @param tol Positive shortcut eigenvalue floor.
#' @param rank_tol Relative SVD tolerance; NULL uses dimension times machine epsilon.
#' @returns Influence object with separate sampling columns, B2 and geometry.
#'   Its `diagnostics` data frame carries the per-cluster hat-invariant monitors
#'   `max_leverage`/`min_hat_eig` (extreme eigenvalues of \eqn{H_{gg}}),
#'   `max_obs_leverage`/`min_obs_leverage` (extreme \eqn{h_{ii}}) and
#'   `min_block_eig_rel`.
#' @keywords internal
pffr_influence_core <- function(
  Xw,
  Vp,
  cluster_id,
  z = NULL,
  adjustment = c("exact", "shortcut", "none"),
  leverage_cap = .999,
  tol = 1e-8,
  rank_tol = NULL
) {
  adjustment <- match.arg(adjustment)
  Xw <- as.matrix(Xw)
  Vp <- as.matrix(Vp)
  if (
    !is.numeric(Xw) ||
      !is.numeric(Vp) ||
      any(!is.finite(Xw)) ||
      any(!is.finite(Vp)) ||
      min(dim(Xw)) < 1L ||
      !identical(dim(Vp), rep(ncol(Xw), 2L))
  )
    stop(
      "Finite conformable numeric Xw and Vp matrices are required.",
      call. = FALSE
    )
  if (
    !is.atomic(cluster_id) ||
      !is.null(dim(cluster_id)) ||
      length(cluster_id) != nrow(Xw) ||
      anyNA(cluster_id)
  )
    stop(
      "cluster_id must identify every working row without missing values.",
      call. = FALSE
    )
  if (
    !is.null(z) &&
      (!is.numeric(z) || length(z) != nrow(Xw) || any(!is.finite(z)))
  )
    stop("z must contain one finite residual per working row.", call. = FALSE)
  if (
    length(leverage_cap) != 1L ||
      !is.finite(leverage_cap) ||
      leverage_cap <= 0 ||
      leverage_cap >= 1 ||
      length(tol) != 1L ||
      !is.finite(tol) ||
      tol <= 0
  )
    stop(
      "leverage_cap must be in (0,1) and tol must be positive.",
      call. = FALSE
    )
  if (
    !is.null(rank_tol) &&
      (length(rank_tol) != 1L ||
        !is.finite(rank_tol) ||
        rank_tol < 0 ||
        rank_tol >= 1)
  )
    stop("rank_tol must be in [0,1).", call. = FALSE)
  groups <- unique(cluster_id)
  G <- length(groups)
  if (G < 2L)
    stop("Cluster inference requires at least two clusters.", call. = FALSE)
  B <- (Vp + t(Vp)) / 2
  C <- 2 * B - crossprod(B, crossprod(Xw) %*% B)
  C <- (C + t(C)) / 2
  p <- ncol(Xw)
  blocks <- vector("list", G)
  K <- if (is.null(z)) NULL else matrix(0, p, G)
  rank <- size <- n_floored <- integer(G)
  min_eig <- max_kappa <- max_hat <- discarded <- numeric(G)
  # Hat-invariant monitors (study LB, claim P-LB5): largest and smallest
  # per-observation leverage, the extreme eigenvalues of H_gg, and the smallest
  # residual-block eigenvalue relative to that block's own scale. All are
  # by-products of geometry computed anyway. The lower bounds matter because a
  # numerically inconsistent bread can make the penalized hat *indefinite*
  # (h_ii < 0 or eigen(H_gg) < 0), which the upper-bound monitors never see.
  max_obs <- min_obs <- min_hat <- min_eig_rel <- numeric(G)
  for (g in seq_len(G)) {
    idx <- which(cluster_id == groups[g])
    Zg <- Xw[idx, , drop = FALSE]
    size[g] <- length(idx)
    dec <- svd(Zg, nu = min(dim(Zg)), nv = min(dim(Zg)))
    threshold <- (rank_tol %||% (max(dim(Zg)) * .Machine$double.eps)) *
      max(dec$d)
    keep <- which(dec$d > threshold)
    rank[g] <- length(keep)
    dropped <- setdiff(seq_along(dec$d), keep)
    discarded[g] <- if (length(dropped)) max(dec$d[dropped]) else 0
    if (!length(keep)) {
      # A rank-0 cluster contributes nothing: H_gg = 0 and the residual block is
      # the identity, so every monitor takes its unproblematic value.
      blocks[[g]] <- list(T = matrix(0, 0L, p), A = matrix(0, 0L, 0L))
      min_eig[g] <- max_kappa[g] <- min_eig_rel[g] <- 1
      next
    }
    T <- dec$d[keep] * t(dec$v[, keep, drop = FALSE])
    Hsmall <- T %*% B %*% t(T)
    Hsmall <- (Hsmall + t(Hsmall)) / 2
    hat_eig <- eigen(Hsmall, symmetric = TRUE, only.values = TRUE)$values
    max_hat[g] <- max(hat_eig)
    min_hat[g] <- min(hat_eig)
    # H_gg = U Hsmall U' with orthonormal U, so h_ii needs no dense hat block.
    Ug <- dec$u[, keep, drop = FALSE]
    h_ii <- rowSums((Ug %*% Hsmall) * Ug)
    max_obs[g] <- max(h_ii)
    min_obs[g] <- min(h_ii)
    if (adjustment == "none") {
      A <- diag(length(keep))
      min_eig[g] <- max_kappa[g] <- 1
    } else {
      small <- diag(length(keep)) -
        T %*% (if (adjustment == "exact") C else B) %*% t(T)
      ee <- eigen((small + t(small)) / 2, symmetric = TRUE)
      values <- if (length(keep) < length(idx)) c(ee$values, 1) else ee$values
      min_eig[g] <- min(values)
      # Upstream's definition |max v| / |min v| (NOT max|v| / min|v|): the two
      # differ exactly when a negative eigenvalue is present, i.e. in the
      # P-LB5 indefinite-bread case this monitor exists to expose.
      max_kappa[g] <- abs(max(values)) / max(abs(min(values)), 1e-300)
      min_eig_rel[g] <- min(values) / max(abs(max(values)), 1e-300)
      floor <- if (adjustment == "exact") (1 - leverage_cap)^2 else
        max(1 - leverage_cap, tol)
      n_floored[g] <- sum(ee$values < floor)
      A <- tcrossprod(
        sweep(ee$vectors, 2L, sqrt(pmax(ee$values, floor)), "/"),
        ee$vectors
      )
    }
    blocks[[g]] <- list(T = T, A = A)
    if (!is.null(z))
      K[, g] <- B %*%
        crossprod(T, A %*% crossprod(dec$u[, keep, drop = FALSE], z[idx]))
  }
  structure(
    list(
      B = B,
      C = C,
      K = K,
      blocks = blocks,
      G = G,
      groups = groups,
      adjustment = adjustment,
      correction = G / (G - 1),
      B2 = matrix(0, p, p),
      diagnostics = data.frame(
        cluster = as.character(groups),
        size = size,
        rank = rank,
        n_floored = n_floored,
        min_block_eig = min_eig,
        max_block_kappa = max_kappa,
        max_leverage = max_hat,
        min_hat_eig = min_hat,
        max_obs_leverage = max_obs,
        min_obs_leverage = min_obs,
        min_block_eig_rel = min_eig_rel,
        max_discarded_singular_value = discarded
      ),
      leverage_cap = leverage_cap,
      tol = tol,
      rank_tol = rank_tol,
      version = "fixed-fit-core-2026-09-09"
    ),
    class = "pffr_influence"
  )
}

#' Assemble covariance from fixed-fit influence columns
#' @param core Influence object with residual columns K.
#' @param freq TRUE excludes B2.
#' @param b2 Include B2 when freq is FALSE.
#' @param center_scores Centre influence columns before forming covariance.
#' @returns Covariance with raw numerical diagnostics.
#' @keywords internal
pffr_influence_vcov <- function(
  core,
  freq = FALSE,
  b2 = TRUE,
  center_scores = FALSE
) {
  if (is.null(core$K))
    stop("Residual influence columns were not computed.", call. = FALSE)
  K <- core$K
  if (isTRUE(center_scores)) K <- K - rowMeans(K)
  V <- core$correction * tcrossprod(K)
  if (!freq && isTRUE(b2)) V <- V + core$B2
  V <- (V + t(V)) / 2
  d <- core$diagnostics
  attr(V, "cl2_adjustment") <- core$adjustment
  attr(V, "n_capped_clusters") <- if (core$adjustment == "shortcut")
    sum(d$n_floored > 0L) else 0L
  attr(V, "n_adjusted") <- if (core$adjustment == "exact")
    sum(d$n_floored > 0L) else 0L
  attr(V, "max_leverage") <- max(d$max_leverage)
  attr(V, "min_block_eig") <- min(d$min_block_eig)
  attr(V, "max_block_kappa") <- max(d$max_block_kappa)
  attr(V, "cluster_rank") <- d$rank
  attr(V, "inference_core_version") <- core$version
  # Study-LB P-LB5 hat-invariant monitors. The upper bound on the residual
  # block is only meaningful for the exact block and the upper bound on
  # eigen(H_gg) only for the shortcut, but the lower bounds (h_ii >= 0,
  # eigen(H_gg) >= 0) and the h_ii upper bound hold on every path, so those are
  # checked regardless of the adjustment.
  attr(V, "max_obs_leverage") <- max(d$max_obs_leverage)
  attr(V, "min_obs_leverage") <- min(d$min_obs_leverage)
  attr(V, "min_hat_eig") <- min(d$min_hat_eig)
  attr(V, "min_block_eig_rel") <- min(d$min_block_eig_rel)
  attr(V, "hat_invariant_violation") <- pffr_hat_invariant_violation(
    max_obs_leverage = max(d$max_obs_leverage),
    max_leverage = if (core$adjustment == "shortcut") max(d$max_leverage) else
      NA_real_,
    min_block_eig_rel = if (core$adjustment == "exact")
      min(d$min_block_eig_rel) else 0,
    min_obs_leverage = min(d$min_obs_leverage),
    min_hat_eig = min(d$min_hat_eig)
  )
  V
}

#' Central Gaussian sampling-variance moment degrees of freedom
#'
#' Gamma_gh = 1(g=h)||q_g||^2 - t_g' C t_h retains fit residualization.
#' This is conditional on weights and smoothing parameters, and is not an
#' exact t law. Noncentral means, B2 and smoothing selection are not covered.
#'
#' `df_gram = "diagonal"` drops the off-diagonal residualization and evaluates
#' the working-iid shortcut `(sum_g ||q_g||^2)^2 / sum_g ||q_g||^4` from **the
#' same `q_g` as the covariance**, i.e. with the resolved `adjustment` of this
#' influence object. It is retained only for re-scoring comparisons against
#' historical Satterthwaite results; it returns about `G` where the residualized
#' moment df returns `G - 1`.
#'
#' The historical (pre-2026-09) df additionally *always* used the shortcut
#' leverage weight \eqn{A_g = (I - H_{gg})^{-1/2}}, whatever covariance was
#' requested. `"diagonal"` therefore reproduces the historical df **exactly only
#' in combination with `cl2_adjustment = "shortcut"`**; on an exact-CL2 fit it
#' is the diagonal moment df of the exact geometry and differs from the
#' historical number (2e-5 to 2e-2 relative on the package fixtures; there is no
#' bound on the difference in general).
#' @param core Fixed-fit influence object.
#' @param Xp Finite full-coefficient contrasts, one per row.
#' @param chunk_size Positive number of contrasts per batch.
#' @param df_gram Full residualized Gram (default) or the historical diagonal
#'   shortcut.
#' @returns df (NA if undefined), G, and expected sampling variance.
#' @keywords internal
pffr_influence_df <- function(
  core,
  Xp,
  chunk_size = 32L,
  df_gram = c("full", "diagonal")
) {
  df_gram <- match.arg(df_gram)
  Xp <- as.matrix(Xp)
  if (!is.numeric(Xp) || ncol(Xp) != ncol(core$B) || any(!is.finite(Xp)))
    stop(
      "Xp must contain finite full-coefficient-space contrasts.",
      call. = FALSE
    )
  if (
    length(chunk_size) != 1L ||
      !is.finite(chunk_size) ||
      chunk_size < 1 ||
      chunk_size > .Machine$integer.max
  )
    stop(
      "chunk_size must be positive and representable as an integer.",
      call. = FALSE
    )
  n <- nrow(Xp)
  out <- rep(NA_real_, n)
  expected <- numeric(n)
  if (!n)
    return(list(df = out, G = core$G, expected_sampling_variance = expected))
  for (start in seq.int(1L, n, by = as.integer(chunk_size))) {
    jj <- seq.int(start, min(n, start + as.integer(chunk_size) - 1L))
    M <- core$B %*% t(Xp[jj, , drop = FALSE])
    q2 <- matrix(0, core$G, length(jj))
    ts <- lapply(seq_len(core$G), function(g) {
      block <- core$blocks[[g]]
      q <- block$A %*% (block$T %*% M)
      q2[g, ] <<- colSums(q^2)
      if (df_gram == "full") crossprod(block$T, q) else NULL
    })
    for (j in seq_along(jj)) {
      if (df_gram == "full") {
        T <- do.call(cbind, lapply(ts, function(x) x[, j]))
        Gamma <- diag(q2[, j], nrow = core$G) - crossprod(T, core$C %*% T)
        Gamma <- (Gamma + t(Gamma)) / 2
      } else {
        # Historical diagonal shortcut: Gamma = diag(||q_g||^2), so
        # tr^2 / tr(Gamma^2) = (sum_g ||q_g||^2)^2 / sum_g ||q_g||^4.
        Gamma <- diag(q2[, j], nrow = core$G)
      }
      tr <- sum(diag(Gamma))
      tr2 <- sum(Gamma^2)
      expected[jj[j]] <- core$correction * tr
      if (
        is.finite(tr) &&
          tr > 100 * .Machine$double.eps * sum(q2[, j]) &&
          is.finite(tr2) &&
          tr2 > 0
      )
        out[jj[j]] <- min(core$G, max(1, tr^2 / tr2))
    }
  }
  list(
    df = out,
    G = core$G,
    expected_sampling_variance = expected,
    reference = "central Gaussian working-model sampling quadratic form"
  )
}

#' Cached fixed-fit influence object for a pffr model
#' @param object Fitted pffr model.
#' @param sandwich NULL inherits the fit; otherwise cluster or cl2.
#' @param cluster Optional per-curve grouping override.
#' @param cl2_adjustment NULL inherits the fit; otherwise auto, exact or shortcut.
#' @param leverage_cap,tol Numerical adjustment settings.
#' @param dof_correction,edf_type CR1 small-sample correction for
#'   `sandwich = "cluster"`; `NULL` inherits the fit. Both enter the cache key
#'   for `sandwich = "cluster"`, because there they scale the cached
#'   `correction`; for `sandwich = "cl2"` they are never applied and are keyed
#'   as `"none"` so the same object is reused whatever they are set to.
#' @returns A `pffr_influence` object: the symmetrized penalized bread `B`, the
#'   residualization matrix `C`, the per-cluster residual influence columns `K`,
#'   the compressed per-cluster geometry `blocks` (`T` and the leverage weight
#'   `A`), the cluster count `G` and labels `groups`, the resolved `adjustment`,
#'   the finite-sample `correction` (`G/(G-1)` times any CR1 dof factor), the
#'   Bayesian smoothing-bias term `B2`, the per-cluster `diagnostics` (see
#'   [pffr_influence_core()]) and the numerical settings. Cached on the fit
#'   unless a `cluster` override is supplied. Internal research prototype:
#'   smoothing-selection uncertainty is absent and the object is conditional on
#'   the fitted smoothing parameters and weights.
#' @keywords internal
pffr_influence <- function(
  object,
  sandwich = NULL,
  cluster = NULL,
  cl2_adjustment = NULL,
  leverage_cap = .999,
  tol = 1e-8,
  dof_correction = NULL,
  edf_type = NULL
) {
  type <- sandwich %||% pffr_canonicalize_cov(object)$fit_type
  if (!type %in% c("cluster", "cl2"))
    stop(
      "Fixed-fit cluster influence requires sandwich='cluster' or 'cl2'.",
      call. = FALSE
    )
  b <- pffr_model_based_gam(object)
  kind <- pffr_score_kind(b$family)
  if (kind == "custom")
    stop("No cluster-robust score path for this family.", call. = FALSE)
  if (kind == "approx") pffr_warn_approx_score(b$family)
  cid <- build_cluster_id(object$pffr, cluster = cluster)
  adjustment <- if (type == "cluster") "none" else
    resolve_cl2_adjustment(
      cl2_adjustment %||% object$pffr$cl2_adjustment %||% "auto",
      G = length(unique(cid)),
      maxDg = max(table(cid)) * if (kind == "gaulss") 2 else 1,
      p = ncol(b$Vp)
    )
  dof_correction <- dof_correction %||% object$pffr$dof_correction %||% "none"
  edf_type <- edf_type %||% object$pffr$edf_type %||% "trace"
  # dof_correction / edf_type scale $correction (and hence every
  # expected_sampling_variance read off this object), so they must be part of
  # the cache key: otherwise an explicit override would silently reuse the
  # entry built for the fit's own setting. They only ever apply on the CR1
  # path (see the type == "cluster" guard below), so cl2 keys them as "none":
  # otherwise an edf-fitted object queried as cl2 would cache a second,
  # byte-identical entry for every (dof_correction, edf_type) combination.
  key_dof <- if (type == "cluster") dof_correction else "none"
  key_edf <- if (type == "cluster") edf_type else "none"
  key <- paste(
    "influence",
    type,
    adjustment,
    leverage_cap,
    tol,
    key_dof,
    key_edf,
    sep = "|"
  )
  cache <- object$pffr$Vsandwich_cache
  if (is.null(cluster) && is.environment(cache) && !is.null(cache[[key]]))
    return(cache[[key]])
  work <- switch(
    kind,
    gaulss = build_cl2_working_gaulss(b, cid),
    scat = build_cl2_working_scat(b, cid),
    build_cl2_working_standard(b, cid)
  )
  core <- pffr_influence_core(
    work$Xw,
    b$Vp,
    work$cluster_id,
    work$z,
    adjustment,
    leverage_cap,
    tol
  )
  core$B2 <- (b$Vp + t(b$Vp)) / 2 - b$Ve
  core$score_kind <- kind
  core$conditional_on_smoothing <- TRUE
  if (type == "cluster")
    core$correction <- core$correction *
      compute_dof_factor(b, cid, dof_correction, edf_type)
  if (is.null(cluster) && is.environment(cache)) cache[[key]] <- core
  core
}
