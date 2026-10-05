# Fixed-fit cluster inference core, GPL (>= 2), as in refund.
# Does not propagate smoothing selection or remove smoothing bias.

#' Compressed fixed-fit cluster influence geometry
#'
#' Uses Z_g = Q_g T_g and R_g = I - Z_g C Z_g', C = 2B - B Z'Z B, i.e. the
#' full-block ("exact") CL2 leverage adjustment
#' \eqn{A_g = \{(I - H)^2\}_{gg}^{-1/2}}{A_g = ((I - H)^2)_gg^(-1/2)}.
#' Numerical rank is measured; no universal response-basis rank bound is used.
#' @param Xw Finite working-likelihood-scaled design.
#' @param Vp Penalized bread in the same scaling.
#' @param cluster_id One nonmissing cluster label per working row.
#' @param z Optional finite working residuals.
#' @param leverage_cap Numerical cap strictly between zero and one; the
#'   eigenvalues of the residual blocks are floored at `(1 - leverage_cap)^2`.
#' @param rank_tol Relative SVD tolerance; NULL uses dimension times machine epsilon.
#' @param df_precompute_bytes Memory budget for the per-cluster residualization
#'   blocks \eqn{R T_g^\top} used by [pffr_influence_df()], where \eqn{C = R^\top R}.
#'   They cost `8 * p * sum(rank_g)` bytes and let the df evaluate its Gram as a
#'   symmetric rank-k update without any \eqn{p \times p} product. Above the
#'   budget, or when `C` is not usably positive definite, the blocks are not
#'   stored and the df falls back to the general product.
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
  leverage_cap = .999,
  rank_tol = NULL,
  df_precompute_bytes = 2^28
) {
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
      leverage_cap >= 1
  )
    stop("leverage_cap must be in (0,1).", call. = FALSE)
  if (
    !is.null(rank_tol) &&
      (length(rank_tol) != 1L ||
        !is.finite(rank_tol) ||
        rank_tol < 0 ||
        rank_tol >= 1)
  )
    stop("rank_tol must be in [0,1).", call. = FALSE)
  if (
    length(df_precompute_bytes) != 1L ||
      !is.numeric(df_precompute_bytes) ||
      is.na(df_precompute_bytes) ||
      df_precompute_bytes < 0
  )
    stop(
      "df_precompute_bytes must be a single nonnegative number.",
      call. = FALSE
    )
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
    small <- diag(length(keep)) - T %*% C %*% t(T)
    ee <- eigen((small + t(small)) / 2, symmetric = TRUE)
    values <- if (length(keep) < length(idx)) c(ee$values, 1) else ee$values
    min_eig[g] <- min(values)
    # Upstream's definition |max v| / |min v| (NOT max|v| / min|v|): the two
    # differ exactly when a negative eigenvalue is present, i.e. in the
    # P-LB5 indefinite-bread case this monitor exists to expose.
    max_kappa[g] <- abs(max(values)) / max(abs(min(values)), 1e-300)
    min_eig_rel[g] <- min(values) / max(abs(max(values)), 1e-300)
    floor <- (1 - leverage_cap)^2
    n_floored[g] <- sum(ee$values < floor)
    A <- tcrossprod(
      sweep(ee$vectors, 2L, sqrt(pmax(ee$values, floor)), "/"),
      ee$vectors
    )
    blocks[[g]] <- list(T = T, A = A)
    if (!is.null(z))
      K[, g] <- B %*%
        crossprod(T, A %*% crossprod(dec$u[, keep, drop = FALSE], z[idx]))
  }
  # Residualization factor for the moment df. pffr_influence_df() needs the
  # Gram t_g' C t_h for every cluster pair and every contrast. For a genuine
  # penalized bread C = B + B S B is positive definite, so with C = R'R that
  # Gram is crossprod(R T): a symmetric rank-k update, half the flops of the
  # general product, and the per-contrast O(p^2 G) product C %*% T disappears
  # because R T_g' is cached here, per cluster. The Gram, not that product, is
  # what dominates the df (87% of it at G = 200), which is why caching C T_g'
  # alone -- the obvious precompute -- buys about 1.1x while this buys about
  # 2.1x; see inst/benchmarks/df-timing.R.
  # The factorization is verified rather than assumed: a contrived or
  # numerically inconsistent bread can make C indefinite, or so ill conditioned
  # that R'R no longer reproduces it, and the df then falls back to the general
  # product. The byte budget bounds the cached blocks (8 * p * sum_g r_g).
  within_budget <- 8 * p * sum(rank) <= df_precompute_bytes
  Rchol <- if (within_budget) tryCatch(chol(C), error = function(e) NULL) else
    NULL
  if (!is.null(Rchol) && max(abs(crossprod(Rchol) - C)) > 1e-10 * max(abs(C)))
    Rchol <- NULL
  df_precompute <- !is.null(Rchol)
  if (df_precompute)
    for (g in seq_len(G)) blocks[[g]]$RT <- Rchol %*% t(blocks[[g]]$T)
  structure(
    list(
      B = B,
      C = C,
      K = K,
      blocks = blocks,
      G = G,
      groups = groups,
      adjustment = "exact",
      df_precompute = df_precompute,
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
      rank_tol = rank_tol,
      version = "fixed-fit-core-2026-09-09"
    ),
    class = "pffr_influence"
  )
}

#' Assemble the CL2 covariance from fixed-fit influence columns
#'
#' Returns the Bayesian form \eqn{G/(G-1) K K^\top + (V_p - V_e)}.
#' @param core Influence object with residual columns K.
#' @returns Covariance with raw numerical diagnostics.
#' @keywords internal
pffr_influence_vcov <- function(core) {
  if (is.null(core$K))
    stop("Residual influence columns were not computed.", call. = FALSE)
  V <- core$correction * tcrossprod(core$K) + core$B2
  V <- (V + t(V)) / 2
  d <- core$diagnostics
  attr(V, "cl2_adjustment") <- "exact"
  attr(V, "n_adjusted") <- sum(d$n_floored > 0L)
  attr(V, "max_leverage") <- max(d$max_leverage)
  attr(V, "min_block_eig") <- min(d$min_block_eig)
  attr(V, "max_block_kappa") <- max(d$max_block_kappa)
  attr(V, "cluster_rank") <- d$rank
  attr(V, "inference_core_version") <- core$version
  # Study-LB P-LB5 hat-invariant monitors: h_ii in [0, 1], eigen(H_gg) >= 0 and
  # a positive semi-definite residual block.
  attr(V, "max_obs_leverage") <- max(d$max_obs_leverage)
  attr(V, "min_obs_leverage") <- min(d$min_obs_leverage)
  attr(V, "min_hat_eig") <- min(d$min_hat_eig)
  attr(V, "min_block_eig_rel") <- min(d$min_block_eig_rel)
  attr(V, "hat_invariant_violation") <- pffr_hat_invariant_violation(
    max_obs_leverage = max(d$max_obs_leverage),
    min_block_eig_rel = min(d$min_block_eig_rel),
    min_obs_leverage = min(d$min_obs_leverage),
    min_hat_eig = min(d$min_hat_eig)
  )
  V
}

#' Central Gaussian sampling-variance moment degrees of freedom
#'
#' Gamma_gh = 1(g=h)||q_g||^2 - t_g' C t_h retains fit residualization, and
#' \eqn{\nu(a) = \{\mathrm{tr}(\Gamma)\}^2/\mathrm{tr}(\Gamma^2)}.
#' This is conditional on weights and smoothing parameters, and is not an
#' exact t law. Noncentral means, B2 and smoothing selection are not covered.
#'
#' Per contrast the residualized Gram \eqn{T^\top C\,T} dominates the cost. When
#' the core carries the cached blocks \eqn{R T_g^\top} with \eqn{C = R^\top R}
#' (see [pffr_influence_core()]'s `df_precompute_bytes`), the Gram is
#' \eqn{V^\top V} for a \eqn{p \times G} matrix \eqn{V}, and
#' \eqn{\mathrm{tr}(\Gamma^2)} needs only \eqn{\lVert V^\top V\rVert_F^2 =
#' \lVert V V^\top\rVert_F^2}, computed from the smaller of the two products
#' (\eqn{G \times G} or \eqn{p \times p}), so the cost per contrast grows
#' with \eqn{\min(G, p)\,G\,p} rather than \eqn{G^2 p}. Cores without the
#' cached blocks use the general product and return the same numbers.
#' @param core Fixed-fit influence object.
#' @param Xp Finite full-coefficient contrasts, one per row.
#' @param chunk_size Positive number of contrasts per batch.
#' @returns df (NA if undefined), G, and expected sampling variance.
#' @keywords internal
pffr_influence_df <- function(core, Xp, chunk_size = 32L) {
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
  p <- ncol(core$B)
  G <- core$G
  # With the cached per-cluster blocks R T_g' (pffr_influence_core()'s
  # df_precompute, C = R'R) the residualized Gram is V'V. Cores without them --
  # a bread whose C is not usably positive definite, or a geometry above the
  # memory budget -- keep the general path.
  factored <- isTRUE(core$df_precompute)
  for (start in seq.int(1L, n, by = as.integer(chunk_size))) {
    jj <- seq.int(start, min(n, start + as.integer(chunk_size) - 1L))
    nj <- length(jj)
    M <- core$B %*% t(Xp[jj, , drop = FALSE])
    q2 <- matrix(0, G, nj)
    # Column j of vmat holds the p x G matrix [v_1 ... v_G] for contrast jj[j],
    # flattened column-major, with v_g = R t_g (factored) or v_g = t_g.
    vmat <- matrix(0, p * G, nj)
    for (g in seq_len(G)) {
      block <- core$blocks[[g]]
      q <- block$A %*% (block$T %*% M)
      q2[g, ] <- colSums(q^2)
      vmat[seq.int((g - 1L) * p + 1L, g * p), ] <- if (factored)
        block$RT %*% q else crossprod(block$T, q)
    }
    for (j in seq_len(nj)) {
      V <- matrix(vmat[, j], p, G)
      d <- q2[, j]
      if (factored) {
        # Gamma = diag(d) - V'V: tr = sum(d) - ||V||^2 and
        # tr(Gamma^2) = sum(d^2) - 2 sum_g d_g ||v_g||^2 + ||V'V||_F^2.
        vnorm2 <- colSums(V^2)
        frob <- if (G <= p) sum(crossprod(V)^2) else sum(tcrossprod(V)^2)
        tr <- sum(d) - sum(vnorm2)
        tr2 <- sum(d^2) - 2 * sum(d * vnorm2) + frob
      } else {
        Gamma <- diag(d, nrow = G) - crossprod(V, core$C %*% V)
        # The general product is symmetric only up to rounding.
        Gamma <- (Gamma + t(Gamma)) / 2
        tr <- sum(diag(Gamma))
        tr2 <- sum(Gamma^2)
      }
      expected[jj[j]] <- core$correction * tr
      if (
        is.finite(tr) &&
          tr > 100 * .Machine$double.eps * sum(d) &&
          is.finite(tr2) &&
          tr2 > 0
      )
        out[jj[j]] <- min(G, max(1, tr^2 / tr2))
    }
  }
  list(
    df = out,
    G = G,
    expected_sampling_variance = expected,
    reference = "central Gaussian working-model sampling quadratic form"
  )
}

#' Cached fixed-fit CL2 influence object for a pffr model
#' @param object Fitted pffr model.
#' @param cluster Optional per-curve grouping override.
#' @param leverage_cap Numerical floor setting, see [pffr_influence_core()].
#' @returns A `pffr_influence` object: the symmetrized penalized bread `B`, the
#'   residualization matrix `C`, the per-cluster residual influence columns `K`,
#'   the compressed per-cluster geometry `blocks` (`T`, the leverage weight `A`
#'   and, when the df precompute applies, the residualization block `RT` =
#'   \eqn{R T_g^\top}), the cluster count `G` and labels `groups`,
#'   the finite-sample `correction` `G/(G-1)`, the Bayesian smoothing-bias term
#'   `B2` \eqn{= V_p - V_e}, the per-cluster `diagnostics` (see
#'   [pffr_influence_core()]) and the numerical settings. Cached on the fit
#'   unless a `cluster` override is supplied. The object is conditional on the
#'   fitted smoothing parameters and weights.
#' @keywords internal
pffr_influence <- function(object, cluster = NULL, leverage_cap = .999) {
  pffr_check_sandwich_ar1(object, TRUE)
  b <- pffr_model_based_gam(object)
  kind <- pffr_score_kind(b$family)
  if (kind == "custom")
    stop("No cluster-robust score path for this family.", call. = FALSE)
  if (kind == "approx") pffr_warn_approx_score(b$family)
  key <- paste("influence", leverage_cap, sep = "|")
  cache <- object$pffr$Vsandwich_cache
  if (is.null(cluster) && is.environment(cache) && !is.null(cache[[key]]))
    return(cache[[key]])
  cid <- build_cluster_id(object$pffr, cluster = cluster)
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
    leverage_cap = leverage_cap
  )
  core$B2 <- (b$Vp + t(b$Vp)) / 2 - b$Ve
  core$score_kind <- kind
  core$conditional_on_smoothing <- TRUE
  if (is.null(cluster) && is.environment(cache)) cache[[key]] <- core
  core
}
