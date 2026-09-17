# The pre-2026-09 Satterthwaite df kernel, kept here VERBATIM as a frozen
# reference (body copied from `git show exact-cl2-default:R/pffr-core.R`, only
# the function name changed). It always used the shortcut leverage weight
# A_g = (I - H_gg)^{-1/2}, whatever covariance was requested, so it pins down
# exactly which of today's `df_gram = "diagonal"` results are backwards
# compatible: those computed with `cl2_adjustment = "shortcut"`.
historical_satterthwaite_df <- function(
  Xw,
  cluster_id,
  Vp,
  Xp,
  use_cl2,
  leverage_cap = 0.999,
  tol = 1e-8
) {
  M <- Vp %*% t(Xp) # p x n_points
  n_pts <- ncol(M)
  s2 <- numeric(n_pts)
  s4 <- numeric(n_pts)
  groups <- unique(cluster_id)
  for (g in groups) {
    idx <- which(cluster_id == g)
    Xwg <- Xw[idx, , drop = FALSE]
    Qg <- Xwg %*% M # D_g x n_points
    if (use_cl2) {
      # Reproduce the shipped CL2 leverage adjustment exactly (same capping as
      # gam_sandwich_cluster_cl2()): A_g = (I - H_gg)^{-1/2}.
      Hgg <- Xwg %*% Vp %*% t(Xwg)
      Hgg <- 0.5 * (Hgg + t(Hgg))
      ee <- eigen(Hgg, symmetric = TRUE)
      if (any(ee$values > leverage_cap, na.rm = TRUE)) {
        ee$values <- pmin(ee$values, leverage_cap)
      }
      Mg <- diag(length(idx)) -
        ee$vectors %*%
          diag(ee$values, nrow = length(ee$values)) %*%
          t(ee$vectors)
      Qg <- sym_inv_sqrt(Mg, tol = tol) %*% Qg
    }
    cn2 <- colSums(Qg^2) # ||q_g||^2 per evaluation point
    s2 <- s2 + cn2
    s4 <- s4 + cn2^2
  }
  G <- length(groups)
  df <- s2^2 / s4
  df[!is.finite(df)] <- NA_real_
  # Bounds: 1 <= df <= G (up to rounding); leave NA (zero-variance) untouched.
  ok <- is.finite(df)
  df[ok] <- pmin(pmax(df[ok], 1), G)
  list(df = df, G = G)
}
