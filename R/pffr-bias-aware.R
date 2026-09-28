# Bias-aware pointwise intervals for pffr fits.
#
# The interval for a linear estimand L theta of a fit (typically selected by
# curve-blocked NCV) is
#   L theta -/+ z * sqrt(se^2 + delta^2),  delta = L (theta - theta_ref),
# where theta_ref comes from a reference fit of the same model (typically REML)
# and se from the chosen covariance (recommended: exact CL2, Bayesian form).

# Helpers ------------------------------------------------------------------------

#' Check that a bias reference fit matches a pffr fit
#'
#' The bias-aware interval compares two fits of the *same* model and data that
#' differ only in how the smoothing parameters were chosen, so that their
#' coefficient vectors live in the same basis. This checks the coefficient
#' layout, the model frame (response, covariates and functional-covariate
#' matrices, in the same row order), prior weights and offsets, the smooth bases
#' (labels, classes, basis dimensions, knots and penalty matrices, which carry the
#' identifiability constraints), `ffpc`/`pcre` metadata, and family and link.
#' Family parameters estimated during fitting (e.g. the negative binomial
#' \eqn{\theta}{theta}) may differ between the fits.
#'
#' @param object The fit whose intervals are computed.
#' @param bias_ref The reference fit.
#' @returns `object$coefficients - bias_ref$coefficients`.
#' @keywords internal
pffr_bias_ref_difference <- function(object, bias_ref) {
  if (!inherits(bias_ref, "pffr")) {
    stop("`bias_ref` must be a fitted pffr model.", call. = FALSE)
  }
  mismatch <- function(what) {
    stop(
      "`bias_ref` is not a fit of the same model and data as `object` (",
      what,
      " differ). Refit the reference with the same formula, bases and data, ",
      "changing only `method`.",
      call. = FALSE
    )
  }
  if (!identical(names(object$coefficients), names(bias_ref$coefficients))) {
    mismatch("coefficient names")
  }
  if (
    !identical(
      pffr_family_base(object$family),
      pffr_family_base(bias_ref$family)
    )
  ) {
    mismatch("families or links")
  }
  same <- function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE))
  meta <- c(
    "nobs",
    "yind",
    "is_sparse",
    "missing_indices",
    "ffpc",
    "pcre_terms"
  )
  if (!same(object$pffr[meta], bias_ref$pffr[meta])) {
    mismatch("response grids, missing-value patterns or ffpc/pcre bases")
  }
  if (!same(object$model, bias_ref$model)) {
    mismatch("model frames (responses or covariates)")
  }
  if (
    !same(object$prior.weights, bias_ref$prior.weights) ||
      !same(object$offset, bias_ref$offset)
  ) {
    mismatch("weights or offsets")
  }
  if (!pffr_same_smooths(object$smooth, bias_ref$smooth)) {
    mismatch("smooth terms")
  }
  object$coefficients - bias_ref$coefficients
}

# Family name without fitted parameters (mgcv writes e.g. "Negative
# Binomial(5.2)" into a fitted nb family), plus the link.
pffr_family_base <- function(family) {
  c(tolower(sub("\\(.*$", "", family$family)), family$link)
}

# Smooth terms match if labels, classes, coefficient ranges, basis dimensions,
# knots and penalty matrices agree.
pffr_same_smooths <- function(smooths, ref_smooths) {
  if (length(smooths) != length(ref_smooths)) return(FALSE)
  same <- function(a, b) {
    identical(a$label, b$label) &&
      identical(class(a), class(b)) &&
      identical(c(a$first.para, a$last.para), c(b$first.para, b$last.para)) &&
      identical(a$bs.dim, b$bs.dim) &&
      isTRUE(all.equal(a$S, b$S)) &&
      isTRUE(all.equal(pffr_smooth_knots(a), pffr_smooth_knots(b)))
  }
  all(mapply(same, smooths, ref_smooths))
}

pffr_smooth_knots <- function(sm) {
  if (!is.null(sm$margin)) return(lapply(sm$margin, pffr_smooth_knots))
  sm$knots
}

#' Combine a variance-only standard error with a bias estimate
#'
#' @param se Standard errors (variance part).
#' @param delta Bias estimates of the same length.
#' @returns `sqrt(se^2 + delta^2)`.
#' @keywords internal
pffr_bias_aware_se <- function(se, delta) {
  sqrt(se^2 + delta^2)
}

# Pointwise intervals for the linear predictor / conditional mean -----------------

#' Pointwise confidence intervals for pffr predictions, optionally bias-aware
#'
#' Pointwise Wald intervals for the linear predictor, or for the conditional
#' mean \eqn{E(Y(t) \mid X)}{E(Y(t)|X)}, of a [pffr()] fit, with a choice of
#' covariance (model-based or sandwich, as in [coef.pffr()]) and an optional
#' bias-aware widening computed from a second fit of the same model.
#'
#' Intervals are always built on the link scale,
#' \deqn{\hat\eta \pm z_{(1+\mathrm{level})/2}\sqrt{\mathrm{se}^2 +
#'   \delta^2},}{eta_hat -/+ z * sqrt(se^2 + delta^2),}
#' with \eqn{\delta = 0} unless `bias_ref` is supplied. For
#' `type = "response"` the estimate and both endpoints are mapped through the
#' inverse link (so the interval is not symmetric around the estimate); `se`
#' and `delta` stay on the link scale.
#'
#' @section Bias-aware intervals:
#' Smoothing-parameter selection by curve-blocked neighbourhood
#' cross-validation (`pffr(method = "NCV")`) is robust to within-curve
#' dependence but smooths more than REML, so its smoothing bias is no longer
#' negligible and a variance-only interval undercovers. The bias-aware interval
#' adds the difference between the two estimates in quadrature:
#' \deqn{\delta = L(\hat\theta_{\mathrm{NCV}} - \hat\theta_{\mathrm{REML}}),}{
#'   delta = L (theta_NCV - theta_REML),}
#' for the rows \eqn{L} of the prediction matrix. The recommended recipe is
#' an NCV fit with curve blocks (the default `ncv_blocks = "cluster"`) as
#' `object`, the exact CL2 sandwich in its Bayesian form
#' (`sandwich = "cl2", cl2_adjustment = "exact", freq = FALSE`) and a REML fit
#' of the same model as `bias_ref`. Fit both with `sandwich = "none"` to avoid
#' computing a fit-time sandwich that is not used.
#'
#' \eqn{\delta} estimates only the part of the smoothing bias in which the two
#' fits differ. Bias that both fits share -- a basis too small for the truth,
#' or both fits oversmoothing a rough truth -- is not covered.
#'
#' @param object A fitted [pffr()] model with a single linear predictor.
#' @param newdata Optional prediction data in the format supplied to [pffr()]
#'   (as in [predict.pffr()]). `NULL` (default) evaluates at the fitted
#'   observation points. Fits with a model offset are supported at the fitted
#'   points only.
#' @param type `"link"` (default) for the linear predictor, `"response"` for the
#'   conditional mean.
#' @param level Confidence level, defaults to `0.95`.
#' @param sandwich,freq,cluster,dof_correction,edf_type,cl2_adjustment
#'   Covariance choice, as in [coef.pffr()]. `sandwich = NULL` inherits the
#'   fit-time choice.
#' @param bias_ref Optional reference fit of the same model and data (typically
#'   the REML fit when `object` is the NCV fit), differing only in how the
#'   smoothing parameters were chosen. If supplied, intervals are bias-aware
#'   (see the section above). Any second fit of the same model is accepted;
#'   \eqn{\delta}{delta} is then the contrast between the two estimators, which
#'   is a smoothing-bias proxy only for the NCV-versus-REML pairing.
#' @returns A data frame with one row per evaluation point, in the row order of
#'   `predict(object, type = "lpmatrix")` (curve-major, index fastest; rows with
#'   missing responses are omitted and sparse responses keep the order of
#'   `ydata` when `newdata = NULL`), with columns `.obs` (curve), `.index`
#'   (value of the response index), `fit` (on the scale of `type`), `se_link`
#'   (variance part, link scale), `delta_link` (link scale; only with
#'   `bias_ref`), `lower` and `upper` (on the scale of `type`). Attributes
#'   `type`, `level`, `crit_value` and `bias_ref_method` (the reference fit's
#'   smoothing-parameter method, or `NA`).
#' @seealso [coef.pffr()] (argument `bias_ref`) for bias-aware intervals of
#'   coefficient functions, [pffr_jackknife_se()].
#' @export
#' @author Fabian Scheipl
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- pffr_simulate(Y ~ ff(X1), n = 30, nxgrid = 15, nygrid = 25)
#' yind <- attr(d, "yindex")
#' fit_ncv <- pffr(Y ~ ff(X1), yind = yind, data = d, method = "NCV",
#'                 sandwich = "none")
#' fit_reml <- pffr(Y ~ ff(X1), yind = yind, data = d, method = "REML",
#'                  sandwich = "none")
#' ci <- pffr_predict_ci(fit_ncv, sandwich = "cl2", cl2_adjustment = "exact",
#'                       bias_ref = fit_reml)
#' head(ci)
#' }
pffr_predict_ci <- function(
  object,
  newdata = NULL,
  type = c("link", "response"),
  level = 0.95,
  sandwich = NULL,
  freq = FALSE,
  cluster = NULL,
  dof_correction = NULL,
  edf_type = NULL,
  cl2_adjustment = NULL,
  bias_ref = NULL
) {
  if (!inherits(object, "pffr")) {
    stop("`object` must be a fitted pffr model.", call. = FALSE)
  }
  type <- match.arg(type)
  if (
    !is.numeric(level) ||
      length(level) != 1 ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1
  ) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (isTRUE((object$family$nlp %||% 1L) > 1L)) {
    stop(
      "pffr_predict_ci() supports single linear-predictor families only.",
      call. = FALSE
    )
  }
  theta_diff <- if (!is.null(bias_ref)) {
    pffr_bias_ref_difference(object, bias_ref)
  }

  lp <- pffr_prediction_design(object, newdata)
  V <- pffr_vcov(
    object,
    sandwich = sandwich,
    freq = freq,
    cluster = cluster,
    dof_correction = dof_correction,
    edf_type = edf_type,
    cl2_adjustment = cl2_adjustment
  )
  se <- sqrt(as.vector(rowSums((lp$X %*% V) * lp$X)))
  delta <- if (!is.null(theta_diff)) as.vector(lp$X %*% theta_diff)
  crit_value <- stats::qnorm((1 + level) / 2)
  half <- crit_value * if (is.null(delta)) se else pffr_bias_aware_se(se, delta)
  lower <- lp$eta - half
  upper <- lp$eta + half
  fit <- lp$eta
  if (type == "response") {
    linkinv <- object$family$linkinv
    ends <- cbind(linkinv(lower), linkinv(upper))
    # A decreasing inverse link (e.g. the inverse link) swaps the endpoints.
    lower <- pmin(ends[, 1], ends[, 2])
    upper <- pmax(ends[, 1], ends[, 2])
    fit <- linkinv(fit)
  }
  out <- lp$points
  out$fit <- fit
  out$se_link <- se
  if (!is.null(delta)) out$delta_link <- delta
  out$lower <- lower
  out$upper <- upper
  attr(out, "type") <- type
  attr(out, "level") <- level
  attr(out, "crit_value") <- crit_value
  attr(out, "bias_ref_method") <- if (is.null(bias_ref)) {
    NA_character_
  } else {
    bias_ref$method %||% NA_character_
  }
  out
}

# Prediction matrix, linear predictor and evaluation points (.obs, .index) at
# the fitted points or at newdata. At the fitted points the stored linear
# predictor includes any offset exactly.
pffr_prediction_design <- function(object, newdata) {
  meta <- object$pffr
  if (is.null(newdata)) {
    X <- predict(object, type = "lpmatrix", reformat = FALSE)
    eta <- as.vector(object$linear.predictors)
    points <- if (isTRUE(meta$is_sparse)) {
      # The model frame omits ydata rows with a missing .value.
      meta$ydata[!is.na(meta$ydata$.value), c(".obs", ".index")]
    } else {
      pffr_grid_points(meta$nobs, meta$yind, meta$missing_indices)
    }
  } else {
    if (!is.null(object$offset) && any(object$offset != 0)) {
      stop(
        "This fit uses a model offset: intervals for `newdata` are not ",
        "supported. Use `newdata = NULL` for the fitted points.",
        call. = FALSE
      )
    }
    X <- predict(object, newdata = newdata, type = "lpmatrix", reformat = FALSE)
    eta <- as.vector(X %*% object$coefficients)
    points <- pffr_grid_points(nrow(X) / length(meta$yind), meta$yind)
  }
  if (!is.null(attr(X, "lpi"))) {
    stop(
      "pffr_predict_ci() supports single linear-predictor families only.",
      call. = FALSE
    )
  }
  if (nrow(points) != nrow(X) || length(eta) != nrow(X)) {
    stop(
      "Could not align evaluation points with the prediction matrix.",
      call. = FALSE
    )
  }
  rownames(points) <- NULL
  list(X = X, eta = eta, points = points)
}

# Evaluation points of a dense response grid in lpmatrix order (curve-major),
# without the rows of missing responses.
pffr_grid_points <- function(nobs, yind, missing_indices = NULL) {
  points <- data.frame(
    .obs = rep(seq_len(nobs), each = length(yind)),
    .index = rep(yind, times = nobs)
  )
  if (length(missing_indices)) points <- points[-missing_indices, ]
  points
}
