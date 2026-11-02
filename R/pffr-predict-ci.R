# Pointwise confidence intervals for pffr predictions.

#' Pointwise confidence intervals for pffr predictions
#'
#' Pointwise Wald intervals for the linear predictor, or for the conditional
#' mean \eqn{E(Y(t) \mid X)}{E(Y(t)|X)}, of a [pffr()] fit, with the same
#' covariance and critical values as [coef.pffr()].
#'
#' Intervals are built on the link scale, \eqn{\hat\eta \pm c\,\mathrm{se}}{
#' eta_hat -/+ c * se}. By default (fits with `sandwich = TRUE`) `se` comes
#' from the CL2 sandwich and \eqn{c} is a Satterthwaite critical value with a
#' separate df for every point (see the section \sQuote{Inference} of
#' [pffr()]); model-based standard errors use Gaussian critical values. For
#' `type = "response"` the estimate and both endpoints are mapped through the
#' inverse link (so the interval is not symmetric around the estimate); `se`
#' stays on the link scale.
#'
#' The full functional intercept \eqn{\alpha(t)}{alpha(t)} (level included)
#' is the linear predictor at covariate values at which all other terms vanish,
#' e.g. `X = 0` for `ff(X)` and `z = 0` for a linear effect of `z`; see the
#' examples and [coef.pffr()].
#'
#' @param object A fitted [pffr()] model with a single linear predictor.
#' @param newdata Optional prediction data in the format supplied to [pffr()]
#'   (as in [predict.pffr()]). `NULL` (default) evaluates at the fitted
#'   observation points. Fits with a model offset are supported at the fitted
#'   points only.
#' @param type `"link"` (default) for the linear predictor, `"response"` for the
#'   conditional mean.
#' @param level Confidence level, defaults to `0.95`.
#' @param sandwich,cluster Covariance choice, as in [coef.pffr()]:
#'   `sandwich = NULL` uses the fit's covariance, `TRUE` the CL2 sandwich and
#'   `FALSE` the model-based covariance.
#' @returns A data frame with one row per evaluation point, in the row order of
#'   `predict(object, type = "lpmatrix")` (curve-major, index fastest; rows with
#'   missing responses are omitted and sparse responses keep the order of
#'   `ydata` when `newdata = NULL`), with columns `.obs` (curve), `.index`
#'   (value of the response index), `fit` (on the scale of `type`), `se_link`
#'   (link scale), `df` (reference df of the critical value: `Inf` for Gaussian
#'   critical values), `lower` and `upper` (on the scale of `type`).
#'   Attributes `type`, `level`, `crit_value` (one value, or one per row for
#'   Satterthwaite critical values) and `crit_used`.
#' @seealso [coef.pffr()], [predict.pffr()].
#' @export
#' @author Fabian Scheipl
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- pffr_simulate(Y ~ ff(X1), n = 30, nxgrid = 15, nygrid = 25)
#' yind <- attr(d, "yindex")
#' fit <- pffr(Y ~ ff(X1), yind = yind, data = d)
#' ci <- pffr_predict_ci(fit)
#' head(ci)
#'
#' # Full functional intercept alpha(t) with its interval: the linear
#' # predictor of a curve with X1 = 0, where the ff() term vanishes.
#' zero <- data.frame(X1 = I(matrix(0, 1, ncol(d$X1))))
#' alpha <- pffr_predict_ci(fit, newdata = zero)
#' head(alpha[, c(".index", "fit", "lower", "upper")])
#' }
pffr_predict_ci <- function(
  object,
  newdata = NULL,
  type = c("link", "response"),
  level = 0.95,
  sandwich = NULL,
  cluster = NULL
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
  use_sandwich <- pffr_sandwich_arg(sandwich) %||%
    !identical(pffr_canonicalize_cov(object)$fit_type, "none")
  pffr_inform_ncv_intervals(object)

  lp <- pffr_prediction_design(object, newdata)
  V <- pffr_vcov(object, sandwich = use_sandwich, cluster = cluster)
  se <- sqrt(as.vector(rowSums((lp$X %*% V) * lp$X)))
  cv <- pffr_pointwise_crit(
    object,
    X = lp$X,
    level = level,
    use_sandwich = use_sandwich,
    cluster = cluster
  )
  half <- cv$crit * se
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
  out$df <- cv$df
  out$lower <- lower
  out$upper <- upper
  attr(out, "type") <- type
  attr(out, "level") <- level
  attr(out, "crit_value") <- cv$crit
  attr(out, "crit_used") <- cv$mode
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
