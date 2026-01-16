# =============================================================================
# Bootstrap Resampling Helpers (pure functions)
# =============================================================================

# Resample grid data by observation indices
# @param modcall Model call object to modify
# @param data Original data frame
# @param indices Bootstrap sample indices
# @returns Modified modcall with resampled data
resample_grid_observations <- function(modcall, data, indices) {
  modcall$data <- data[indices, , drop = FALSE]
  modcall
}

# Resample sparse/irregular data by observation indices
# Handles bootstrap duplicates by replicating ydata rows appropriately
# @param modcall Model call object to modify
# @param data Original covariate data frame
# @param ydata Original sparse response data frame
# @param indices Bootstrap sample indices
# @returns Modified modcall with resampled data and ydata
resample_sparse_observations <- function(modcall, data, ydata, indices) {
  data_boot <- data[indices, , drop = FALSE]
  # Set rownames to consecutive integers (pffr validates ydata$.obs against rownames)
  rownames(data_boot) <- seq_len(nrow(data_boot))

  # For each bootstrap index, get all ydata rows for that observation
  # This correctly handles when the same observation is sampled multiple times
  ydata_list <- lapply(seq_along(indices), \(i) {
    obs_i <- indices[i]
    ydata_rows <- ydata[ydata$.obs == obs_i, , drop = FALSE]
    ydata_rows$.obs <- i
    ydata_rows
  })
  ydata_boot <- do.call(rbind, ydata_list)
  rownames(ydata_boot) <- NULL

  modcall$data <- data_boot
  modcall$ydata <- ydata_boot
  modcall
}

# Resample using residual bootstrap (fitted + resampled residuals)
# @param modcall Model call object to modify
# @param data Original data frame
# @param indices Bootstrap sample indices
# @param fitted_values Matrix of fitted values from original model
# @param residuals Matrix of residuals from original model
# @param response_name Name of response variable in data
# @param center Whether to center residuals at each time point
# @returns Modified modcall with synthetic response
resample_residuals <- function(
  modcall,
  data,
  indices,
  fitted_values,
  residuals,
  response_name,
  center = FALSE
) {
  resid_boot <- residuals[indices, , drop = FALSE]
  if (center) {
    resid_boot <- resid_boot -
      matrix(
        colMeans(resid_boot),
        nrow = nrow(resid_boot),
        ncol = ncol(resid_boot),
        byrow = TRUE
      )
  }
  data_boot <- data
  data_boot[[response_name]] <- fitted_values + resid_boot
  modcall$data <- data_boot
  modcall
}

# =============================================================================
# Bootstrap CI Computation Helpers
# =============================================================================

# Extract coefficients from a fitted pffr model as a vector
# @param model Fitted pffr model (or try-error)
# @param n1, n2, n3 Grid sizes for coefficient evaluation
# @param fallback Vector of NAs to return on failure
# @returns Numeric vector of coefficient values
extract_boot_coefficients <- function(model, n1, n2, n3, fallback) {
  if (inherits(model, "try-error")) {
    return(fallback)
  }

  coefs <- coef(model, se = FALSE, n1 = n1, n2 = n2, n3 = n3)
  coef_vec <- c(
    coefs$pterms[, "value"],
    lapply(coefs$smterms, \(x) x$value)
  )
  unlist(coef_vec)
}

# Compute bootstrap confidence intervals for a single coefficient
# @param boot_result Object returned by boot::boot
# @param index Which coefficient index
# @param conf Confidence levels (numeric vector)
# @param type Type of CI ("percent", "bca", etc.)
# @returns Matrix of CI bounds
compute_single_boot_ci <- function(boot_result, index, conf, type) {
  # boot.ci uses "perc" but returns "percent"
  type_arg <- if (type == "percent") "perc" else type
  ci <- boot::boot.ci(boot_result, conf = conf, index = index, type = type_arg)
  bounds <- ci[[type]][, 4:5, drop = FALSE]
  if (is.null(bounds)) {
    bounds <- matrix(NA, nrow = length(conf), ncol = 2)
  }
  as.vector(bounds)
}

# Compute bootstrap CIs for multiple coefficients
# @param boot_result Object returned by boot::boot
# @param indices Which coefficient indices to compute CIs for
# @param conf Confidence levels
# @param type Type of CI
# @param apply_fn Function to use for iteration (lapply, mclapply, or parLapply wrapper)
# @returns Matrix with CI bounds as rows
compute_boot_cis <- function(
  boot_result,
  indices,
  conf,
  type,
  apply_fn = lapply
) {
  ci_list <- apply_fn(indices, \(i) {
    compute_single_boot_ci(boot_result, i, conf, type)
  })
  ci_matrix <- do.call(cbind, ci_list)
  rownames(ci_matrix) <- paste0(
    c((1 - conf) / 2, 1 - (1 - conf) / 2) * 100,
    "%"
  )
  t(ci_matrix)
}

# =============================================================================
# Main Exported Functions
# =============================================================================

#' Simple bootstrap CIs for pffr
#'
#' This function resamples observations in the data set to obtain approximate CIs for different
#' terms and coefficient functions that correct for the effects of dependency and heteroskedasticity
#' of the residuals along the index of the functional response, i.e., it aims for correct inference
#' if the residuals along the index of the functional response are not i.i.d.
#'
#' @param object a fitted \code{\link{pffr}}-model
#' @param n1 see \code{\link{coef.pffr}}
#' @param n2 see \code{\link{coef.pffr}}
#' @param n3 see \code{\link{coef.pffr}}
#' @param B  number of bootstrap replicates, defaults to (a measly) 100
#' @param parallel see \code{\link[boot]{boot}}
#' @param cl see \code{\link[boot]{boot}}
#' @param ncpus see \code{\link[boot]{boot}}. Defaults to \code{getOption("boot.ncpus", 1L)} (like \code{boot}).
#' @param conf desired levels of bootstrap CIs, defaults to 0.90 and 0.95
#' @param type type of bootstrap interval, see \code{\link[boot]{boot.ci}}. Defaults to "percent" for percentile-based CIs.
#' @param method either "resample" (default) to resample response trajectories, or "residual" to resample responses as fitted values
#' plus residual trajectories or "residual.c" to resample responses as fitted values
#' plus residual trajectories that are centered at zero for each gridpoint.
#' @param showProgress TRUE/FALSE
#' @param ... further arguments handed to \code{\link[boot]{boot}} like e.g. \code{strata} for models with random effects.
#' @return a list with similar structure as the return value of \code{\link{coef.pffr}}, containing the
#'         original point estimates of the various terms along with their bootstrap CIs.
#' @author Fabian Scheipl
#' @importFrom boot boot boot.ci
#' @importFrom parallel makePSOCKcluster clusterSetRNGStream parLapply mclapply stopCluster
#' @export
pffr_coefboot <- function(
  object,
  n1 = 100,
  n2 = 40,
  n3 = 20,
  B = 100,
  ncpus = getOption("boot.ncpus", 1),
  parallel = c("no", "multicore", "snow"),
  cl = NULL,
  conf = c(.9, .95),
  type = "percent",
  method = c("resample", "residual", "residual.c"),
  showProgress = TRUE,
  ...
) {
  # Validate inputs

  method <- match.arg(method)
  parallel <- match.arg(parallel)
  stopifnot(B > 0, conf > 0, conf < 1)
  ncpus <- ncpus %||% getOption("boot.ncpus", 1L)

  # Extract model call components, evaluating any unevaluated language objects
  modcall <- prepare_modcall_for_bootstrap(object)

  # Get original coefficients (template for return structure)
  coef_template <- coef(object, se = FALSE, n1 = n1, n2 = n2, n3 = n3)

  # Compute fallback return vector (all NAs) for failed bootstrap replicates
  n_coefs <- length(unlist(c(
    coef_template$pterms[, "value"],
    lapply(coef_template$smterms, \(x) x$value)
  )))
  fallback_vec <- rep(NA_real_, n_coefs)

  # Build resampling function based on data type and method
  resample_fn <- build_resample_function(
    object = object,
    modcall = modcall,
    method = method
  )

  # Bootstrap statistic function
  boot_statistic <- function(
    data,
    indices,
    modcall,
    n1,
    n2,
    n3,
    fallback,
    resample_fn,
    show_progress
  ) {
    modcall_boot <- resample_fn(modcall, data, indices)
    model_boot <- try(eval(as.call(modcall_boot)), silent = TRUE)
    coefs <- extract_boot_coefficients(model_boot, n1, n2, n3, fallback)
    if (show_progress) cat(".")
    coefs
  }

  # Run bootstrap
  message("starting bootstrap (", B, " replications).\n")
  boot_result <- boot::boot(
    data = modcall$data,
    statistic = boot_statistic,
    R = B,
    modcall = modcall,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    fallback = fallback_vec,
    resample_fn = resample_fn,
    show_progress = showProgress,
    parallel = parallel,
    ncpus = ncpus,
    ...
  )
  if (showProgress) message("done.\n")

  # Handle failed replicates
  boot_result <- handle_failed_replicates(boot_result, B)

  # Build parallel apply function for CI computation
  apply_fn <- build_parallel_apply(parallel, ncpus, cl)

  # Compute and attach CIs
  message("calculating bootstrap CIs....")
  result <- attach_bootstrap_cis(
    coef_template = coef_template,
    boot_result = boot_result,
    conf = conf,
    type = type,
    apply_fn = apply_fn
  )

  structure(result, call = match.call())
}


#' Simple bootstrap CIs for pffr (deprecated)
#'
#' @description
#' **Deprecated**
#'
#' `coefboot.pffr()` was renamed to [pffr_coefboot()] for consistency with the
#' package naming conventions.
#'
#' @inheritParams pffr_coefboot
#' @export
#' @keywords internal
coefboot.pffr <- function(
  object,
  n1 = 100,
  n2 = 40,
  n3 = 20,
  B = 100,
  ncpus = getOption("boot.ncpus", 1),
  parallel = c("no", "multicore", "snow"),
  cl = NULL,
  conf = c(.9, .95),
  type = "percent",
  method = c("resample", "residual", "residual.c"),
  showProgress = TRUE,
  ...
) {
  .Deprecated("pffr_coefboot")
  pffr_coefboot(
    object = object,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    B = B,
    ncpus = ncpus,
    parallel = parallel,
    cl = cl,
    conf = conf,
    type = type,
    method = method,
    showProgress = showProgress,
    ...
  )
}


# Prepare model call for bootstrap resampling
# @param object Fitted pffr model
# @returns Modified modcall with evaluated components
prepare_modcall_for_bootstrap <- function(object) {
  modcall <- object$pffr$call
  modcall$formula <- object$pffr$formula
  frml_env <- environment(object$pffr$formula)

  if (is.language(modcall$data)) {
    modcall$data <- eval(modcall$data, frml_env)
  }

  if (!is.null(modcall$ydata) && is.language(modcall$ydata)) {
    modcall$ydata <- eval(modcall$ydata, frml_env)
  }
  if (is.language(modcall$yind)) {
    modcall$yind <- eval(modcall$yind, frml_env)
  }
  modcall$fit <- TRUE
  modcall
}

# Build resampling function based on data type and method
# @param object Fitted pffr model
# @param modcall Prepared model call
# @param method Resampling method
# @returns Function(modcall, data, indices) for resampling
build_resample_function <- function(object, modcall, method) {
  is_sparse <- object$pffr$is_sparse
  ydata <- modcall$ydata

  if (method == "resample") {
    if (is_sparse) {
      # Closure captures ydata for sparse resampling
      function(mc, data, indices) {
        resample_sparse_observations(mc, data, ydata, indices)
      }
    } else {
      resample_grid_observations
    }
  } else if (method %in% c("residual", "residual.c")) {
    if (is_sparse) {
      stop("residual bootstrap not implemented for sparse data")
    }
    # Capture fitted values and residuals from original model
    fitted_vals <- fitted(object)
    resids <- residuals(object)
    response_name <- deparse(object$formula[[2]])
    center <- (method == "residual.c")

    function(mc, data, indices) {
      resample_residuals(
        mc,
        data,
        indices,
        fitted_vals,
        resids,
        response_name,
        center
      )
    }
  }
}

# Handle failed bootstrap replicates
# @param boot_result Object from boot::boot
# @param B Total number of replicates requested
# @returns Modified boot_result with failed replicates removed
handle_failed_replicates <- function(boot_result, B) {
  failed_idx <- which(rowSums(is.na(boot_result$t)) == ncol(boot_result$t))
  n_failed <- length(failed_idx)

  if (n_failed == B) {
    stop("All bootstrap replicates failed.")
  }
  if (n_failed > 0) {
    warning(n_failed, " out of ", B, " bootstrap replicates failed!")
    boot_result$t <- boot_result$t[-failed_idx, , drop = FALSE]
    boot_result$R <- boot_result$R - n_failed
  }
  boot_result
}

# Build parallel apply function for CI computation
# @param parallel Parallelization method
# @param ncpus Number of CPUs
# @param cl Optional cluster object
# @returns Apply function with appropriate parallelization
build_parallel_apply <- function(parallel, ncpus, cl) {
  if (parallel == "multicore") {
    parallel::mclapply
  } else if (parallel == "snow") {
    if (is.null(cl)) {
      cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
    }
    if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
      parallel::clusterSetRNGStream(cl)
    }
    # Return wrapper that manages cluster lifecycle
    function(X, FUN) {
      res <- parallel::parLapply(cl, X, FUN)
      parallel::stopCluster(cl)
      res
    }
  } else {
    lapply
  }
}

# Attach bootstrap CIs to coefficient template
# @param coef_template Coefficient structure from coef.pffr
# @param boot_result Bootstrap results
# @param conf Confidence levels
# @param type CI type
# @param apply_fn Parallel apply function
# @returns coef_template with CIs attached
attach_bootstrap_cis <- function(
  coef_template,
  boot_result,
  conf,
  type,
  apply_fn
) {
  # Get term structure
  pterm_names <- rownames(coef_template$pterms)
  n_pterms <- length(pterm_names)

  smterm_names <- names(coef_template$smterms)
  smterm_lengths <- vapply(
    coef_template$smterms,
    \(x) length(x$value),
    integer(1)
  )

  # Remove empty smooth terms
  nonempty <- smterm_lengths > 0
  if (any(!nonempty)) {
    smterm_names <- smterm_names[nonempty]
    smterm_lengths <- smterm_lengths[nonempty]
    coef_template$smterms <- coef_template$smterms[nonempty]
  }

  # Compute start/end indices for smooth terms (offset by n_pterms)
  smterm_starts <- c(1, cumsum(head(smterm_lengths, -1)) + 1)
  smterm_ends <- cumsum(smterm_lengths)

  # Attach CIs to parametric terms
  pterm_cis <- compute_boot_cis(
    boot_result,
    seq_len(n_pterms),
    conf,
    type,
    apply_fn
  )
  coef_template$pterms <- cbind(coef_template$pterms, pterm_cis)

  # Attach CIs to smooth terms
  for (i in seq_along(smterm_starts)) {
    indices <- n_pterms + (smterm_starts[i]:smterm_ends[i])
    smterm_cis <- compute_boot_cis(boot_result, indices, conf, type, apply_fn)
    coef_template$smterms[[i]] <- cbind(
      coef_template$smterms[[i]]$coef,
      smterm_cis
    )
  }
  names(coef_template$smterms) <- smterm_names

  coef_template
}

# =============================================================================
# GLS Decorrelation Helpers
# =============================================================================

# Compute square root of inverse covariance matrix for GLS decorrelation
# @param hat_sigma Estimated covariance matrix
# @param cond_cutoff Condition number cutoff for positive definiteness correction
# @returns Square root of inverse covariance matrix
compute_sqrt_sigma_inv <- function(hat_sigma, cond_cutoff = 500) {
  e_sigma <- eigen(hat_sigma, symmetric = TRUE)
  cond <- max(e_sigma$values) / min(e_sigma$values)

  if (cond > cond_cutoff) {
    hat_sigma_pd <- Matrix::nearPD(
      hat_sigma,
      keepDiag = TRUE,
      ensureSymmetry = TRUE,
      do2eigen = TRUE,
      posd.tol = 1 / cond_cutoff
    )
    warning(
      "Supplied <hatSigma> had condition number ",
      round(cond),
      "\n   -- projected further into pos.-definite cone (new condition number: ",
      cond_cutoff,
      ")."
    )
    e_sigma <- eigen(as.matrix(hat_sigma_pd$mat), symmetric = TRUE)
  }

  e_sigma$vectors %*% diag(1 / sqrt(e_sigma$values)) %*% t(e_sigma$vectors)
}

# Decorrelate design matrix and response for GLS fitting
# @param G GAM G object from gam(..., fit = FALSE)
# @param sqrt_sigma_inv Square root of inverse covariance
# @param nobs Number of observations
# @param nyindex Number of time points per observation
# @returns Modified G object with decorrelated X and y
decorrelate_gam_matrices <- function(G, sqrt_sigma_inv, nobs, nyindex) {
  # Decorrelate response (reshape, transform, flatten)
  y_matrix <- matrix(G$y, nrow = nobs, ncol = nyindex, byrow = TRUE)
  y_decorr <- t(sqrt_sigma_inv %*% t(y_matrix))
  G$y <- as.vector(t(y_decorr))

  # Decorrelate design matrix block-wise
  for (i in seq_len(nobs)) {
    rows <- ((i - 1) * nyindex + 1):(i * nyindex)
    G$X[rows, ] <- sqrt_sigma_inv %*% G$X[rows, ]
  }

  G
}

# =============================================================================
# pffr Label Map Building Helpers
# =============================================================================

# Build label map connecting original terms to mgcv smooth labels
# Extracted from pffr/pffrGLS for reusability
# @param term_map Named vector of transformed term strings
# @param m_smooth List of smooth objects from fitted model
# @param where_specials List of term indices by type
# @param formula_env Environment with transformed variables
# @param ffpc_terms List of ffpc term objects (or NULL)
# @returns Named list mapping original terms to mgcv labels
build_pffr_label_map <- function(
  term_map,
  m_smooth,
  where_specials,
  formula_env,
  ffpc_terms = NULL
) {
  label_map <- as.list(term_map)
  smooth_labels <- vapply(m_smooth, \(x) x$label, character(1))

  has_par <- length(where_specials$par) > 0
  has_ffpc <- length(where_specials$ffpc) > 0

  if (has_par || has_ffpc) {
    # Handle parametric terms
    if (has_par) {
      for (w in where_specials$par) {
        var <- get(names(label_map)[w], envir = formula_env)
        if (is.factor(var)) {
          where <- vapply(m_smooth, \(x) x$by, character(1)) ==
            names(label_map)[w]
          label_map[[w]] <- vapply(m_smooth[where], \(x) x$label, character(1))
        } else {
          label_map[[w]] <- paste0(
            "s(",
            names(label_map)[w],
            "):",
            names(label_map)[w]
          )
        }
      }
    }

    # Handle ffpc terms
    if (has_ffpc && !is.null(ffpc_terms)) {
      for (ind in seq_along(where_specials$ffpc)) {
        w <- where_specials$ffpc[ind]
        where <- vapply(m_smooth, \(x) x$id %||% NA_character_, character(1)) ==
          ffpc_terms[[ind]]$id
        label_map[[w]] <- vapply(m_smooth[where], \(x) x$label, character(1))
      }
    }

    # Match remaining terms
    other_indices <- setdiff(
      seq_along(label_map),
      c(where_specials$par, where_specials$ffpc)
    )
    if (length(other_indices) > 0) {
      label_map[other_indices] <- smooth_labels[pmatch(
        vapply(label_map[other_indices], extract_smooth_label, character(1)),
        smooth_labels
      )]
    }
  } else {
    # Simple case: match all terms
    label_map[seq_along(label_map)] <- smooth_labels[pmatch(
      vapply(label_map, extract_smooth_label, character(1)),
      smooth_labels
    )]
  }

  # Fill in NA labels with original term map
  na_labels <- vapply(
    label_map,
    \(x) any(is.null(x)) || any(is.na(x)),
    logical(1)
  )
  if (any(na_labels)) {
    label_map[na_labels] <- term_map[na_labels]
  }

  label_map
}

# Extract smooth label from term string or call
# @param x Term string or evaluated call
# @returns Label string
extract_smooth_label <- function(x) {
  parsed <- parse(text = x)[[1]]
  if (length(parsed) != 1) {
    evaled <- eval(parsed)
    evaled$label
  } else {
    x
  }
}

#' Penalized function-on-function regression with non-i.i.d. residuals
#'
#' Implements additive regression for functional and scalar covariates and functional responses.
#' This function is a wrapper for \code{mgcv}'s \code{\link[mgcv]{gam}} and its siblings to fit models of the general form \cr
#' \eqn{Y_i(t) = \mu(t) + \int X_i(s)\beta(s,t)ds + f(z_{1i}, t) + f(z_{2i}) + z_{3i} \beta_3(t) + \dots  + E_i(t))}\cr
#' with a functional (but not necessarily continuous) response \eqn{Y(t)},
#' (optional) smooth intercept \eqn{\mu(t)}, (multiple) functional covariates \eqn{X(t)} and scalar covariates
#' \eqn{z_1}, \eqn{z_2}, etc. The residual functions \eqn{E_i(t) \sim GP(0, K(t,t'))} are assumed to be i.i.d.
#' realizations of a Gaussian process. An estimate of the covariance operator \eqn{K(t,t')} evaluated on \code{yind}
#' has to be supplied in the \code{hatSigma}-argument.
#'
#' @section Details:
#' Note that \code{hatSigma} has to be positive definite. If \code{hatSigma} is close to positive \emph{semi-}definite or badly conditioned,
#' estimated standard errors become unstable (typically much too small). \code{pffr_gls} will try to diagnose this and issue a warning.
#' The danger is especially big if the number of functional observations is smaller than the number of gridpoints
#' (i.e, \code{length(yind)}), since the raw covariance estimate will not have full rank.\cr
#' Please see \code{\link[refund]{pffr}} for details on model specification and
#' implementation. \cr THIS IS AN EXPERIMENTAL VERSION AND NOT WELL TESTED YET -- USE AT YOUR OWN RISK.
#'
#' @param formula a formula with special terms as for \code{\link[mgcv]{gam}}, with additional special terms \code{\link{ff}()} and \code{c()}. See \code{\link[refund]{pffr}}.
#' @param yind a vector with length equal to the number of columns of the matrix of functional responses giving the vector of evaluation points \eqn{(t_1, \dots ,t_{G})}.
#'   see \code{\link[refund]{pffr}}
#' @param data an (optional) \code{data.frame} or named list containing the data.
#' @param ydata an (optional) \code{data.frame} for sparse/irregular functional responses. See \code{\link[refund]{pffr}}.
#' @param algorithm the name of the function used to estimate the model. Defaults to \code{\link[mgcv]{gam}} if the matrix of functional responses has less than \code{2e5} data points
#' 	 and to \code{\link[mgcv]{bam}} if not. "gamm" (see \code{\link[mgcv]{gamm}}) and "gamm4" (see \code{\link[gamm4]{gamm4}}) are valid options as well.
#' @param hatSigma (an estimate of) the within-observation covariance (along the responses' index), evaluated at \code{yind}. See Details.
#' @param method See \code{\link[refund]{pffr}}
#' @param bs.yindex See \code{\link[refund]{pffr}}
#' @param bs.int See \code{\link[refund]{pffr}}
#' @param tensortype See \code{\link[refund]{pffr}}
#' @param cond.cutoff if the condition number of \code{hatSigma} is greater than this, \code{hatSigma} is
#'    made ``more'' positive-definite via \code{\link[Matrix]{nearPD}} to ensure a condition number equal to cond.cutoff. Defaults to 500.
#' @param ... additional arguments that are valid for \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}. See \code{\link[refund]{pffr}}.
#'
#' @return a fitted \code{pffr}-object, see \code{\link[refund]{pffr}}.
#' @seealso \code{\link[refund]{pffr}}, \code{\link[refund]{fpca.sc}}
#'
#' @examples
#' \dontrun{
#' # Simulate data with correlated residuals
#' set.seed(1234)
#' n <- 50
#' nygrid <- 30
#' yind <- seq(0, 1, length.out = nygrid)
#'
#' # Create covariate and response
#' xlin <- rnorm(n)
#' beta_t <- sin(2 * pi * yind)
#' Y <- outer(xlin, beta_t) + matrix(rnorm(n * nygrid, sd = 0.5), n, nygrid)
#' dat <- data.frame(Y = I(Y), xlin = xlin)
#'
#' # Fit initial model to estimate residual covariance
#' m1 <- pffr(Y ~ xlin, yind = yind, data = dat)
#' resid_mat <- residuals(m1)
#' hat_sigma <- cov(resid_mat)
#'
#' # Fit GLS model with estimated covariance
#' m2 <- pffr_gls(Y ~ xlin, yind = yind, data = dat, hatSigma = hat_sigma)
#' summary(m2)
#' }
#'
#' @export
#' @importFrom Matrix nearPD
#' @importFrom mgcv gam bam gamm gam.fit3
#' @importFrom gamm4 gamm4
#' @importFrom lme4 lmer
#' @importFrom stats terms.formula
#' @author Fabian Scheipl
pffr_gls <- function(
  formula,
  yind,
  hatSigma,
  data = NULL,
  ydata = NULL,
  algorithm = NA,
  method = "REML",
  tensortype = c("ti", "t2"),
  bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
  bs.int = list(bs = "ps", k = 20, m = c(2, 1)),
  cond.cutoff = 5e2,
  ...
) {
  call <- match.call()
  tensortype <- as.symbol(match.arg(tensortype))

  ## warn if any entries in ... are not arguments for gam/gam.fit3 or gamm4/lmer
  dots <- list(...)
  if (length(dots)) {
    validDots <- if (!is.na(algorithm) && algorithm == "gamm4") {
      c(names(formals(gamm4)), names(formals(lmer)))
    } else {
      c(names(formals(gam)), names(formals(gam.fit3)))
    }
    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if (length(notUsed))
      warning(
        "Arguments <",
        paste(notUsed, collapse = ", "),
        "> supplied but not used."
      )
  }

  is_sparse <- !is.null(ydata)
  if (is_sparse) {
    stopifnot(ncol(ydata) == 3)
    stopifnot(c(".obs", ".index", ".value") == colnames(ydata))
  }

  # Parse formula using modular parser (same as pffr)
  parsed <- parse_pffr_model_formula(formula, data, ydata)
  tf <- parsed$tf
  trmstrings <- parsed$trmstrings
  terms <- parsed$terms
  frmlenv <- parsed$frmlenv
  where.specials <- parsed$where.specials
  responsename <- parsed$responsename

  #start new formula
  newfrml <- paste(responsename, "~", sep = "")
  formula_env <- new.env()
  evalenv <- if ("data" %in% names(call)) eval.parent(call$data) else NULL

  if (is_sparse) {
    nobs <- length(unique(ydata$.obs))
    stopifnot(all(ydata$.obs %in% rownames(data)))
    stopifnot(all(ydata$.obs %in% 1:nobs))

    nobs.data <- nrow(as.matrix(data[[1]]))
    stopifnot(nobs == nobs.data)
    ntotal <- nrow(ydata)

    # For GLS with sparse data, yind must be provided
    if (missing(yind)) {
      stop("yind must be provided for pffrGLS with sparse data")
    }
    nyindex <- length(yind)
  } else {
    nobs <- nrow(eval(responsename, envir = evalenv, enclos = frmlenv))
    nyindex <- ncol(eval(responsename, envir = evalenv, enclos = frmlenv))
    ntotal <- nobs * nyindex
  }

  if (missing(algorithm) || is.na(algorithm)) {
    algorithm <- ifelse(ntotal > 1e5, "bam", "gam")
  }
  algorithm <- as.symbol(algorithm)
  if (as.character(algorithm) == "bam" && !("chunk.size" %in% names(call))) {
    call$chunk.size <- 10000
  }

  # GLS not implemented for gamm4
  if (as.character(algorithm) == "gamm4") {
    stop("pffrGLS not implemented for gamm4")
  }

  # Validate and compute decorrelation matrix
  stopifnot(
    is.matrix(hatSigma),
    nrow(hatSigma) == nyindex,
    ncol(hatSigma) == nyindex
  )
  sqrt_sigma_inv <- compute_sqrt_sigma_inv(hatSigma, cond.cutoff)

  if (!is_sparse) {
    # Handle yind for regular grid data
    if (missing(yind)) {
      if (length(c(where.specials$ff, where.specials$sff))) {
        if (length(where.specials$ff)) {
          ffcall <- expand.call(ff, as.call(terms[where.specials$ff][1])[[1]])
        } else {
          ffcall <- expand.call(sff, as.call(terms[where.specials$sff][1])[[1]])
        }
        if (!is.null(ffcall$yind)) {
          yind <- eval(ffcall$yind, envir = evalenv, enclos = frmlenv)
          yindname <- deparse(ffcall$yind)
        } else {
          yind <- 1:nyindex
          yindname <- "yindex"
        }
      } else {
        yind <- 1:nyindex
        yindname <- "yindex"
      }
    } else {
      if (is.symbol(substitute(yind)) || is.character(yind)) {
        yindname <- deparse(substitute(yind))
        if (!is.null(data) && !is.null(data[[yindname]])) {
          yind <- data[[yindname]]
        }
      } else {
        yindname <- "yindex"
      }
      stopifnot(is.vector(yind), is.numeric(yind), length(yind) == nyindex)
    }
    if (length(yindname) > 1) yindname <- "yindex"
    stopifnot(all.equal(order(yind), 1:nyindex))

    yindvec <- rep(yind, times = nobs)
    yindex_vec_name <- as.symbol(paste(yindname, ".vec", sep = ""))
    assign(x = deparse(yindex_vec_name), value = yindvec, envir = formula_env)

    # Store original response for later restoration
    original_response <- as.vector(t(eval(
      responsename,
      envir = evalenv,
      enclos = frmlenv
    )))

    # Assign response in long format to formula_env
    assign(
      x = deparse(responsename),
      value = original_response,
      envir = formula_env
    )

    missing_indices <- if (any(is.na(original_response))) {
      which(is.na(original_response))
    } else NULL

    if (!is.null(missing_indices)) {
      stop("pffrGLS does not support missing values in response")
    }

    obs_indices <- rep(1:nobs, each = nyindex)
  } else {
    # Sparse/irregular data
    yindname <- "yindex"
    yindvec <- ydata$.index
    yindex_vec_name <- as.symbol(paste(yindname, ".vec", sep = ""))
    assign(
      x = deparse(yindex_vec_name),
      value = ydata$.index,
      envir = formula_env
    )

    original_response <- ydata$.value
    assign(x = deparse(responsename), value = ydata$.value, envir = formula_env)

    missing_indices <- NULL
    obs_indices <- ydata$.obs

    # For sparse data, GLS decorrelation is more complex
    # We need to handle irregular grid points
    stop(
      "pffrGLS with sparse data not yet fully implemented - decorrelation requires interpolation"
    )
  }

  ##################################################################################
  # Modify formula terms using modular helpers (same as pffr)
  newtrmstrings <- attr(tf, "term.labels")

  # Transform intercept
  if (parsed$has_intercept) {
    int_result <- transform_intercept_term(yindex_vec_name, bs.int, yindname)
    intstring <- int_result$term_string

    newfrml <- paste(newfrml, intstring, sep = " ")
    addFint <- TRUE
  } else {
    newfrml <- paste(newfrml, "0", sep = "")
    addFint <- FALSE
    intstring <- NULL
  }

  # Transform c() terms
  if (length(where.specials$c)) {
    newtrmstrings[where.specials$c] <- sapply(
      trmstrings[where.specials$c],
      transform_c_term
    )
  }

  # Prep function-on-function terms
  if (length(c(where.specials$ff, where.specials$sff))) {
    ffterms <- lapply(
      terms[c(where.specials$ff, where.specials$sff)],
      function(x) eval(x, envir = evalenv, enclos = frmlenv)
    )

    newtrmstrings[c(where.specials$ff, where.specials$sff)] <- sapply(
      ffterms,
      function(x) safeDeparse(x$call)
    )

    # Make ff data and assign to formula_env (same logic as pffr)
    makeff <- function(x) {
      tmat <- matrix(yindvec, nrow = length(yindvec), ncol = length(x$xind))
      smat <- matrix(
        x$xind,
        nrow = length(yindvec),
        ncol = length(x$xind),
        byrow = TRUE
      )
      if (!is.null(x[["LX"]])) {
        LStacked <- x$LX[obs_indices, ]
      } else {
        LStacked <- x$L[obs_indices, ]
        XStacked <- x$X[obs_indices, ]
      }
      if (!is.null(x$limits)) {
        use <- x$limits(smat, tmat)
        LStacked <- LStacked * use

        # Find windows and reduce matrix size if possible
        windows <- compute_integration_windows(use)
        max_width <- max(windows[, 3])
        if (max_width < ncol(smat)) {
          eff_windows <- expand_windows_to_maxwidth(windows, ncol(smat))
          smat <- shift_and_shorten_matrix(smat, eff_windows)
          tmat <- shift_and_shorten_matrix(tmat, eff_windows)
          LStacked <- shift_and_shorten_matrix(LStacked, eff_windows)
          if (is.null(x$LX)) {
            XStacked <- shift_and_shorten_matrix(XStacked, eff_windows)
          }
        }
      }
      assign(x = x$yindname, value = tmat, envir = formula_env)
      assign(x = x$xindname, value = smat, envir = formula_env)
      assign(x = x$LXname, value = LStacked, envir = formula_env)
      if (is.null(x[["LX"]])) {
        assign(x = x$xname, value = XStacked, envir = formula_env)
      }
      invisible(NULL)
    }
    lapply(ffterms, makeff)
  } else {
    ffterms <- NULL
  }

  # Prep ffpc terms
  if (length(where.specials$ffpc)) {
    ffpcterms <- lapply(terms[where.specials$ffpc], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })

    lapply(ffpcterms, function(trm) {
      lapply(colnames(trm$data), function(nm) {
        assign(x = nm, value = trm$data[obs_indices, nm], envir = formula_env)
        invisible(NULL)
      })
      invisible(NULL)
    })

    getFfpcFormula <- function(trm) {
      frmls <- lapply(colnames(trm$data), function(pc) {
        arglist <- c(
          name = "s",
          x = as.symbol(yindex_vec_name),
          by = as.symbol(pc),
          id = trm$id,
          trm$splinepars
        )
        call <- do.call("call", arglist, envir = formula_env)
        call$x <- as.symbol(yindex_vec_name)
        call$by <- as.symbol(pc)
        safeDeparse(call)
      })
      return(paste(unlist(frmls), collapse = " + "))
    }
    newtrmstrings[where.specials$ffpc] <- sapply(ffpcterms, getFfpcFormula)
    ffpcterms <- lapply(ffpcterms, function(x) x[names(x) != "data"])
  } else {
    ffpcterms <- NULL
  }

  # Prep pcre terms
  if (length(where.specials$pcre)) {
    pcreterms <- lapply(terms[where.specials$pcre], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })
    lapply(pcreterms, function(trm) {
      if (!is_sparse && all(trm$yind == yind)) {
        lapply(colnames(trm$efunctions), function(nm) {
          assign(
            x = nm,
            value = trm$efunctions[rep(1:nyindex, times = nobs), nm],
            envir = formula_env
          )
          invisible(NULL)
        })
      } else {
        stopifnot(min(trm$yind) <= min(yind))
        stopifnot(max(trm$yind) >= max(yind))
        lapply(colnames(trm$efunctions), function(nm) {
          tmp <- approx(
            x = trm$yind,
            y = trm$efunctions[, nm],
            xout = yindvec,
            method = "linear"
          )$y
          assign(x = nm, value = tmp, envir = formula_env)
          invisible(NULL)
        })
      }
      assign(x = trm$idname, value = trm$id[obs_indices], envir = formula_env)
      invisible(NULL)
    })
    newtrmstrings[where.specials$pcre] <- sapply(
      pcreterms,
      function(x) safeDeparse(x$call)
    )
  } else {
    pcreterms <- NULL
  }

  # Transform smooth terms (s, te, t2)
  if (length(c(where.specials$s, where.specials$te, where.specials$t2))) {
    newtrmstrings[c(where.specials$s, where.specials$te, where.specials$t2)] <-
      sapply(
        terms[c(where.specials$s, where.specials$te, where.specials$t2)],
        function(x)
          transform_smooth_term(
            x,
            yindex_vec_name,
            bs.yindex,
            tensortype,
            algorithm
          )
      )
  }

  # Transform parametric terms
  if (length(where.specials$par)) {
    newtrmstrings[where.specials$par] <- sapply(
      terms[where.specials$par],
      function(x) transform_par_term(x, yindex_vec_name, bs.yindex)
    )
  }

  # Assign expanded variables to formula_env
  where.specials$notff <- c(
    where.specials$c,
    where.specials$par,
    where.specials$s,
    where.specials$te,
    where.specials$t2
  )
  if (length(where.specials$notff)) {
    evalenv <- if ("data" %in% names(call)) {
      list2env(eval.parent(call$data))
    } else frmlenv
    lapply(terms[where.specials$notff], function(x) {
      isC <- safeDeparse(x) %in% sapply(terms[where.specials$c], safeDeparse)
      if (isC) {
        x <- formula(paste(
          "~",
          gsub("\\)$", "", gsub("^c\\(", "", deparse(x)))
        ))[[2]]
      }
      nms <- if (!is.null(names(x))) {
        all.vars(x[names(x) %in% c("", "by")])
      } else all.vars(x)

      sapply(nms, function(nm) {
        var <- get(nm, envir = evalenv)
        if (is.matrix(var)) {
          stopifnot(!is_sparse || ncol(var) == nyindex)
          assign(x = nm, value = as.vector(t(var)), envir = formula_env)
        } else {
          stopifnot(length(var) == nobs)
          assign(x = nm, value = var[obs_indices], envir = formula_env)
        }
        invisible(NULL)
      })
      invisible(NULL)
    })
  }

  # Build mgcv formula
  newfrml <- build_mgcv_formula(
    responsename = responsename,
    intercept_string = if (addFint) intstring else NULL,
    term_strings = newtrmstrings,
    has_intercept = addFint,
    formula_env = formula_env
  )

  # Build mgcv data
  pffrdata <- build_mgcv_data(formula_env)

  # Build the call for gam/bam
  newcall <- expand.call(pffr, call)
  newcall$yind <- newcall$tensortype <- newcall$bs.int <-
    newcall$bs.yindex <- newcall$algorithm <- newcall$ydata <-
      newcall$hatSigma <- newcall$cond.cutoff <- NULL
  newcall$formula <- newfrml
  newcall$data <- quote(pffrdata)
  newcall[[1]] <- algorithm

  # Handle ... args
  dotargs <- names(newcall)[names(newcall) %in% names(dots)]
  newcall[dotargs] <- dots[dotargs]

  if ("weights" %in% dotargs) {
    if (length(dots$weights) == nobs) {
      newcall$weights <- dots$weights[obs_indices]
    } else if (
      !is.null(dim(dots$weights)) && all(dim(dots$weights) == c(nobs, nyindex))
    ) {
      newcall$weights <- as.vector(t(dots$weights))
    } else {
      stop(
        "weights must be vector of length nobs or matrix of dim (nobs, nyindex)"
      )
    }
  }

  if ("offset" %in% dotargs) {
    if (length(dots$offset) == nobs) {
      newcall$offset <- dots$offset[obs_indices]
    } else if (
      !is.null(dim(dots$offset)) && all(dim(dots$offset) == c(nobs, nyindex))
    ) {
      newcall$offset <- as.vector(t(dots$offset))
    } else {
      stop(
        "offset must be vector of length nobs or matrix of dim (nobs, nyindex)"
      )
    }
  }

  # GLS fitting: get G object, decorrelate, refit
  newcall$fit <- FALSE
  G <- eval(newcall)
  G <- decorrelate_gam_matrices(G, sqrt_sigma_inv, nobs, nyindex)
  m <- gam(G = G, fit = TRUE, method = method)

  ###--- Post-fit corrections: restore original scale ---###
  m$y <- original_response
  m$fitted.values <- predict(m)
  m$residuals <- m$y - m$fitted.values
  ###--- End post-fit corrections ---###

  # Build term mapping (same as pffr)
  m.smooth <- if (as.character(algorithm) %in% c("gamm4", "gamm")) {
    m$gam$smooth
  } else {
    m$smooth
  }

  trmmap <- newtrmstrings
  names(trmmap) <- names(terms)
  if (addFint) trmmap <- c(trmmap, intstring)

  # Build label map (same as pffr)
  labelmap <- as.list(trmmap)
  lbls <- sapply(m.smooth, function(x) x$label)
  if (length(c(where.specials$par, where.specials$ffpc))) {
    if (length(where.specials$par)) {
      for (w in where.specials$par) {
        if (is.factor(get(names(labelmap)[w], envir = formula_env))) {
          labelmap[[w]] <- {
            where <- sapply(m.smooth, function(x) x$by) == names(labelmap)[w]
            sapply(m.smooth[where], function(x) x$label)
          }
        } else {
          labelmap[[w]] <- paste0(
            "s(",
            yindex_vec_name,
            "):",
            names(labelmap)[w]
          )
        }
      }
    }
    if (length(where.specials$ffpc)) {
      ind <- 1
      for (w in where.specials$ffpc) {
        labelmap[[w]] <- {
          where <- sapply(m.smooth, function(x) x$id) == ffpcterms[[ind]]$id
          sapply(m.smooth[where], function(x) x$label)
        }
        ind <- ind + 1
      }
    }
    labelmap[-c(where.specials$par, where.specials$ffpc)] <- lbls[pmatch(
      sapply(
        labelmap[-c(where.specials$par, where.specials$ffpc)],
        function(x) {
          if (length(parse(text = x)[[1]]) != 1) {
            tmp <- eval(parse(text = x))
            return(tmp$label)
          } else {
            return(x)
          }
        }
      ),
      lbls
    )]
  } else {
    labelmap[1:length(labelmap)] <- lbls[pmatch(
      sapply(labelmap[1:length(labelmap)], function(x) {
        if (length(parse(text = x)[[1]]) != 1) {
          tmp <- eval(parse(text = x))
          return(tmp$label)
        } else {
          return(x)
        }
      }),
      lbls
    )]
  }
  nalbls <- sapply(
    labelmap,
    function(x) any(is.null(x)) | any(is.na(x[!is.null(x)]))
  )
  if (any(nalbls)) {
    labelmap[nalbls] <- trmmap[nalbls]
  }

  names(m.smooth) <- lbls
  if (as.character(algorithm) %in% c("gamm4", "gamm")) {
    m$gam$smooth <- m.smooth
  } else {
    m$smooth <- m.smooth
  }

  # Create shortlabels mapping (same as pffr)
  shortlabels <- create_shortlabels(
    labelmap = labelmap,
    m.smooth = m.smooth,
    yindname = yindname,
    where.specials = where.specials
  )

  # Build return structure (complete $pffr slot)
  ret <- list(
    call = call,
    formula = formula,
    termmap = trmmap,
    labelmap = labelmap,
    shortlabels = shortlabels,
    responsename = responsename,
    nobs = nobs,
    nyindex = nyindex,
    yindname = yindname,
    yind = yind,
    where = where.specials,
    ff = ffterms,
    ffpc = ffpcterms,
    pcreterms = pcreterms,
    missing_indices = missing_indices,
    is_sparse = is_sparse,
    ydata = ydata,
    # GLS-specific
    hatSigma = hatSigma,
    sqrtSigmaInv = sqrt_sigma_inv
  )

  if (as.character(algorithm) %in% c("gamm4", "gamm")) {
    m$gam$pffr <- ret
    class(m$gam) <- c("pffr", class(m$gam))
  } else {
    m$pffr <- ret
    class(m) <- c("pffr", class(m))
  }

  return(m)
} # end pffr_gls()


#' Penalized function-on-function regression with non-i.i.d. residuals (deprecated)
#'
#' @description
#' **Deprecated**
#'
#' `pffrGLS()` was renamed to [pffr_gls()] for consistency with the
#' package naming conventions.
#'
#' @inheritParams pffr_gls
#' @export
#' @keywords internal
pffrGLS <- function(
  formula,
  yind,
  hatSigma,
  data = NULL,
  ydata = NULL,
  algorithm = NA,
  method = "REML",
  tensortype = c("ti", "t2"),
  bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
  bs.int = list(bs = "ps", k = 20, m = c(2, 1)),
  cond.cutoff = 5e2,
  ...
) {
  .Deprecated("pffr_gls")
  pffr_gls(
    formula = formula,
    yind = yind,
    hatSigma = hatSigma,
    data = data,
    ydata = ydata,
    algorithm = algorithm,
    method = method,
    tensortype = tensortype,
    bs.yindex = bs.yindex,
    bs.int = bs.int,
    cond.cutoff = cond.cutoff,
    ...
  )
}
