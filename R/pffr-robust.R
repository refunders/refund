#--------------------------------------
# Bootstrap Resampling Helpers (pure functions)
#--------------------------------------

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

#--------------------------------------
# Bootstrap CI Computation Helpers
#--------------------------------------

# Build flatten/reconstruction layout for coef.pffr-like objects
# @param coef_template Object returned by coef(..., ci = "none")
# @returns List with coefficient index/range metadata
build_coefboot_layout <- function(coef_template) {
  n_pterms <- nrow(as.matrix(coef_template$pterms))
  pterm_names <- rownames(coef_template$pterms)
  if (is.null(pterm_names)) pterm_names <- paste0("pterm_", seq_len(n_pterms))

  smterm_names <- names(coef_template$smterms)
  smterm_lengths <- vapply(
    coef_template$smterms,
    \(sm) nrow(as.data.frame(sm$coef)),
    integer(1)
  )

  smterm_ranges <- vector("list", length(smterm_lengths))
  next_idx <- n_pterms + 1L
  for (i in seq_along(smterm_lengths)) {
    len_i <- smterm_lengths[[i]]
    if (len_i > 0L) {
      smterm_ranges[[i]] <- next_idx:(next_idx + len_i - 1L)
      next_idx <- next_idx + len_i
    } else {
      smterm_ranges[[i]] <- integer(0)
    }
  }

  coef_names <- character(next_idx - 1L)
  if (n_pterms > 0L) {
    coef_names[seq_len(n_pterms)] <- pterm_names
  }
  for (i in seq_along(smterm_ranges)) {
    idx <- smterm_ranges[[i]]
    if (!length(idx)) next
    coef_names[idx] <- paste0(smterm_names[[i]], "[", seq_along(idx), "]")
  }

  list(
    n_pterms = n_pterms,
    pterm_names = pterm_names,
    smterm_names = smterm_names,
    smterm_lengths = smterm_lengths,
    smterm_ranges = smterm_ranges,
    n_coefs = length(coef_names),
    coef_names = coef_names
  )
}

# Flatten coef.pffr-like object to a numeric vector with stable ordering
# @param coef_object Object returned by coef(..., ci = "none")
# @param layout Layout returned by build_coefboot_layout()
# @returns Numeric coefficient vector
flatten_coefboot_coefficients <- function(coef_object, layout) {
  out <- numeric(layout$n_coefs)

  pvals <- as.numeric(coef_object$pterms[, "value", drop = TRUE])
  if (length(pvals) != layout$n_pterms) {
    stop("Mismatch in number of parametric coefficients.")
  }
  if (layout$n_pterms > 0L) {
    out[seq_len(layout$n_pterms)] <- pvals
  }

  for (i in seq_along(layout$smterm_ranges)) {
    idx <- layout$smterm_ranges[[i]]
    if (!length(idx)) next
    sval <- as.numeric(coef_object$smterms[[i]]$coef[, "value", drop = TRUE])
    if (length(sval) != length(idx)) {
      stop("Mismatch in smooth-term coefficient length for term ", i, ".")
    }
    out[idx] <- sval
  }

  out
}

# Build fixed evaluation grids for coef.pffr extraction
# @param object Fitted pffr model
# @param n1,n2,n3 Grid sizes
# @returns Named list of per-smooth data grids
build_coefboot_eval_grid <- function(object, n1, n2, n3) {
  pffr_info <- list(
    yind_name = object$pffr$yind_name,
    pcre_terms = object$pffr$pcre_terms
  )
  grid_sizes <- list(n1 = n1, n2 = n2, n3 = n3)

  grids <- lapply(seq_along(object$smooth), \(i) {
    trm <- object$smooth[[i]]
    is_pcre <- "pcre.random.effect" %in% class(trm)
    if (trm$dim > 3 && !is_pcre) return(NULL)
    coef_make_data_grid(
      trm = trm,
      model_data = object$model,
      pffr_info = pffr_info,
      grid_sizes = grid_sizes,
      is_pcre = is_pcre
    )
  })
  names(grids) <- names(object$smooth)
  grids
}

# Extract coefficients from a fitted pffr model as a vector
# @param model Fitted pffr model (or try-error)
# @param n1, n2, n3 Grid sizes for coefficient evaluation
# @param fallback Vector of NAs to return on failure
# @param layout Layout from build_coefboot_layout()
# @param eval_grid Optional fixed evaluation grids for coef.pffr
# @returns Numeric vector of coefficient values
extract_boot_coefficients <- function(
  model,
  n1,
  n2,
  n3,
  fallback,
  layout,
  eval_grid = NULL
) {
  if (inherits(model, "try-error")) {
    return(fallback)
  }

  coefs <- try(
    coef(
      model,
      se = FALSE,
      sandwich = "none",
      ci = "none",
      n1 = n1,
      n2 = n2,
      n3 = n3,
      eval_grid = eval_grid
    ),
    silent = TRUE
  )
  if (inherits(coefs, "try-error")) {
    return(fallback)
  }

  vec <- try(flatten_coefboot_coefficients(coefs, layout), silent = TRUE)
  if (inherits(vec, "try-error") || length(vec) != length(fallback)) {
    return(fallback)
  }
  vec
}

# Convert confidence level to stable column suffix (e.g., 0.95 -> "95")
# @param conf Numeric confidence level in (0, 1)
# @returns Character label
format_boot_conf_label <- function(conf) {
  if (!is.numeric(conf) || length(conf) != 1L || !is.finite(conf)) {
    stop("'conf' must be a single finite numeric value.")
  }
  label <- sprintf("%.8f", 100 * conf)
  label <- sub("0+$", "", label)
  label <- sub("\\.$", "", label)
  label
}

# Compute percentile CI bounds from a bootstrap sample vector
# @param x Numeric bootstrap draws
# @param conf Confidence levels
# @returns Numeric vector lower/upper pairs for each level
compute_percentile_ci <- function(x, conf) {
  as.vector(vapply(conf, \(c_level) {
    probs <- c((1 - c_level) / 2, 1 - (1 - c_level) / 2)
    stats::quantile(
      x,
      probs = probs,
      type = 8,
      na.rm = TRUE,
      names = FALSE
    )
  }, numeric(2)))
}

# Compute bootstrap confidence intervals for a single coefficient
# @param boot_result Object returned by boot::boot
# @param index Which coefficient index
# @param conf Confidence levels (numeric vector)
# @param type Type of CI ("percent", "bca", etc.)
# @returns Numeric vector of CI bounds
compute_single_boot_ci <- function(boot_result, index, conf, type) {
  if (type == "percent") {
    return(compute_percentile_ci(boot_result$t[, index], conf))
  }

  type_arg <- if (type == "percent") "perc" else type
  out <- tryCatch(
    {
      ci <- boot::boot.ci(
        boot_result,
        conf = conf,
        index = index,
        type = type_arg
      )
      bounds <- ci[[type]]
      if (is.null(bounds)) stop("CI bounds unavailable")
      lower <- bounds[, ncol(bounds) - 1L]
      upper <- bounds[, ncol(bounds)]
      as.vector(rbind(lower, upper))
    },
    warning = function(w) rep(NA_real_, 2L * length(conf)),
    error = function(e) rep(NA_real_, 2L * length(conf))
  )
  out
}

# Compute bootstrap CIs for multiple coefficients
# @param boot_result Object returned by boot::boot
# @param indices Which coefficient indices to compute CIs for
# @param conf Confidence levels
# @param type Type of CI
# @param apply_fn Function to use for iteration
# @returns Matrix with CI bounds as rows
compute_boot_cis <- function(
  boot_result,
  indices,
  conf,
  type,
  apply_fn = function(X, FUN) lapply(X, FUN)
) {
  ci_list <- apply_fn(indices, \(i) {
    compute_single_boot_ci(boot_result, i, conf, type)
  })
  ci_matrix <- do.call(rbind, ci_list)

  lower_names <- paste0("lower_", vapply(conf, format_boot_conf_label, character(1)))
  upper_names <- paste0("upper_", vapply(conf, format_boot_conf_label, character(1)))
  colnames(ci_matrix) <- as.vector(rbind(lower_names, upper_names))

  ci_matrix
}

#--------------------------------------
# Main Exported Functions
#--------------------------------------

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
  type <- match.arg(type, c("norm", "basic", "stud", "percent", "bca"))

  if (!is.numeric(B) || length(B) != 1 || !is.finite(B) || B < 1) {
    stop("'B' must be a single positive integer.")
  }
  B <- as.integer(B)
  if (!is.numeric(conf) || any(!is.finite(conf)) || any(conf <= 0 | conf >= 1)) {
    stop("'conf' must contain values strictly between 0 and 1.")
  }
  conf <- as.numeric(conf)
  conf <- conf[!duplicated(conf)]
  primary_conf <- conf[1L]

  if (is.null(ncpus)) {
    ncpus <- getOption("boot.ncpus", 1L)
  }
  if (!is.numeric(ncpus) || length(ncpus) != 1 || !is.finite(ncpus) || ncpus < 1) {
    stop("'ncpus' must be a single positive integer.")
  }
  ncpus <- as.integer(ncpus)

  # Extract model call components, evaluating any unevaluated language objects
  modcall <- prepare_modcall_for_bootstrap(object)

  # Get original coefficients (template for return structure)
  coef_template <- coef(
    object,
    se = FALSE,
    sandwich = "none",
    ci = "none",
    n1 = n1,
    n2 = n2,
    n3 = n3
  )
  layout <- build_coefboot_layout(coef_template)
  eval_grid <- build_coefboot_eval_grid(object, n1, n2, n3)

  # Compute fallback return vector (all NAs) for failed bootstrap replicates
  fallback_vec <- rep(NA_real_, layout$n_coefs)

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
    layout,
    eval_grid,
    resample_fn,
    show_progress
  ) {
    modcall_boot <- resample_fn(modcall, data, indices)
    if (!is.call(modcall_boot)) {
      modcall_boot <- as.call(modcall_boot)
    }
    model_boot <- try(eval(modcall_boot), silent = TRUE)
    coefs <- extract_boot_coefficients(
      model_boot,
      n1,
      n2,
      n3,
      fallback,
      layout,
      eval_grid = eval_grid
    )
    if (show_progress) cat(".")
    coefs
  }

  # Run bootstrap
  if (showProgress) message("starting bootstrap (", B, " replications).\n")
  boot_args <- list(
    data = modcall$data,
    statistic = boot_statistic,
    R = B,
    modcall = modcall,
    n1 = n1,
    n2 = n2,
    n3 = n3,
    fallback = fallback_vec,
    layout = layout,
    eval_grid = eval_grid,
    resample_fn = resample_fn,
    show_progress = showProgress,
    parallel = parallel,
    ncpus = ncpus,
    ...
  )
  if (parallel == "snow" && !is.null(cl)) {
    boot_args$cl <- cl
  }
  boot_result <- do.call(boot::boot, boot_args, quote = TRUE)
  if (showProgress) message("done.\n")

  # Handle failed replicates
  boot_result <- handle_failed_replicates(boot_result, B)
  n_failed <- attr(boot_result, "n_failed") %||% 0L
  failure_rate <- attr(boot_result, "failure_rate") %||%
    if (B > 0L) n_failed / B else NA_real_

  # Build parallel apply function for CI computation
  parallel_ctx <- build_parallel_apply(parallel, ncpus, cl)
  on.exit(parallel_ctx$cleanup(), add = TRUE)
  apply_fn <- parallel_ctx$apply

  # Compute and attach CIs
  if (showProgress) message("calculating bootstrap CIs....")
  result <- attach_bootstrap_cis(
    coef_template = coef_template,
    layout = layout,
    boot_result = boot_result,
    conf = conf,
    primary_conf = primary_conf,
    type = type,
    apply_fn = apply_fn,
    method = method,
    B_requested = B,
    n_failed = n_failed,
    failure_rate = failure_rate
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
  # Ensure bootstrap refits do not trigger robust covariance recalculation.
  fit_fun <- safeDeparse(modcall[[1]])
  if (fit_fun %in% c("pffr", "refund::pffr")) {
    modcall$sandwich <- "none"
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
    gaussian_identity <- identical(object$family$family, "gaussian") &&
      identical(object$family$link[[1]], "identity")
    if (!gaussian_identity) {
      stop(
        "residual bootstrap is only implemented for gaussian identity models. ",
        "Use method = \"resample\"."
      )
    }
    # Capture fitted values and residuals from original model
    fitted_vals <- fitted(object)
    resids <- residuals(object, reformat = TRUE)
    response_name <- deparse(object$pffr$response_name)
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
  failure_rate <- if (B > 0L) n_failed / B else NA_real_

  if (n_failed == B) {
    stop("All bootstrap replicates failed.")
  }
  if (n_failed > 0) {
    warning(
      n_failed, " out of ", B, " bootstrap replicates failed (",
      round(100 * failure_rate, 1), "%). ",
      "Using ", B - n_failed, " successful replicates for CI computation. ",
      "Consider increasing B or simplifying the model if this persists.",
      call. = FALSE
    )
    boot_result$t <- boot_result$t[-failed_idx, , drop = FALSE]
    boot_result$R <- boot_result$R - n_failed
  }
  attr(boot_result, "n_failed") <- n_failed
  attr(boot_result, "B_requested") <- B
  attr(boot_result, "failure_rate") <- failure_rate
  boot_result
}

# Build parallel apply function for CI computation
# @param parallel Parallelization method
# @param ncpus Number of CPUs
# @param cl Optional cluster object
# @returns Apply function with appropriate parallelization
build_parallel_apply <- function(parallel, ncpus, cl) {
  if (parallel == "multicore") {
    return(list(
      apply = function(X, FUN) parallel::mclapply(X, FUN, mc.cores = ncpus),
      cleanup = function() invisible(NULL)
    ))
  }

  if (parallel == "snow") {
    own_cluster <- is.null(cl)
    cl_use <- cl
    if (own_cluster) {
      cl_use <- parallel::makePSOCKcluster(rep("localhost", ncpus))
      if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
        parallel::clusterSetRNGStream(cl_use)
      }
    }
    return(list(
      apply = function(X, FUN) parallel::parLapply(cl_use, X, FUN),
      cleanup = function() {
        if (own_cluster) parallel::stopCluster(cl_use)
        invisible(NULL)
      }
    ))
  }

  list(
    apply = function(X, FUN) lapply(X, FUN),
    cleanup = function() invisible(NULL)
  )
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
  layout,
  boot_result,
  conf,
  primary_conf,
  type,
  apply_fn,
  method,
  B_requested,
  n_failed,
  failure_rate
) {
  ci_matrix <- compute_boot_cis(
    boot_result = boot_result,
    indices = seq_len(layout$n_coefs),
    conf = conf,
    type = type,
    apply_fn = apply_fn
  )
  rownames(ci_matrix) <- layout$coef_names

  primary_label <- format_boot_conf_label(primary_conf)
  lower_col <- paste0("lower_", primary_label)
  upper_col <- paste0("upper_", primary_label)

  p_idx <- seq_len(layout$n_pterms)
  if (layout$n_pterms > 0L) {
    coef_template$pterms <- cbind(
      coef_template$pterms,
      lower = ci_matrix[p_idx, lower_col],
      upper = ci_matrix[p_idx, upper_col]
    )
  }

  for (i in seq_along(layout$smterm_ranges)) {
    idx <- layout$smterm_ranges[[i]]
    if (!length(idx)) {
      coef_template$smterms[[i]]$coef <- cbind(
        coef_template$smterms[[i]]$coef,
        lower = numeric(0),
        upper = numeric(0)
      )
      next
    }
    coef_template$smterms[[i]]$coef <- cbind(
      coef_template$smterms[[i]]$coef,
      lower = ci_matrix[idx, lower_col],
      upper = ci_matrix[idx, upper_col]
    )
  }

  boot_ci <- list(
    type = type,
    conf = conf,
    primary_conf = primary_conf,
    coefficient_bounds = ci_matrix,
    pterms = if (layout$n_pterms > 0L) ci_matrix[p_idx, , drop = FALSE] else
      matrix(numeric(0), nrow = 0, ncol = ncol(ci_matrix)),
    smterms = lapply(layout$smterm_ranges, \(idx) {
      if (!length(idx)) {
        matrix(numeric(0), nrow = 0, ncol = ncol(ci_matrix))
      } else {
        ci_matrix[idx, , drop = FALSE]
      }
    })
  )
  names(boot_ci$smterms) <- layout$smterm_names

  coef_template$boot_ci <- boot_ci
  coef_template$ci_meta <- list(
    type = "bootstrap",
    level = primary_conf,
    levels = conf,
    boot_type = type,
    B_requested = B_requested,
    B_used = boot_result$R,
    n_failed = n_failed,
    failure_rate = failure_rate,
    method = method
  )

  coef_template
}

#--------------------------------------
# GLS Decorrelation Helpers
#--------------------------------------

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

#--------------------------------------
# pffr Label Map Building Helpers
#--------------------------------------

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
  yind_missing <- missing(yind)

  # Step 1: Validate arguments (no AR checking for GLS)
  dots <- list(...)
  validated <- pffr_validate_dots(call, algorithm, dots, check_ar = FALSE)
  dots <- validated$dots

  # Validate ydata format
  pffr_validate_ydata(ydata)
  is_sparse <- !is.null(ydata)

  # Step 2: Parse formula
  parsed <- parse_pffr_model_formula(formula, data, ydata)
  tf <- parsed$tf
  term_strings <- parsed$term_strings
  terms <- parsed$terms
  frml_env <- parsed$frml_env
  where_specials <- parsed$where_specials
  response_name <- parsed$response_name

  # Step 3: Detect data dimensions
  formula_env <- new.env()
  eval_env <- if ("data" %in% names(call)) eval.parent(call$data) else NULL

  if (is_sparse) {
    nobs <- length(unique(ydata$.obs))
    stopifnot(all(ydata$.obs %in% rownames(data)))
    stopifnot(all(ydata$.obs %in% 1:nobs))
    nobs_data <- nrow(as.matrix(data[[1]]))
    stopifnot(nobs == nobs_data)
    ntotal <- nrow(ydata)

    if (yind_missing) {
      stop("yind must be provided for pffr_gls with sparse data")
    }
    nyindex <- length(yind)
    yind_name <- "yindex"
  } else {
    nobs <- nrow(eval(response_name, envir = eval_env, enclos = frml_env))
    nyindex <- ncol(eval(response_name, envir = eval_env, enclos = frml_env))
    ntotal <- nobs * nyindex

    # Setup yind for dense (same logic as pffr)
    if (yind_missing) {
      if (length(c(where_specials$ff, where_specials$sff))) {
        if (length(where_specials$ff)) {
          ffcall <- expand.call(ff, as.call(terms[where_specials$ff][1])[[1]])
        } else {
          ffcall <- expand.call(sff, as.call(terms[where_specials$sff][1])[[1]])
        }
        if (!is.null(ffcall$yind)) {
          yind <- eval(ffcall$yind, envir = eval_env, enclos = frml_env)
          yind_name <- deparse(ffcall$yind)
        } else {
          yind <- 1:nyindex
          yind_name <- "yindex"
        }
      } else {
        yind <- 1:nyindex
        yind_name <- "yindex"
      }
    } else {
      if (is.symbol(substitute(yind)) | is.character(yind)) {
        yind_name <- deparse(substitute(yind))
        if (!is.null(data) && !is.null(data[[yind_name]])) {
          yind <- data[[yind_name]]
        }
      } else {
        yind_name <- "yindex"
      }
      stopifnot(is.vector(yind), is.numeric(yind), length(yind) == nyindex)
    }
    if (length(yind_name) > 1) yind_name <- "yindex"
    stopifnot(all.equal(order(yind), 1:nyindex))
  }

  # Step 4: Configure algorithm
  if (is.na(algorithm)) {
    algorithm <- ifelse(ntotal > 1e5, "bam", "gam")
  }
  algorithm <- as.symbol(algorithm)
  if (as.character(algorithm) == "bam" && !("chunk.size" %in% names(call))) {
    call$chunk.size <- 10000
  }
  if (as.character(algorithm) == "gamm4") {
    stop("pffr_gls not implemented for gamm4")
  }

  # Step 5: Validate and compute decorrelation matrix
  stopifnot(
    is.matrix(hatSigma),
    nrow(hatSigma) == nyindex,
    ncol(hatSigma) == nyindex
  )
  sqrt_sigma_inv <- compute_sqrt_sigma_inv(hatSigma, cond.cutoff)

  # Step 6: Setup response data
  if (is_sparse) {
    # Sparse data not fully supported
    yind_vec <- ydata$.index
    yindex_vec_name <- as.symbol(paste0(yind_name, ".vec"))
    assign(deparse(yindex_vec_name), yind_vec, envir = formula_env)

    original_response <- ydata$.value
    assign(deparse(response_name), ydata$.value, envir = formula_env)
    missing_indices <- NULL
    obs_indices <- ydata$.obs

    stop(
      "pffr_gls with sparse data not yet fully implemented - decorrelation requires interpolation"
    )
  } else {
    resp_info <- pffr_setup_response(
      is_sparse = FALSE,
      ydata = NULL,
      yind = yind,
      yind_name = yind_name,
      nobs = nobs,
      nyindex = nyindex,
      response_name = response_name,
      eval_env = eval_env,
      frml_env = frml_env,
      formula_env = formula_env
    )
    yind_vec <- resp_info$yind_vec
    yindex_vec_name <- resp_info$yindex_vec_name
    obs_indices <- resp_info$obs_indices
    missing_indices <- resp_info$missing_indices

    # Store original response for post-fit restoration
    original_response <- get(deparse(response_name), formula_env)

    if (!is.null(missing_indices)) {
      stop("pffr_gls does not support missing values in response")
    }
  }

  # Step 7: Transform formula terms (same as pffr)
  new_term_strings <- attr(tf, "term.labels")

  if (parsed$has_intercept) {
    int_result <- transform_intercept_term(yindex_vec_name, bs.int, yind_name)
    int_string <- int_result$term_string
    add_f_int <- TRUE
  } else {
    add_f_int <- FALSE
    int_string <- NULL
  }

  if (length(where_specials$c)) {
    new_term_strings[where_specials$c] <- sapply(
      term_strings[where_specials$c],
      transform_c_term
    )
  }

  # Process ff/sff terms
  if (length(c(where_specials$ff, where_specials$sff))) {
    ff_terms <- lapply(
      terms[c(where_specials$ff, where_specials$sff)],
      \(x) eval(x, envir = eval_env, enclos = frml_env)
    )
    new_term_strings[c(where_specials$ff, where_specials$sff)] <- sapply(
      ff_terms,
      \(x) safeDeparse(x$call)
    )

    makeff <- function(x) {
      tmat <- matrix(yind_vec, nrow = length(yind_vec), ncol = length(x$xind))
      smat <- matrix(
        x$xind,
        nrow = length(yind_vec),
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
        windows <- compute_integration_windows(use)
        max_width <- max(windows[, 3])
        if (max_width < ncol(smat)) {
          eff_windows <- expand_windows_to_maxwidth(windows, ncol(smat))
          smat <- shift_and_shorten_matrix(smat, eff_windows)
          tmat <- shift_and_shorten_matrix(tmat, eff_windows)
          LStacked <- shift_and_shorten_matrix(LStacked, eff_windows)
          if (is.null(x$LX))
            XStacked <- shift_and_shorten_matrix(XStacked, eff_windows)
        }
      }
      assign(x$yindname, tmat, envir = formula_env)
      assign(x$xindname, smat, envir = formula_env)
      assign(x$LXname, LStacked, envir = formula_env)
      if (is.null(x[["LX"]])) assign(x$xname, XStacked, envir = formula_env)
      invisible(NULL)
    }
    lapply(ff_terms, makeff)
  } else {
    ff_terms <- NULL
  }

  # Process ffpc terms
  if (length(where_specials$ffpc)) {
    ffpc_terms <- lapply(
      terms[where_specials$ffpc],
      \(x) eval(x, envir = eval_env, enclos = frml_env)
    )
    lapply(ffpc_terms, \(trm) {
      lapply(
        colnames(trm$data),
        \(nm) assign(nm, trm$data[obs_indices, nm], envir = formula_env)
      )
    })
    getFfpcFormula <- function(trm) {
      frmls <- lapply(colnames(trm$data), \(pc) {
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
      paste(unlist(frmls), collapse = " + ")
    }
    new_term_strings[where_specials$ffpc] <- sapply(ffpc_terms, getFfpcFormula)
    ffpc_terms <- lapply(ffpc_terms, \(x) x[names(x) != "data"])
  } else {
    ffpc_terms <- NULL
  }

  # Process pcre terms
  if (length(where_specials$pcre)) {
    pcre_terms <- lapply(
      terms[where_specials$pcre],
      \(x) eval(x, envir = eval_env, enclos = frml_env)
    )
    lapply(pcre_terms, \(trm) {
      if (!is_sparse && all(trm$yind == yind)) {
        lapply(
          colnames(trm$efunctions),
          \(nm)
            assign(
              nm,
              trm$efunctions[rep(1:nyindex, times = nobs), nm],
              envir = formula_env
            )
        )
      } else {
        stopifnot(min(trm$yind) <= min(yind), max(trm$yind) >= max(yind))
        lapply(colnames(trm$efunctions), \(nm) {
          tmp <- approx(
            x = trm$yind,
            y = trm$efunctions[, nm],
            xout = yind_vec,
            method = "linear"
          )$y
          assign(nm, tmp, envir = formula_env)
        })
      }
      assign(trm$idname, trm$id[obs_indices], envir = formula_env)
    })
    new_term_strings[where_specials$pcre] <- sapply(
      pcre_terms,
      \(x) safeDeparse(x$call)
    )
  } else {
    pcre_terms <- NULL
  }

  # Transform smooth terms
  if (length(c(where_specials$s, where_specials$te, where_specials$t2))) {
    new_term_strings[c(
      where_specials$s,
      where_specials$te,
      where_specials$t2
    )] <-
      sapply(
        terms[c(where_specials$s, where_specials$te, where_specials$t2)],
        \(x)
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
  if (length(where_specials$par)) {
    new_term_strings[where_specials$par] <- sapply(
      terms[where_specials$par],
      \(x) transform_par_term(x, yindex_vec_name, bs.yindex)
    )
  }

  # Assign expanded variables
  where_specials$notff <- c(
    where_specials$c,
    where_specials$par,
    where_specials$s,
    where_specials$te,
    where_specials$t2
  )
  if (length(where_specials$notff)) {
    eval_env <- if ("data" %in% names(call))
      list2env(eval.parent(call$data)) else frml_env
    lapply(terms[where_specials$notff], \(x) {
      isC <- safeDeparse(x) %in% sapply(terms[where_specials$c], safeDeparse)
      if (isC)
        x <- formula(paste(
          "~",
          gsub("\\)$", "", gsub("^c\\(", "", deparse(x)))
        ))[[2]]
      nms <- if (!is.null(names(x))) all.vars(x[names(x) %in% c("", "by")]) else
        all.vars(x)
      sapply(nms, \(nm) {
        var <- get(nm, envir = eval_env)
        if (is.matrix(var)) {
          stopifnot(!is_sparse || ncol(var) == nyindex)
          assign(nm, as.vector(t(var)), envir = formula_env)
        } else {
          stopifnot(length(var) == nobs)
          assign(nm, var[obs_indices], envir = formula_env)
        }
      })
    })
  }

  # Step 8: Build mgcv formula and data
  new_formula <- build_mgcv_formula(
    response_name = response_name,
    intercept_string = if (add_f_int) int_string else NULL,
    term_strings = new_term_strings,
    has_intercept = add_f_int,
    formula_env = formula_env
  )
  pffr_data <- build_mgcv_data(formula_env)

  # Step 9: Build call (GLS-specific)
  new_call <- expand.call(pffr, call)
  new_call$yind <- new_call$tensortype <- new_call$bs.int <-
    new_call$bs.yindex <- new_call$algorithm <- new_call$ydata <-
      new_call$hatSigma <- new_call$cond.cutoff <- NULL
  new_call$formula <- new_formula
  new_call$data <- quote(pffr_data)
  new_call[[1]] <- algorithm

  dotargs <- names(new_call)[names(new_call) %in% names(dots)]
  new_call[dotargs] <- dots[dotargs]

  if ("weights" %in% dotargs) {
    if (length(dots$weights) == nobs) {
      new_call$weights <- dots$weights[obs_indices]
    } else if (
      !is.null(dim(dots$weights)) && all(dim(dots$weights) == c(nobs, nyindex))
    ) {
      new_call$weights <- as.vector(t(dots$weights))
    } else {
      stop(
        "weights must be vector of length nobs or matrix of dim (nobs, nyindex)"
      )
    }
  }

  if ("offset" %in% dotargs) {
    if (length(dots$offset) == nobs) {
      new_call$offset <- dots$offset[obs_indices]
    } else if (
      !is.null(dim(dots$offset)) && all(dim(dots$offset) == c(nobs, nyindex))
    ) {
      new_call$offset <- as.vector(t(dots$offset))
    } else {
      stop(
        "offset must be vector of length nobs or matrix of dim (nobs, nyindex)"
      )
    }
  }

  # Step 10: GLS fitting
  new_call$fit <- FALSE
  G <- eval(new_call)
  G <- decorrelate_gam_matrices(G, sqrt_sigma_inv, nobs, nyindex)
  m <- gam(G = G, fit = TRUE, method = method)

  # Post-fit corrections
  m$y <- original_response
  m$fitted.values <- predict(m)
  m$residuals <- m$y - m$fitted.values

  # Step 11: Post-processing using shared functions
  m_smooth <- if (as.character(algorithm) %in% c("gamm4", "gamm"))
    m$gam$smooth else m$smooth

  label_map_result <- pffr_build_label_map(
    new_term_strings = new_term_strings,
    terms = terms,
    add_f_int = add_f_int,
    int_string = int_string,
    m_smooth = m_smooth,
    where_specials = where_specials,
    ffpc_terms = ffpc_terms,
    formula_env = formula_env,
    yindex_vec_name = yindex_vec_name
  )
  label_map <- label_map_result$label_map
  term_map <- label_map_result$term_map

  names(m_smooth) <- sapply(m_smooth, \(x) x$label)
  if (as.character(algorithm) %in% c("gamm4", "gamm")) {
    m$gam$smooth <- m_smooth
  } else {
    m$smooth <- m_smooth
  }

  short_labels <- create_shortlabels(
    label_map = label_map,
    m_smooth = m_smooth,
    yind_name = yind_name,
    where_specials = where_specials
  )

  # Build metadata (with GLS-specific fields)
  ret <- pffr_build_metadata(
    call = call,
    formula = formula,
    term_map = term_map,
    label_map = label_map,
    short_labels = short_labels,
    response_name = response_name,
    nobs = nobs,
    nyindex = nyindex,
    yind_name = yind_name,
    yind = yind,
    where_specials = where_specials,
    ff_terms = ff_terms,
    ffpc_terms = ffpc_terms,
    pcre_terms = pcre_terms,
    missing_indices = missing_indices,
    is_sparse = is_sparse,
    ydata = ydata,
    sandwich = FALSE
  )
  # Add GLS-specific fields
  ret$hatSigma <- hatSigma
  ret$sqrtSigmaInv <- sqrt_sigma_inv

  m <- pffr_attach_metadata(m, algorithm, ret)
  m
}


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
