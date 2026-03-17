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
  as.vector(vapply(
    conf,
    \(c_level) {
      probs <- c((1 - c_level) / 2, 1 - (1 - c_level) / 2)
      stats::quantile(
        x,
        probs = probs,
        type = 8,
        na.rm = TRUE,
        names = FALSE
      )
    },
    numeric(2)
  ))
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

  lower_names <- paste0(
    "lower_",
    vapply(conf, format_boot_conf_label, character(1))
  )
  upper_names <- paste0(
    "upper_",
    vapply(conf, format_boot_conf_label, character(1))
  )
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
  if (
    !is.numeric(conf) || any(!is.finite(conf)) || any(conf <= 0 | conf >= 1)
  ) {
    stop("'conf' must contain values strictly between 0 and 1.")
  }
  conf <- as.numeric(conf)
  conf <- conf[!duplicated(conf)]
  primary_conf <- conf[1L]

  if (is.null(ncpus)) {
    ncpus <- getOption("boot.ncpus", 1L)
  }
  if (
    !is.numeric(ncpus) || length(ncpus) != 1 || !is.finite(ncpus) || ncpus < 1
  ) {
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
#' @return A list with similar structure as the return value of
#'   \code{\link{coef.pffr}}, containing the original point estimates along
#'   with bootstrap CIs.
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
      n_failed,
      " out of ",
      B,
      " bootstrap replicates failed (",
      round(100 * failure_rate, 1),
      "%). ",
      "Using ",
      B - n_failed,
      " successful replicates for CI computation. ",
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
#' (deprecated)
#'
#' @description
#' **Deprecated**
#'
#' `pffr_gls()` is deprecated. Its GLS-based covariance correction produced
#' poorly calibrated inference in practice. Use \code{\link{pffr}()} with
#' \code{sandwich = "cluster"} (the new default) or \code{sandwich = "cl2"}
#' for robust covariance estimation instead.
#'
#' @param formula,yind,hatSigma,data,ydata,algorithm,method,tensortype,bs.yindex,bs.int,cond.cutoff,...
#'   Ignored. The function always errors.
#'
#' @return Does not return; always errors with a deprecation message.
#' @seealso \code{\link{pffr}} for the recommended approach with sandwich CIs.
#'
#' @export
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
  .Deprecated(
    new = "pffr",
    msg = paste(
      "pffr_gls() is deprecated and now errors.",
      "Its GLS-based covariance correction produced poorly calibrated inference.",
      "Use pffr() with sandwich = \"cluster\" (the default) or",
      "sandwich = \"cl2\" for robust covariance estimation instead.",
      "See ?pffr for details."
    )
  )
  stop(
    "pffr_gls() has been removed. ",
    "Use pffr() with sandwich = \"cluster\" or sandwich = \"cl2\" instead.",
    call. = FALSE
  )
}


#' Penalized function-on-function regression with non-i.i.d. residuals
#' (deprecated)
#'
#' @description
#' **Deprecated**
#'
#' `pffrGLS()` is deprecated. Use \code{\link{pffr}()} with
#' \code{sandwich = "cluster"} or \code{sandwich = "cl2"} instead.
#'
#' @inheritParams pffr_gls
#' @return Does not return; always errors with a deprecation message.
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
  pffr_gls()
}
