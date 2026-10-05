#' Core Modular Functions for pffr
#'
#' These functions implement the core pipeline for pffr() model fitting,
#' extracted for better maintainability and testing. See PLAN.md Phase 7
#' for architectural design rationale.
#'
#' @name pffr-core
#' @keywords internal
NULL


#--------------------------------------
# Step 1: Argument Validation
#--------------------------------------

#' Validate pffr arguments and process dots
#'
#' Validates the arguments passed to pffr() and processes the `...` arguments.
#' Returns structured information about AR settings, family, and validated dots.
#'
#' @param call The matched call from pffr().
#' @param algorithm The algorithm argument (may be NA).
#' @param dots The list of additional arguments (...).
#' @param check_ar Whether to check AR-related arguments.
#' @returns A list with:
#'   - `dots`: Validated and possibly modified dots
#'   - `use_ar`: Logical, whether AR(1) errors are requested
#'   - `rho`: The rho value (or NULL)
#'   - `gaulss`: Logical, whether using gaulss family
#' @keywords internal
pffr_validate_dots <- function(call, algorithm, dots, check_ar = TRUE) {
  rho_arg <- dots[["rho"]]
  use_ar <- FALSE
  gaulss <- FALSE

  # Check for unsupported AR.start
  if (check_ar && "AR.start" %in% names(dots)) {
    stop(
      "Please do not supply `AR.start` directly; pffr constructs it automatically when `rho` is specified."
    )
  }

  # Determine valid arguments based on algorithm
  valid_dots <- if (!is.na(algorithm) && algorithm == "gamm4") {
    c(names(formals(gamm4::gamm4)), names(formals(lme4::lmer)))
  } else {
    c(
      names(formals(mgcv::gam)),
      names(formals(mgcv::bam)),
      names(formals(mgcv::gam.fit3)),
      names(formals(mgcv::jagam))
    )
  }

  # Check for gaulss family
  if (!is.null(dots$family)) {
    if (
      (is.character(dots$family) && dots$family == "gaulss") ||
        (is.list(dots$family) && dots$family$family == "gaulss")
    ) {
      valid_dots <- c(valid_dots, "varformula")
      gaulss <- TRUE
    }
  }

  # Warn about unused arguments
  not_used <- names(dots)[!(names(dots) %in% valid_dots)]
  if (length(not_used)) {
    warning(
      "Arguments <",
      paste(not_used, collapse = ", "),
      "> supplied but not used."
    )
  }

  # Validate rho argument
  if (check_ar && !is.null(rho_arg)) {
    if (!(is.numeric(rho_arg) && length(rho_arg) == 1 && !is.na(rho_arg))) {
      stop("`rho` must be a single numeric value.")
    }
    if (abs(rho_arg) >= 1) {
      stop("`rho` must have absolute value strictly less than 1.")
    }
    use_ar <- abs(rho_arg) > 0
  }

  # Validate family for AR(1) errors
  if (use_ar) {
    family_obj <- gaussian()
    if ("family" %in% names(dots) && !is.null(dots$family)) {
      fam <- dots$family
      if (is.character(fam)) {
        if (length(fam) != 1) {
          stop("Character `family` specifications must have length 1.")
        }
        fam <- match.fun(fam)
      }
      if (is.function(fam)) {
        fam <- fam()
      }
      if (!is.list(fam) || is.null(fam$family) || is.null(fam$link)) {
        stop("Unable to interpret `family` argument when `rho` is supplied.")
      }
      family_obj <- fam
    }

    gaussian_identity <- identical(family_obj$family, "gaussian") &&
      identical(family_obj$link, "identity")
    discrete_specified <- "discrete" %in% names(dots)
    discrete_requested <- isTRUE(dots$discrete)

    if (!gaussian_identity) {
      if (discrete_specified && !discrete_requested) {
        stop(
          "Autocorrelated errors (via `rho`) require either a Gaussian identity model or `discrete = TRUE` (see ?mgcv::bam)."
        )
      }
      if (!discrete_specified) {
        dots$discrete <- TRUE
      }
    }
  }

  list(
    dots = dots,
    use_ar = use_ar,
    rho = rho_arg,
    gaulss = gaulss
  )
}


#' Validate sparse data format
#'
#' @param ydata The ydata argument.
#' @returns TRUE if valid, stops with error otherwise.
#' @keywords internal
pffr_validate_ydata <- function(ydata) {
  if (!is.null(ydata)) {
    stopifnot(ncol(ydata) == 3)
    stopifnot(c(".obs", ".index", ".value") == colnames(ydata))
  }
  invisible(TRUE)
}


#--------------------------------------
# Step 2: Data Structure Detection
#--------------------------------------

#' Detect data dimensions for pffr
#'
#' @param ydata The ydata argument (NULL for dense data).
#' @param data The data argument.
#' @param response_name The response variable name (symbol).
#' @param frml_env The formula environment.
#' @param eval_env Evaluation environment/list used to resolve variables.
#' @returns A list with nobs, nyindex, ntotal, is_sparse.
#' @keywords internal
pffr_get_dimensions <- function(
  ydata,
  data,
  response_name,
  frml_env,
  eval_env
) {
  is_sparse <- !is.null(ydata)

  if (is_sparse) {
    nobs <- length(unique(ydata$.obs))
    stopifnot(all(ydata$.obs %in% rownames(data)))
    stopifnot(all(ydata$.obs %in% 1:nobs))

    nobs_data <- nrow(as.matrix(data[[1]]))
    stopifnot(nobs == nobs_data)
    ntotal <- nrow(ydata)
    nyindex <- NA_integer_ # Set later
  } else {
    nobs <- nrow(eval(response_name, envir = eval_env, enclos = frml_env))
    nyindex <- ncol(eval(response_name, envir = eval_env, enclos = frml_env))
    ntotal <- nobs * nyindex
  }

  list(
    is_sparse = is_sparse,
    nobs = nobs,
    nyindex = nyindex,
    ntotal = ntotal
  )
}


#' Select and configure algorithm
#'
#' @param algorithm User-specified algorithm (may be NA).
#' @param ntotal Total number of data points.
#' @param call The matched call (modified in place for method).
#' @param where_specials List of special term indices.
#' @param use_ar Whether AR(1) errors requested.
#' @returns Symbol for the algorithm.
#' @keywords internal
pffr_configure_algorithm <- function(
  algorithm,
  ntotal,
  call,
  where_specials,
  use_ar
) {
  algorithm_specified <- !is.na(algorithm)

  if (!algorithm_specified) {
    # Default: bam for large data or when AR(1) requested, gam otherwise
    algorithm <- if (use_ar || ntotal > 1e5) "bam" else "gam"
  }

  algorithm <- as.symbol(algorithm)

  # No te-terms possible in gamm4
  if (as.character(algorithm) == "gamm4") {
    stopifnot(length(unlist(where_specials[c("te", "ti")])) < 1)
  }

  # AR(1) errors only supported for bam - error only if explicitly specified non-bam
  if (use_ar && algorithm_specified && as.character(algorithm) != "bam") {
    stop(
      "Autocorrelated errors via `rho` are currently supported only when `algorithm = \"bam\"`."
    )
  }

  algorithm
}


#--------------------------------------
# Step 3: Y-index Handling
#--------------------------------------

#' Setup y-index for sparse data
#'
#' @param ydata The ydata data.frame.
#' @returns A list with yind, yind_name, nyindex.
#' @keywords internal
pffr_setup_yind_sparse <- function(ydata) {
  yind_name <- "yindex"
  yind <- if (length(unique(ydata$.index)) > 100) {
    seq(min(ydata$.index), max(ydata$.index), length.out = 100)
  } else {
    sort(unique(ydata$.index))
  }
  nyindex <- length(yind)

  list(yind = yind, yind_name = yind_name, nyindex = nyindex)
}


#' Setup y-index for dense data
#'
#' @param yind The yind argument (may be missing).
#' @param yind_missing Logical, whether yind was missing in call.
#' @param nyindex Number of y-index points.
#' @param where_specials List of special term indices.
#' @param terms List of parsed terms.
#' @param eval_env Evaluation environment.
#' @param frml_env Formula environment.
#' @param data The data argument.
#' @param yind_expr The yind expression from the original `pffr()` call
#'   (used only for naming/lookup; optional).
#' @returns A list with yind, yind_name.
#' @keywords internal
pffr_setup_yind_dense <- function(
  yind,
  yind_missing,
  nyindex,
  where_specials,
  terms,
  eval_env,
  frml_env,
  data,
  yind_expr = NULL
) {
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
    if (is.null(yind_expr)) yind_expr <- quote(yind)
    if (is.symbol(yind_expr) || is.character(yind)) {
      yind_name <- deparse(yind_expr)
      if (!is.null(data) && !is.null(data[[yind_name]])) {
        yind <- data[[yind_name]]
      } else if (
        is.character(yind) && length(yind) == 1 && !is.null(data[[yind]])
      ) {
        yind_name <- yind
        yind <- data[[yind_name]]
      } else {
        yind_name <- "yindex"
      }
    } else {
      yind_name <- "yindex"
    }
    stopifnot(is.vector(yind), is.numeric(yind), length(yind) == nyindex)
  }

  if (length(yind_name) > 1) yind_name <- "yindex"
  stopifnot(all.equal(order(yind), 1:nyindex))

  list(yind = yind, yind_name = yind_name)
}


#--------------------------------------
# Step 4: Response Data Setup
#--------------------------------------

#' Setup response and indices for formula environment
#'
#' @param is_sparse Whether data is sparse.
#' @param ydata The ydata (or NULL).
#' @param yind The y-index vector.
#' @param yind_name Name of y-index.
#' @param nobs Number of observations.
#' @param nyindex Number of y-index points.
#' @param response_name Response variable name (symbol).
#' @param eval_env Evaluation environment.
#' @param frml_env Formula environment.
#' @param formula_env Environment to assign to.
#' @returns A list with yind_vec, yindex_vec_name, obs_indices, missing_indices.
#' @keywords internal
pffr_setup_response <- function(
  is_sparse,
  ydata,
  yind,
  yind_name,
  nobs,
  nyindex,
  response_name,
  eval_env,
  frml_env,
  formula_env
) {
  if (is_sparse) {
    yind_vec <- ydata$.index
    yindex_vec_name <- as.symbol(paste0(yind_name, ".vec"))
    assign(deparse(yindex_vec_name), yind_vec, envir = formula_env)

    assign(deparse(response_name), ydata$.value, envir = formula_env)

    missing_indices <- NULL
    obs_indices <- ydata$.obs
  } else {
    yind_vec <- rep(yind, times = nobs)
    yindex_vec_name <- as.symbol(paste0(yind_name, ".vec"))
    assign(deparse(yindex_vec_name), yind_vec, envir = formula_env)

    response_vec <- as.vector(t(eval(
      response_name,
      envir = eval_env,
      enclos = frml_env
    )))
    assign(deparse(response_name), response_vec, envir = formula_env)

    missing_indices <- if (anyNA(response_vec)) which(is.na(response_vec)) else
      NULL
    obs_indices <- rep(1:nobs, each = nyindex)
  }

  list(
    yind_vec = yind_vec,
    yindex_vec_name = yindex_vec_name,
    obs_indices = obs_indices,
    missing_indices = missing_indices
  )
}


#--------------------------------------
# Step 5: AR.start Setup
#--------------------------------------

#' Build AR.start indicator for AR(1) errors
#'
#' @param response_name Response variable name.
#' @param formula_env Formula environment.
#' @param obs_indices Observation indices.
#' @keywords internal
pffr_build_ar_start <- function(response_name, formula_env, obs_indices) {
  resp_long <- get(as.character(response_name), envir = formula_env)
  valid_idx <- which(!is.na(resp_long))

  if (!length(valid_idx)) {
    stop("Cannot build AR.start because all responses are missing.")
  }

  start_idx <- rep(FALSE, length(resp_long))
  splits <- split(valid_idx, obs_indices[valid_idx])
  start_positions <- as.integer(vapply(splits, \(idx) idx[1], integer(1)))
  start_idx[start_positions] <- TRUE

  assign("AR.start", start_idx, envir = formula_env)
  invisible(TRUE)
}


#--------------------------------------
# Step 6: Call Building
#--------------------------------------

#' Build and configure the mgcv call
#'
#' @param call Original pffr call.
#' @param algorithm Algorithm symbol.
#' @param new_formula Transformed formula.
#' @param pffr_data Data frame for mgcv.
#' @param dots Validated dots.
#' @param use_ar Whether AR(1) errors requested.
#' @param nobs Number of observations.
#' @param nyindex Number of y-index points.
#' @param obs_indices Observation indices.
#' @returns The configured call object.
#' @keywords internal
pffr_build_call <- function(
  call,
  algorithm,
  new_formula,
  pffr_data,
  dots,
  use_ar,
  nobs,
  nyindex,
  obs_indices
) {
  newcall <- expand.call(pffr, call)
  newcall$yind <- newcall$tensortype <- newcall$bs.int <-
    newcall$bs.yindex <- newcall$algorithm <- newcall$ydata <- NULL
  newcall$sandwich <- NULL
  newcall$cluster <- NULL
  newcall$ncv_blocks <- NULL
  newcall$formula <- new_formula
  newcall$data <- quote(pffr_data)
  newcall[[1]] <- algorithm

  if (as.character(algorithm) == "gamm4" && "method" %in% names(newcall)) {
    # gamm4 deprecated `method` in favor of REML; map legacy pffr `method`.
    method_value <- tryCatch(
      eval(newcall$method, envir = parent.frame()),
      error = function(e) newcall$method
    )
    method_chr <- as.character(method_value)[1]

    if (!("REML" %in% names(newcall)) && method_chr %in% c("REML", "ML")) {
      newcall$REML <- identical(method_chr, "REML")
    }
    if (!(method_chr %in% c("REML", "ML"))) {
      warning(
        "For algorithm = \"gamm4\", only method = \"REML\" or \"ML\" are interpreted. ",
        "Use REML = TRUE/FALSE for direct control."
      )
    }
    newcall$method <- NULL
  }

  if (use_ar) {
    newcall$AR.start <- quote(pffr_data$AR.start)
    if (isTRUE(dots$discrete)) {
      newcall$discrete <- TRUE
    }
  }

  # Transfer dot args
  dotargs <- names(newcall)[names(newcall) %in% names(dots)]
  newcall[dotargs] <- dots[dotargs]

  # Validate subset
  if ("subset" %in% dotargs) {
    stop("<subset>-argument is not supported.")
  }

  # Handle weights
  if ("weights" %in% dotargs) {
    wtsdone <- FALSE
    if (length(dots$weights) == nobs) {
      newcall$weights <- dots$weights[obs_indices]
      wtsdone <- TRUE
    }
    if (
      !is.null(dim(dots$weights)) && all(dim(dots$weights) == c(nobs, nyindex))
    ) {
      newcall$weights <- as.vector(t(dots$weights))
      wtsdone <- TRUE
    }
    if (!wtsdone) {
      stop(
        "weights have to be supplied as a vector with length=rows(data) or a matrix with the same dimensions as the response."
      )
    }
  }

  # Handle offset
  if ("offset" %in% dotargs) {
    ofstdone <- FALSE
    if (length(dots$offset) == nobs) {
      newcall$offset <- dots$offset[obs_indices]
      ofstdone <- TRUE
    }
    if (
      !is.null(dim(dots$offset)) && all(dim(dots$offset) == c(nobs, nyindex))
    ) {
      newcall$offset <- as.vector(t(dots$offset))
      ofstdone <- TRUE
    }
    if (!ofstdone) {
      stop(
        "offsets have to be supplied as a vector with length=rows(data) or a matrix with the same dimensions as the response."
      )
    }
  }

  # Handle jagam
  if (as.character(algorithm) == "jagam") {
    newcall <- newcall[names(newcall) %in% c("", names(formals(mgcv::jagam)))]
    if (is.null(newcall$file)) {
      newcall$file <- tempfile(
        "pffr2jagam",
        tmpdir = getwd(),
        fileext = ".jags"
      )
    }
  }

  newcall
}


#--------------------------------------
# Step 7: Post-Processing
#--------------------------------------

#' Build the label map from terms to smooth labels
#'
#' @param new_term_strings Transformed term strings.
#' @param terms Original parsed terms.
#' @param add_f_int Whether functional intercept was added.
#' @param int_string Intercept string.
#' @param m_smooth List of smooth objects from fitted model.
#' @param where_specials List of special term indices.
#' @param ffpc_terms Processed ffpc terms.
#' @param formula_env Formula environment.
#' @param yindex_vec_name Y-index vector name.
#' @returns The label_map list.
#' @keywords internal
pffr_build_label_map <- function(
  new_term_strings,
  terms,
  add_f_int,
  int_string,
  m_smooth,
  where_specials,
  ffpc_terms,
  formula_env,
  yindex_vec_name
) {
  term_map <- new_term_strings
  names(term_map) <- names(terms)
  if (add_f_int) term_map <- c(term_map, int_string)

  label_map <- as.list(term_map)
  lbls <- sapply(m_smooth, \(x) x$label)

  if (length(c(where_specials$par, where_specials$ffpc))) {
    # Handle parametric terms
    if (length(where_specials$par)) {
      for (w in where_specials$par) {
        if (is.factor(get(names(label_map)[w], envir = formula_env))) {
          where <- sapply(m_smooth, \(x) x$by) == names(label_map)[w]
          label_map[[w]] <- sapply(m_smooth[where], \(x) x$label)
        } else {
          label_map[[w]] <- paste0(
            "s(",
            yindex_vec_name,
            "):",
            names(label_map)[w]
          )
        }
      }
    }

    # Handle ffpc terms
    if (length(where_specials$ffpc)) {
      ind <- 1
      for (w in where_specials$ffpc) {
        where <- sapply(m_smooth, \(x) x$id) == ffpc_terms[[ind]]$id
        label_map[[w]] <- sapply(m_smooth[where], \(x) x$label)
        ind <- ind + 1
      }
    }

    # Match remaining terms
    other_idx <- setdiff(
      seq_along(label_map),
      c(where_specials$par, where_specials$ffpc)
    )
    if (length(other_idx)) {
      label_map[other_idx] <- lbls[pmatch(
        sapply(label_map[other_idx], get_smooth_label_from_term),
        lbls
      )]
    }
  } else {
    label_map[seq_along(label_map)] <- lbls[pmatch(
      sapply(label_map, get_smooth_label_from_term),
      lbls
    )]
  }

  # Fix any NA labels
  nalbls <- sapply(
    label_map,
    \(x) any(is.null(x)) || any(is.na(x[!is.null(x)]))
  )
  if (any(nalbls)) {
    label_map[nalbls] <- term_map[nalbls]
  }

  list(label_map = label_map, term_map = term_map)
}


#' Get smooth label from a term string
#' @keywords internal
get_smooth_label_from_term <- function(x) {
  parsed <- parse(text = x)[[1]]
  if (length(parsed) != 1) {
    tmp <- eval(parsed)
    tmp$label
  } else {
    x
  }
}


#' Build the pffr metadata list
#'
#' @param call The matched call to `pffr()`.
#' @param formula The model formula.
#' @param term_map Mapping of formula terms to their internal term strings.
#' @param label_map Mapping of formula terms to fitted smooth labels.
#' @param short_labels Short labels for fitted smooths.
#' @param response_name Name of the response variable.
#' @param nobs Number of observations.
#' @param nyindex Number of response-index values per observation.
#' @param yind_name Name of the response-index variable.
#' @param yind Response-index values.
#' @param where_specials Locations of special terms in the formula.
#' @param ff_terms Evaluated function-on-function terms.
#' @param ffpc_terms Evaluated function-on-principal-components terms.
#' @param pcre_terms Evaluated principal-component regression terms.
#' @param missing_indices Indices of missing response values.
#' @param is_sparse Whether the response data use the sparse format.
#' @param ydata Sparse response data, or `NULL` for dense data.
#' @param sandwich Sandwich covariance option used for the fit.
#' @keywords internal
pffr_build_metadata <- function(
  call,
  formula,
  term_map,
  label_map,
  short_labels,
  response_name,
  nobs,
  nyindex,
  yind_name,
  yind,
  where_specials,
  ff_terms,
  ffpc_terms,
  pcre_terms,
  missing_indices,
  is_sparse,
  ydata,
  sandwich
) {
  list(
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
    where = where_specials,
    ff = ff_terms,
    ffpc = ffpc_terms,
    pcre_terms = pcre_terms,
    missing_indices = missing_indices,
    is_sparse = is_sparse,
    ydata = ydata,
    sandwich = sandwich,
    # Covariance storage contract version: format 2 keeps $Vp/$Vc/$Ve
    # model-based ALWAYS; the robust covariance lives in $pffr$Vsandwich.
    cov_format = PFFR_COV_STORAGE_FORMAT,
    # Cache for on-demand sandwich recomputation via pffr_vcov(); an
    # environment so results persist on the fit across accessor calls.
    Vsandwich_cache = new.env(parent = emptyenv())
  )
}


#' Attach pffr metadata to model
#'
#' @param m Fitted model.
#' @param algorithm Algorithm symbol.
#' @param ret pffr metadata list.
#' @returns Model with pffr class and metadata.
#' @keywords internal
pffr_attach_metadata <- function(m, algorithm, ret) {
  if (as.character(algorithm) %in% c("gamm4", "gamm")) {
    m$gam$pffr <- ret
    class(m$gam) <- c("pffr", class(m$gam))
  } else {
    m$pffr <- ret
    class(m) <- c("pffr", class(m))
  }
  m
}


#--------------------------------------
# Step 8: Formula Term Processing
#--------------------------------------

#' Process ff/sff terms and assign to formula environment
#'
#' @param ff_terms List of evaluated ff/sff terms.
#' @param yind_vec Y-index vector (stacked).
#' @param obs_indices Observation indices.
#' @param formula_env Formula environment to assign to.
#' @keywords internal
pffr_process_ff_terms <- function(
  ff_terms,
  yind_vec,
  obs_indices,
  formula_env
) {
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
        if (is.null(x$LX)) {
          XStacked <- shift_and_shorten_matrix(XStacked, eff_windows)
        }
      }
    }
    assign(x$yindname, tmat, envir = formula_env)
    assign(x$xindname, smat, envir = formula_env)
    assign(x$LXname, LStacked, envir = formula_env)
    if (is.null(x[["LX"]])) {
      assign(x$xname, XStacked, envir = formula_env)
    }
    invisible(NULL)
  }
  lapply(ff_terms, makeff)
  invisible(NULL)
}


#' Process ffpc terms and return formula strings
#'
#' @param ffpc_terms List of evaluated ffpc terms.
#' @param obs_indices Observation indices.
#' @param yindex_vec_name Y-index vector name (symbol).
#' @param formula_env Formula environment to assign to.
#' @returns Character vector of formula strings for ffpc terms.
#' @keywords internal
pffr_process_ffpc_terms <- function(
  ffpc_terms,
  obs_indices,
  yindex_vec_name,
  formula_env
) {
  # Assign ffpc data to formula environment

  lapply(ffpc_terms, \(trm) {
    lapply(colnames(trm$data), \(nm) {
      assign(nm, trm$data[obs_indices, nm], envir = formula_env)
      invisible(NULL)
    })
    invisible(NULL)
  })

  # Build formula strings
  get_ffpc_formula <- function(trm) {
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

  sapply(ffpc_terms, get_ffpc_formula)
}


#' Process pcre terms and return formula strings
#'
#' @param pcre_terms List of evaluated pcre terms.
#' @param is_sparse Whether data is sparse.
#' @param yind Y-index vector.
#' @param yind_vec Y-index vector (stacked).
#' @param nyindex Number of y-index points.
#' @param nobs Number of observations.
#' @param obs_indices Observation indices.
#' @param formula_env Formula environment to assign to.
#' @returns Character vector of formula strings for pcre terms.
#' @keywords internal
pffr_process_pcre_terms <- function(
  pcre_terms,
  is_sparse,
  yind,
  yind_vec,
  nyindex,
  nobs,
  obs_indices,
  formula_env
) {
  lapply(pcre_terms, \(trm) {
    if (!is_sparse && all(trm$yind == yind)) {
      lapply(colnames(trm$efunctions), \(nm) {
        assign(
          nm,
          trm$efunctions[rep(1:nyindex, times = nobs), nm],
          envir = formula_env
        )
        invisible(NULL)
      })
    } else {
      if (min(trm$yind) > min(yind) || max(trm$yind) < max(yind)) {
        stop("pcre term yind must span at least the range of the response yind")
      }
      lapply(colnames(trm$efunctions), \(nm) {
        tmp <- approx(
          x = trm$yind,
          y = trm$efunctions[, nm],
          xout = yind_vec,
          method = "linear"
        )$y
        assign(nm, tmp, envir = formula_env)
        invisible(NULL)
      })
    }
    assign(trm$idname, trm$id[obs_indices], envir = formula_env)
    invisible(NULL)
  })

  sapply(pcre_terms, \(x) safeDeparse(x$call))
}


#' Expand variables to formula environment for smooth/parametric terms
#'
#' @param terms List of parsed terms.
#' @param where_specials List of term indices by type.
#' @param is_sparse Whether data is sparse.
#' @param nyindex Number of y-index points.
#' @param nobs Number of observations.
#' @param obs_indices Observation indices.
#' @param eval_env Evaluation environment.
#' @param formula_env Formula environment to assign to.
#' @keywords internal
pffr_expand_variables <- function(
  terms,
  where_specials,
  is_sparse,
  nyindex,
  nobs,
  obs_indices,
  eval_env,
  formula_env
) {
  notff_indices <- c(
    where_specials$c,
    where_specials$par,
    where_specials$s,
    where_specials$te,
    where_specials$t2
  )

  if (!length(notff_indices)) return(invisible(NULL))

  lapply(terms[notff_indices], \(x) {
    isC <- safeDeparse(x) %in% sapply(terms[where_specials$c], safeDeparse)
    if (isC) {
      x <- formula(paste(
        "~",
        gsub("\\)$", "", gsub("^c\\(", "", deparse(x)))
      ))[[2]]
    }
    nms <- if (!is.null(names(x))) {
      all.vars(x[names(x) %in% c("", "by")])
    } else {
      all.vars(x)
    }
    sapply(nms, \(nm) {
      var <- get(nm, envir = eval_env)
      if (is.matrix(var)) {
        if (is_sparse && ncol(var) != nyindex) {
          stop("Matrix covariate '", nm, "' must have ", nyindex, " columns")
        }
        assign(nm, as.vector(t(var)), envir = formula_env)
      } else {
        if (length(var) != nobs) {
          stop("Covariate '", nm, "' must have length ", nobs)
        }
        assign(nm, var[obs_indices], envir = formula_env)
      }
      invisible(NULL)
    })
    invisible(NULL)
  })
  invisible(NULL)
}


#--------------------------------------
# Step 9: Sandwich Correction
#--------------------------------------

#' Build cluster ID vector from pffr metadata
#'
#' Maps each row of the vectorized model matrix back to its curve.
#'
#' @param pffr_meta The `pffr` metadata list from a fitted model.
#' @param cluster Optional user-supplied grouping with one entry per curve,
#'   mapping each curve to its independent unit (e.g. a subject id for repeated
#'   measures). Expanded to one entry per vectorized observation so the
#'   cluster-robust sandwich clusters at that level instead of by curve.
#'   `NULL` inherits the fit-time grouping, otherwise clusters by curve.
#'   Custom grouping is supported only for dense responses.
#' @returns Integer vector of length equal to the number of fitted rows.
#' @keywords internal
build_cluster_id <- function(pffr_meta, cluster = NULL) {
  cluster <- cluster %||% pffr_meta$cluster
  if (
    !is.null(cluster) &&
      (!is.atomic(cluster) || !is.null(dim(cluster)) || anyNA(cluster))
  )
    stop(
      "cluster must be an atomic vector without missing values.",
      call. = FALSE
    )
  if (!is.null(cluster)) {
    # User-supplied grouping: one entry per curve (functional observation),
    # mapping each curve to its independent unit (e.g. subject for repeated
    # measures). Expanded to one entry per vectorized observation. This lets
    # the cluster-robust sandwich cluster at the correct level when curves are
    # nested in higher-level units, rather than the default by-curve clustering.
    if (isTRUE(pffr_meta$is_sparse)) {
      stop(
        "Custom `cluster` is only supported for densely-observed (gridded) ",
        "responses, not sparse/irregular fits.",
        call. = FALSE
      )
    }
    if (length(cluster) != pffr_meta$nobs) {
      stop(
        sprintf(
          "`cluster` must have one entry per curve (length %d); got %d.",
          pffr_meta$nobs,
          length(cluster)
        ),
        call. = FALSE
      )
    }
    cluster_id <- rep(cluster, each = pffr_meta$nyindex)
  } else if (isTRUE(pffr_meta$is_sparse)) {
    # The model frame omits ydata rows with a missing .value.
    cluster_id <- pffr_meta$ydata$.obs[!is.na(pffr_meta$ydata$.value)]
  } else {
    cluster_id <- rep(seq_len(pffr_meta$nobs), each = pffr_meta$nyindex)
  }
  if (length(pffr_meta$missing_indices) > 0L) {
    cluster_id <- cluster_id[-pffr_meta$missing_indices]
  }
  cluster_id
}

#' Classify a family's sandwich-score path
#'
#' Determines how the CL2 sandwich builds per-observation scores for a given
#' family. The returned label drives the dispatch in
#' [gam_sandwich_cluster_cl2()] and [pffr_influence()].
#'
#' \describe{
#'   \item{`"gaulss"`}{Gaussian location-scale; exact two-block Fisher-whitened
#'     score ([build_cl2_working_gaulss()]).}
#'   \item{`"scat"`}{Scaled-t (`mgcv::scat`); exact single-block score
#'     ([build_cl2_working_scat()]).}
#'   \item{`"exact"`}{An ordinary exponential-dispersion family (`gaussian`,
#'     `poisson`, `binomial`, `Gamma`, `inverse.gaussian`, quasi-families).
#'     Their generic working-residual score
#'     \eqn{(y-\mu)\,(\mathrm{d}\mu/\mathrm{d}\eta)/(\phi V(\mu))} is the exact
#'     log-likelihood score, so no approximation is involved. For the
#'     quasi-families the estimated \eqn{\hat\phi} (`fit$sig2`) enters this
#'     score once and the bread \eqn{V_p = \hat\phi (X'WX + S)^{-1}} once, so it
#'     cancels from the sampling core \eqn{V_p M V_p} (which is therefore
#'     identical to the fixed-dispersion fit's at the same \eqn{\lambda}) while
#'     the additive \eqn{B_2 = V_p - V_e} allowance scales with \eqn{\hat\phi},
#'     consistently with the \eqn{S/\hat\phi} penalty convention; see
#'     `tests/testthat/test-pffr-quasi-score.R`.}
#'   \item{`"approx"`}{An extended family (`nb`, `tw`, `betar`, `ocat`, ...)
#'     that is neither scaled-t nor location-scale. The generic
#'     working-residual score is only an exponential-family approximation to the
#'     true score; the sandwich builders emit [pffr_warn_approx_score()].}
#'   \item{`"custom"`}{A family defining its own `family$sandwich` (other than
#'     gaulss, e.g. `multinom`); no cluster score factorization is
#'     implemented, so [pffr()] falls back to model-based intervals.}
#' }
#'
#' @param family A family object (ordinary `stats::family` or `mgcv` extended /
#'   general family).
#' @returns One of `"gaulss"`, `"scat"`, `"exact"`, `"approx"`, `"custom"`.
#' @keywords internal
pffr_score_kind <- function(family) {
  fam <- tolower(as.character(family$family))
  if (fam == "gaulss") {
    return("gaulss")
  }
  # scat's family string is "scaled t" (unfitted) or "Scaled t(nu,sig)" (fitted)
  if (grepl("^scaled t", fam)) {
    return("scat")
  }
  if (!is.null(family$sandwich)) {
    return("custom")
  }
  if (inherits(family, "extended.family")) {
    return("approx")
  }
  "exact"
}

#' Warn once that a family's sandwich uses the working-residual approximation
#'
#' For families classified `"approx"` by [pffr_score_kind()] the cluster
#' sandwich builds scores from the exponential-family working residual, which is
#' only an approximation to the true log-likelihood score. This helper emits the
#' review-mandated disclosure, at most once per session per family (keyed via
#' [pffr_warn_once()], so repeated `coef()`/`predict()` recomputations do not
#' spam). Fitted extended families embed estimated parameters in the family
#' string (e.g. `"Negative Binomial(2.403)"`, `"Tweedie(p=1.17)"`), so the
#' warn-once key strips the parenthetical (and case), keeping the semantics
#' once-per-family rather than once-per-theta. Tests reset by removing the
#' `approx_score_*` keys from `refund:::.pffr_state`.
#'
#' @param family A family object.
#' @returns Invisibly `NULL`; called for the warning side effect.
#' @keywords internal
pffr_warn_approx_score <- function(family) {
  fam <- as.character(family$family)
  key_fam <- tolower(trimws(sub("\\(.*$", "", fam)))
  pffr_warn_once(
    paste0("approx_score_", key_fam),
    paste0(
      "sandwich scores for family '",
      fam,
      "' use the exponential-family working-residual approximation"
    )
  )
}

#' Build CL2 working representation for standard single-LP families
#'
#' Factorizes per-observation scores as `Xw_i * z_i` to avoid constructing
#' a dense score matrix.
#'
#' @param b Fitted GAM object.
#' @param cluster_id Cluster vector.
#' @returns List with `Xw`, `z`, and `cluster_id`.
#' @keywords internal
build_cl2_working_standard <- function(b, cluster_id) {
  X <- model.matrix(b)
  y <- as.vector(b$y)
  mu <- as.vector(b$fitted.values)
  eta <- as.vector(b$linear.predictors)

  var_mu <- as.vector(b$family$variance(mu))
  mu_eta <- as.vector(b$family$mu.eta(eta))
  sig2 <- b$sig2
  if (is.null(sig2) || !is.finite(sig2) || sig2 <= 0) sig2 <- 1

  pw <- b$prior.weights
  if (is.null(pw)) pw <- rep(1, length(y))
  pw <- as.vector(pw)

  denom <- sig2 * var_mu
  sqrt_common <- sqrt(pw / denom)
  sqrt_common[!is.finite(sqrt_common)] <- 0

  x_scale <- mu_eta * sqrt_common
  x_scale[!is.finite(x_scale)] <- 0

  z <- (y - mu) * sqrt_common
  z[!is.finite(z)] <- 0

  Xw <- X * x_scale
  Xw[!is.finite(Xw)] <- 0

  list(Xw = Xw, z = z, cluster_id = cluster_id)
}

#' Build CL2 working representation for gaulss family
#'
#' Uses a Fisher-weighted two-block IRLS whitening, with one pseudo-row per
#' linear predictor component (location and scale). For `gaulss` the expected
#' Fisher information is block-diagonal in the location (`eta1`) and scale
#' (`eta2`) predictors,
#' \deqn{I_{11} = \tau^2,\quad I_{22} = 2 (\mathrm{d}\tau/\mathrm{d}\eta_2)^2 /
#' \tau^2,\quad I_{12} = 0,}
#' and each score block factorizes as (Fisher weight) times (residual). The two
#' blocks therefore decouple: each pseudo-row is whitened with its own Fisher
#' working weight \eqn{W_k}, putting `Xtilde_k = sqrt(W_k) X_k` into the design
#' and `z_k = w_k / sqrt(W_k)` into the residual, where `w_k` are the
#' (prior-weighted) gaulss score weights. Both blocks
#' reconstruct the exact gaulss score (`Xtilde_k^T z_k = w_k X_k`), but because
#' `Xtilde` carries the Fisher weight, the per-cluster hat block
#' `H_gg = Xtilde_g Vp Xtilde_g^T` is the genuine penalized hat: its total trace
#' equals the model EDF.
#'
#' This replaces an earlier `sign(w) sqrt(|w|)` sign-split factorization, whose
#' hat block scaled with the residual magnitude `|r|` rather than the
#' (dimensionless) leverage and so collapsed the small-cluster BRL leverage
#' inflation.
#'
#' @param b Fitted GAM object with `family = gaulss`.
#' @param cluster_id Cluster vector for original observations.
#' @returns List with `Xw`, `z`, and expanded `cluster_id`.
#' @keywords internal
build_cl2_working_gaulss <- function(b, cluster_id) {
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

  mu_eta1 <- b$family$linfo[[1]]$mu.eta(eta1) # location link derivative
  mu_eta2 <- b$family$linfo[[2]]$mu.eta(eta2) # dtau/deta2

  # Score weights w.r.t. the location and scale linear predictors
  w1 <- (tau^2 * r) * mu_eta1
  w2 <- (1 / tau - tau * r^2) * mu_eta2

  pw <- b$prior.weights
  if (is.null(pw)) pw <- rep(1, n_obs)
  if (any(pw != 1)) {
    w1 <- pw * w1
    w2 <- pw * w2
  }

  # Expected Fisher working weights (block-diagonal location/scale).
  # General (link-robust) location weight uses mu_eta1; W2 uses dtau/deta2.
  W1 <- pw * tau^2 * mu_eta1^2 # location
  W2 <- pw * 2 * mu_eta2^2 / tau^2 # scale

  make_block <- function(W, w, lp_idx) {
    s <- sqrt(W)
    s[!is.finite(s)] <- 0

    # z_k = w_k / sqrt(W_k); guard zero/non-finite Fisher weights.
    z <- w / s
    z[!is.finite(z)] <- 0

    Xw_block <- matrix(0, nrow = n_obs, ncol = ncol(X))
    if (!is.null(lp_idx) && length(lp_idx) > 0) {
      Xw_block[, lp_idx] <- X[, lp_idx, drop = FALSE] * s
      Xw_block[!is.finite(Xw_block)] <- 0
    }
    list(Xw = Xw_block, z = z)
  }

  block1 <- make_block(W1, w1, lpi[[1]])
  block2 <- make_block(W2, w2, lpi[[2]])

  list(
    Xw = rbind(block1$Xw, block2$Xw),
    z = c(block1$z, block2$z),
    cluster_id = rep(cluster_id, times = 2L)
  )
}

#' Build CL2 working representation for the scaled-t (scat) family
#'
#' Factorizes the exact scaled-t location score `s_i = z_i * Xw_i` using the
#' Fisher-whitened design the CL2 leverage adjustment needs. With location
#' Fisher information (in \eqn{\mu}-space)
#' \deqn{w = (\nu+1) / ((\nu+3)\,\sigma^2)} (the standard t-location result,
#' equal to `0.5 * family$Dd()$EDmu2`), the eta-space whitening weight is
#' \eqn{W_i = (\mathrm{d}\mu/\mathrm{d}\eta)_i^2\, w} WITHOUT the prior weights
#' \eqn{\omega_i}: `mgcv::scat`'s `Dd()$EDmu2` deliberately omits `wt`, so the
#' Fisher weights `gam.fit4` stores (`wf = pmax(0, EDeta2/2)` = `fit$weights`)
#' and uses for `Vp`/EDF are prior-weight-free (verified empirically:
#' `fit$weights == mu_eta^2 * w` and `trace(Vp X'diag(fit$weights)X) ==
#' sum(edf)` to 1e-10 on a weighted fit). The prior weights enter the residual
#' instead,
#' \deqn{\tilde x_i = \sqrt{W_i}\, x_i, \qquad z_i = \omega_i\,
#' (\partial\ell/\partial\mu)_i\, (\mathrm{d}\mu/\mathrm{d}\eta)_i / \sqrt{W_i},}
#' so `Xw^T z` still reconstructs the exact (prior-weighted) score, and the
#' per-cluster hat block `H_gg = Xw_g Vp Xw_g^T` is the genuine penalized hat of
#' the fit (its total trace equals the model EDF), mirroring
#' [build_cl2_working_gaulss()].
#'
#' @param b Fitted GAM object with `family = scat()`.
#' @param cluster_id Cluster vector.
#' @returns List with `Xw`, `z`, and `cluster_id`.
#' @keywords internal
build_cl2_working_scat <- function(b, cluster_id) {
  X <- model.matrix(b)
  th <- b$family$getTheta(TRUE)
  nu <- th[1]
  sig <- th[2]

  mu <- as.vector(b$fitted.values)
  eta <- as.vector(b$linear.predictors)
  r <- as.vector(b$y) - mu
  mu_eta <- as.vector(b$family$mu.eta(eta))

  pw <- b$prior.weights
  if (is.null(pw)) pw <- rep(1, length(mu))
  pw <- as.vector(pw)

  dl_dmu <- (nu + 1) * r / (nu * sig^2 + r^2)
  w_fisher <- (nu + 1) / ((nu + 3) * sig^2) # Fisher info in mu-space (scalar)

  # eta-space expected Fisher whitening weight, WITHOUT prior weights: mgcv's
  # scat EDmu2 omits wt, so gam.fit4's stored Fisher weights (fit$weights) and
  # hence Vp/EDF are prior-weight-free; including pw here would break the
  # trace(H) = EDF hat identity on weighted fits (measured +21% at pw~U(0.5,2)).
  # The prior weights enter the residual z below, keeping Xw^T z exact.
  W <- (mu_eta^2) * w_fisher
  s <- sqrt(W)
  s[!is.finite(s)] <- 0

  score_eta <- pw * dl_dmu * mu_eta # d l / d eta, prior-weighted
  z <- score_eta / s
  z[!is.finite(z)] <- 0

  Xw <- X * s
  Xw[!is.finite(Xw)] <- 0

  list(Xw = Xw, z = z, cluster_id = cluster_id)
}

#' Cluster-robust CL2 sandwich covariance estimator
#'
#' Computes the curve-clustered CL2 covariance with the full-block ("exact")
#' Bell--McCaffrey leverage adjustment
#' \eqn{A_g = \{(I - H)^2\}_{gg}^{-1/2}}{A_g = ((I - H)^2)_gg^(-1/2)} in its
#' Bayesian form, i.e. \eqn{G/(G-1)\,V_p \sum_g \tilde X_g^\top A_g z_g
#' z_g^\top A_g \tilde X_g V_p + (V_p - V_e)}, see [pffr_influence_core()].
#'
#' @param b Fitted GAM object (must not have class `"pffr"`) with model-based
#'   covariance matrices.
#' @param cluster_id Vector mapping each vectorized observation to a curve (or
#'   cluster). Ignored if `influence` is supplied.
#' @param leverage_cap Eigenvalue floor setting, see [pffr_influence_core()].
#' @param influence Optional precomputed influence object
#'   ([pffr_influence()]).
#' @returns A p x p covariance matrix with leverage diagnostics as attributes:
#'   `n_adjusted` (blocks floored at `(1 - leverage_cap)^2`), `min_block_eig`,
#'   `max_block_kappa`, `min_block_eig_rel`, `max_leverage`, and the
#'   hat-invariant monitors `max_obs_leverage`, `min_obs_leverage`,
#'   `min_hat_eig` and `hat_invariant_violation` (`NULL` when the penalized hat
#'   respects its bounds, otherwise a description of the violation; see
#'   [pffr_hat_invariant_violation()]). A violation is reported with one
#'   warning.
#' @keywords internal
gam_sandwich_cluster_cl2 <- function(
  b,
  cluster_id,
  leverage_cap = 0.999,
  influence = NULL
) {
  if (is.null(influence)) {
    kind <- pffr_score_kind(b$family)
    if (kind == "custom")
      stop(
        "No cluster-robust covariance is available for this family.",
        call. = FALSE
      )
    if (kind == "approx") pffr_warn_approx_score(b$family)
    work <- switch(
      kind,
      gaulss = build_cl2_working_gaulss(b, cluster_id),
      scat = build_cl2_working_scat(b, cluster_id),
      build_cl2_working_standard(b, cluster_id)
    )
    influence <- pffr_influence_core(
      work$Xw,
      b$Vp,
      work$cluster_id,
      work$z,
      leverage_cap = leverage_cap
    )
    influence$B2 <- (b$Vp + t(b$Vp)) / 2 - b$Ve
  }
  V <- pffr_influence_vcov(influence)
  # Study-LB P-LB5: a penalized hat that has broken its own bounds means the
  # bread and the weighted design are numerically inconsistent. The (1 - cap)^2
  # floor on the residual-block eigenvalues is a routine numerical safeguard
  # and stays silent.
  hat_violation <- attr(V, "hat_invariant_violation")
  if (!is.null(hat_violation)) {
    warning(
      "Cluster-robust covariance is NOT trustworthy for this fit: ",
      hat_violation,
      " The penalized hat matrix satisfies 0 <= h_ii <= 1 and ",
      "0 <= eigen(H_gg) <= 1 exactly, so this indicates that the model-based ",
      "bread and the weighted design have become numerically inconsistent -- ",
      "typically an ill-conditioned or barely converged fit (extreme fitted ",
      "values, a near-singular penalized Hessian, or a basis far too rich for ",
      "the data). Inspect and refit the model.",
      call. = FALSE
    )
  }
  V
}

#' Numerical-sanity check on the per-cluster hat blocks
#'
#' The penalized hat \eqn{H = X_w V_p X_w'} obeys \eqn{0 \le h_{ii} \le 1} and,
#' since any principal block of a symmetric matrix has its eigenvalues inside
#' the parent's range, \eqn{0 \le \mathrm{eigen}(H_{gg}) \le 1}. The exact
#' Bell--McCaffrey block \eqn{B_g = ((I - H)^2)_{gg}} is positive semi-definite
#' for the same reason. When a fit is numerically degenerate --- an
#' ill-conditioned or barely converged fit, e.g. a Poisson fit whose fitted
#' means span many orders of magnitude --- the bread \eqn{V_p} stops being the
#' inverse of the same weighted cross-product, and these invariants break by
#' orders of magnitude rather than by round-off (study LB, claim P-LB5).
#'
#' @param max_obs_leverage Largest per-observation leverage \eqn{h_{ii}} seen,
#'   or `NA`.
#' @param min_block_eig_rel Smallest eigenvalue of any \eqn{B_g} relative to
#'   that block's largest eigenvalue.
#' @param min_obs_leverage Smallest per-observation leverage \eqn{h_{ii}} seen,
#'   or `NA`. A negative value means the penalized hat is indefinite.
#' @param min_hat_eig Smallest eigenvalue of any \eqn{H_{gg}}, or `NA`.
#' @param tol Relative slack allowed before an invariant counts as violated.
#' @returns `NULL` when every invariant holds, otherwise a one-sentence
#'   character description of the violation.
#' @keywords internal
pffr_hat_invariant_violation <- function(
  max_obs_leverage = NA_real_,
  min_block_eig_rel = 0,
  min_obs_leverage = NA_real_,
  min_hat_eig = NA_real_,
  tol = 1e-6
) {
  msgs <- character(0)
  if (isTRUE(is.finite(max_obs_leverage) && max_obs_leverage > 1 + tol)) {
    msgs <- c(
      msgs,
      sprintf(
        "the largest per-observation leverage is %.3g, above the bound 1;",
        max_obs_leverage
      )
    )
  }
  if (isTRUE(is.finite(min_obs_leverage) && min_obs_leverage < -tol)) {
    msgs <- c(
      msgs,
      sprintf(
        "the smallest per-observation leverage is %.3g, below the bound 0;",
        min_obs_leverage
      )
    )
  }
  if (isTRUE(is.finite(min_hat_eig) && min_hat_eig < -tol)) {
    msgs <- c(
      msgs,
      sprintf(
        paste0(
          "the smallest per-cluster hat eigenvalue is %.3g, but H_gg is ",
          "positive semi-definite by construction;"
        ),
        min_hat_eig
      )
    )
  }
  if (isTRUE(is.finite(min_block_eig_rel) && min_block_eig_rel < -tol)) {
    msgs <- c(
      msgs,
      sprintf(
        paste0(
          "the exact Bell-McCaffrey block has a relative eigenvalue of %.3g, ",
          "but it is positive semi-definite by construction;"
        ),
        min_block_eig_rel
      )
    )
  }
  if (length(msgs) == 0) return(NULL)
  paste(msgs, collapse = " ")
}

#' Resolve a pffr fit's covariance matrices to a canonical representation
#'
#' The current storage contract (format 2) keeps `Vp`/`Vc`/`Ve` model-based
#' ALWAYS and stores the robust covariance in `object$pffr$Vsandwich`. Fits
#' made with intermediate development versions (format 1) instead overwrote
#' `Vp`/`Vc`/`Ve` with the robust matrices and stashed the model-based
#' originals in `object$pffr$model_cov`; fits made with refund 0.1-40 and
#' earlier overwrote them without a stash. This helper resolves each layout
#' into a single shape, emitting a session-scoped one-time back-compat warning
#' for old fits. All covariance recomputation must use the model-based
#' matrices as the (penalized) bread, or the sandwich is applied on top of
#' itself.
#'
#' @param object A fitted pffr model.
#' @returns A list with `model` (list of model-based `Vp`/`Vc`/`Ve`),
#'   `fit_type` (the fit-time sandwich type, `"none"` for model-based fits),
#'   `Vsandwich` (the fit-time robust matrix or `NULL`), and `format`
#'   (`"new"`/`"old"`/`"ancient"`).
#' @keywords internal
pffr_canonicalize_cov <- function(object) {
  meta <- object$pffr
  fit_type <- normalize_sandwich_type(meta$sandwich)

  # Format 2 (current): model-based Vp/Vc/Ve, robust in $pffr$Vsandwich. The
  # cov_format stamp is set at metadata-build time, so this branch is also
  # taken during fitting, before apply_sandwich_correction() has stored
  # Vsandwich (Vp/Vc/Ve are model-based there too).
  if (
    isTRUE((meta$cov_format %||% 0L) >= 2L) ||
      !is.null(meta$sandwich_info) ||
      !is.null(meta[["Vsandwich"]])
  ) {
    return(list(
      model = list(Vp = object$Vp, Vc = object$Vc, Ve = object$Ve),
      fit_type = meta$sandwich_info$type %||% fit_type,
      Vsandwich = meta[["Vsandwich"]],
      format = "new"
    ))
  }

  # Plain fit: no sandwich applied, Vp/Vc/Ve are model-based.
  if (identical(fit_type, "none")) {
    return(list(
      model = list(Vp = object$Vp, Vc = object$Vc, Ve = object$Ve),
      fit_type = "none",
      Vsandwich = NULL,
      format = "new"
    ))
  }

  # Format 1: Vp/Vc/Ve overwritten with robust, model-based in model_cov.
  mc <- meta$model_cov
  if (!is.null(mc)) {
    pffr_warn_once(
      "oldformat_modelcov",
      paste0(
        "This pffr fit was created by an older refund version: its $Vp/$Vc/$Ve ",
        "hold a sandwich-adjusted covariance and the model-based matrices are ",
        "stashed in $pffr$model_cov. Intervals are recomputed from the ",
        "model-based matrices. Refit with the current refund version to ",
        "store the current format."
      )
    )
    return(list(
      model = list(Vp = mc$Vp, Vc = mc$Vc, Ve = mc$Ve),
      fit_type = fit_type,
      Vsandwich = object$Vp,
      format = "old"
    ))
  }

  # refund <= 0.1-40: robust matrices, model-based irrecoverably lost.
  pffr_warn_once(
    "oldformat_nostash",
    paste0(
      "This pffr fit was created with sandwich = \"",
      fit_type,
      "\" by refund 0.1-40 or earlier and does not carry the model-based ",
      "covariance matrices. Its stored covariance is used with Gaussian ",
      "critical values; refit the model with the current refund version for ",
      "CL2 intervals with Satterthwaite critical values."
    )
  )
  list(
    model = list(Vp = object$Vp, Vc = object$Vc, Ve = object$Ve),
    fit_type = fit_type,
    Vsandwich = object$Vp,
    format = "ancient"
  )
}

#' Return a fit's underlying gam with model-based covariance matrices
#'
#' Strips the `"pffr"` class and guarantees `Vp`/`Vc`/`Ve` are the model-based
#' matrices, so downstream mgcv machinery and the sandwich estimator use the
#' penalized bread rather than an already-robustified matrix.
#'
#' @param object A fitted pffr model.
#' @returns The underlying gam object with model-based `Vp`/`Vc`/`Ve`.
#' @keywords internal
pffr_model_based_gam <- function(object) {
  canon <- pffr_canonicalize_cov(object)
  object$Vp <- canon$model$Vp
  object$Vc <- canon$model$Vc
  object$Ve <- canon$model$Ve
  class(object) <- setdiff(class(object), "pffr")
  object
}

#' Does a fit use an AR(1) working correlation?
#'
#' `TRUE` when the fit was estimated by [mgcv::bam()] with a nonzero `rho`
#' (stored by mgcv as `$AR1.rho`), i.e. on AR(1)-whitened residuals.
#'
#' @param object A fitted pffr model or its stripped gam version.
#' @returns A single logical.
#' @keywords internal
pffr_is_ar1_fit <- function(object) {
  rho <- object$AR1.rho
  is.numeric(rho) && length(rho) == 1 && is.finite(rho) && rho != 0
}

#' Refuse the sandwich covariance for AR(1) working-correlation fits
#'
#' The CL2 sandwich builds its scores from the working-independence design and
#' residuals. On a fit with an AR(1) working correlation (`rho` in [pffr()])
#' the estimating equations are the AR-whitened ones, so the sandwich would be
#' silently wrong. This is the single guard shared by every sandwich path
#' ([pffr_influence()], [pffr_vcov()] and the fit-time check in [pffr()]).
#'
#' @param object A fitted pffr model or its stripped gam version, or `NULL`
#'   when `rho` is given directly.
#' @param sandwich Is the sandwich requested?
#' @param rho Optional AR(1) coefficient; defaults to the fit's `$AR1.rho`.
#' @returns `invisible(TRUE)`; errors for a sandwich request on an AR(1) fit.
#' @keywords internal
pffr_check_sandwich_ar1 <- function(object, sandwich, rho = object$AR1.rho) {
  if (!isTRUE(sandwich) || !pffr_is_ar1_fit(list(AR1.rho = rho))) {
    return(invisible(TRUE))
  }
  stop(
    "sandwich = TRUE is not supported for fits with an AR(1) working ",
    "correlation (rho = ",
    format(rho, digits = 3),
    "): the sandwich assumes working-independence scores and would be wrong ",
    "on the AR-whitened fit. Use sandwich = FALSE with the AR(1) working ",
    "model, or drop `rho` to use the CL2 sandwich.",
    call. = FALSE
  )
}

#' Model-based covariance of a fit
#'
#' \eqn{V_c}, which adds mgcv's smoothing-parameter uncertainty correction,
#' falling back to \eqn{V_p}. For NCV fits mgcv's \eqn{V_c} is not a valid
#' covariance (it is much smaller than \eqn{V_p}), so they use \eqn{V_p}.
#'
#' @param object A fitted pffr model.
#' @param model List of model-based `Vp`/`Vc`/`Ve`.
#' @returns A covariance matrix.
#' @keywords internal
pffr_model_based_cov <- function(object, model) {
  if (identical(object$method, "NCV")) return(model$Vp)
  model$Vc %||% model$Vp
}

#' Resolve a user-supplied `sandwich` argument
#'
#' `TRUE` requests the CL2 sandwich, `FALSE` the model-based covariance and
#' `NULL` (only where documented) inherits the fit. The character values of
#' refund 0.1-40 are deprecated: `"none"` maps to `FALSE`, `"cl2"` to `TRUE`,
#' and `"cluster"` and `"hc"`, which are no longer available, to `TRUE` (CL2).
#'
#' @param sandwich The argument as supplied.
#' @returns `TRUE`, `FALSE` or `NULL`.
#' @keywords internal
pffr_sandwich_arg <- function(sandwich) {
  if (is.null(sandwich)) return(NULL)
  if (is.logical(sandwich) && length(sandwich) == 1L && !is.na(sandwich)) {
    return(sandwich)
  }
  legacy <- c("cluster", "cl2", "hc", "none")
  if (
    is.character(sandwich) && length(sandwich) == 1L && sandwich %in% legacy
  ) {
    use <- sandwich != "none"
    .Deprecated(
      msg = paste0(
        "sandwich = \"",
        sandwich,
        "\" is deprecated; use sandwich = ",
        use,
        ".",
        if (sandwich %in% c("cluster", "hc")) {
          paste0(
            " The \"",
            sandwich,
            "\" covariance is no longer available; the CL2 sandwich is ",
            "used instead."
          )
        } else {
          ""
        }
      )
    )
    return(use)
  }
  stop(
    "`sandwich` must be TRUE (CL2 sandwich, the default) or FALSE ",
    "(model-based covariance).",
    call. = FALSE
  )
}

#' Error on arguments removed from the pffr inference interface
#'
#' Arguments of development versions of refund that were removed rather than
#' deprecated (they were never released) would otherwise disappear silently
#' into `...`.
#'
#' @param dots Named list of the `...` arguments.
#' @param fun Name of the calling function, for the message.
#' @returns `invisible(NULL)`; errors if a removed argument is present.
#' @keywords internal
pffr_check_removed_args <- function(dots, fun) {
  removed <- c(
    "freq",
    "cl2_adjustment",
    "dof_correction",
    "edf_type",
    "crit",
    "ci_ref",
    "df_gram",
    "bias_ref",
    "se_method"
  )
  found <- intersect(names(dots), removed)
  if (length(found)) {
    stop(
      fun,
      "(): argument(s) ",
      paste0("`", found, "`", collapse = ", "),
      " no longer exist. Intervals use the CL2 sandwich with Satterthwaite ",
      "critical values (sandwich = TRUE, the default) or the model-based ",
      "covariance with Gaussian critical values (sandwich = FALSE); see ",
      "?pffr.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Resolve the covariance matrix for a pffr fit (single accessor)
#'
#' The one internal accessor through which every refund consumer
#' (`coef.pffr`, `plot.pffr`, `predict.pffr`, `pffr_predict_ci`, simultaneous
#' bands) obtains a covariance matrix, so that `$Vp`/`$Vc`/`$Ve` stay
#' model-based and robust matrices are never double-applied.
#'
#' @param object A fitted pffr model.
#' @param sandwich `NULL` (default) uses the fit-time choice; `TRUE` the CL2
#'   sandwich (exact leverage adjustment, Bayesian form); `FALSE` the
#'   model-based covariance ([pffr_model_based_cov()]).
#' @param cluster Optional per-curve grouping forcing recomputation at a custom
#'   cluster level.
#' @returns A covariance matrix.
#' @keywords internal
pffr_vcov <- function(object, sandwich = NULL, cluster = NULL) {
  canon <- pffr_canonicalize_cov(object)
  use_sandwich <- sandwich %||% !identical(canon$fit_type, "none")
  if (!use_sandwich) {
    return(pffr_model_based_cov(object, canon$model))
  }
  pffr_check_sandwich_ar1(object, TRUE)
  if (identical(canon$format, "ancient")) {
    return(canon$Vsandwich)
  }
  # Serve the fit-time matrix when it is the exact CL2 sandwich.
  stored_exact <- identical(canon$fit_type, "cl2") &&
    identical(object$pffr$sandwich_info$cl2_adjustment, "exact") &&
    !is.null(canon$Vsandwich)
  if (is.null(cluster) && stored_exact) {
    return(canon$Vsandwich)
  }
  # Recompute from the model-based bread (safe by construction) and cache it
  # in fit$pffr$Vsandwich_cache -- an environment, so the cache persists
  # across coef()/predict() calls on the same stored object. Custom `cluster`
  # requests are never cached.
  cache <- object$pffr$Vsandwich_cache
  cacheable <- is.null(cluster) && is.environment(cache)
  if (cacheable && !is.null(cache$cl2)) {
    return(cache$cl2)
  }
  V <- gam_sandwich_cluster_cl2(
    pffr_model_based_gam(object),
    cluster_id = NULL,
    influence = pffr_influence(object, cluster = cluster)
  )
  if (cacheable) cache$cl2 <- V
  V
}

#' Critical-value setup for pointwise intervals
#'
#' Intervals from the CL2 sandwich use Satterthwaite critical values, with a
#' separate working-model moment df for every evaluation point (see
#' [pffr_influence_df()]); model-based intervals, and fits of refund 0.1-40
#' and earlier whose model-based bread is lost, use Gaussian ones.
#'
#' @param object Fitted pffr model.
#' @param use_sandwich Do the standard errors come from the CL2 sandwich?
#' @param cluster Optional per-curve grouping, as in [coef.pffr()].
#' @returns List with `mode` (`"satterthwaite"` or `"z"`) and the influence
#'   object `core` (or `NULL`).
#' @keywords internal
pffr_crit_setup <- function(object, use_sandwich, cluster = NULL) {
  if (
    !isTRUE(use_sandwich) ||
      identical(pffr_canonicalize_cov(object)$format, "ancient")
  ) {
    return(list(mode = "z", core = NULL))
  }
  list(
    mode = "satterthwaite",
    core = pffr_influence(object, cluster = cluster)
  )
}

#' Pointwise critical values and reference df
#'
#' @param setup A [pffr_crit_setup()] result.
#' @param level Confidence level.
#' @param Xp Contrast matrix (`n_points x p`, full coefficient space), only
#'   needed for Satterthwaite critical values.
#' @returns A list with `crit` (scalar or per-point vector) and `df` (per-point
#'   vector; `Inf` for Gaussian critical values). Where the moment df is
#'   undefined (a contrast with zero sampling variance), the Gaussian critical
#'   value is used, with a message.
#' @keywords internal
pffr_pointwise_crit_values <- function(setup, level, Xp) {
  prob <- (1 + level) / 2
  n <- nrow(Xp)
  if (setup$mode == "z" || n == 0L) {
    return(list(crit = stats::qnorm(prob), df = rep(Inf, n)))
  }
  df <- pffr_influence_df(setup$core, Xp)$df
  if (any(!is.finite(df))) {
    pffr_note_undefined_df()
    df[!is.finite(df)] <- Inf
  }
  list(crit = stats::qt(prob, df), df = df)
}

#' Pointwise critical values for prediction contrasts
#'
#' Applies [pffr_crit_setup()] and [pffr_pointwise_crit_values()] to the rows
#' of a prediction matrix, so predictions get the same critical values as
#' [coef.pffr()] gives for the same contrasts. The df cost grows with the
#' number of rows and with the number of clusters (about 1 ms per row at
#' \eqn{G = 100}); for more than
#' `getOption("refund.pffr_satterthwaite_max_points", 1e4)` rows Gaussian
#' critical values are used, with a message.
#'
#' @param object Fitted pffr model.
#' @param X Prediction matrix (full coefficient space), one row per point.
#' @param level Confidence level.
#' @param use_sandwich Do the standard errors come from the CL2 sandwich?
#' @param cluster Optional per-curve grouping.
#' @returns List with `mode`, `crit` (scalar, or one per row) and `df` (one
#'   per row; `Inf` for Gaussian critical values).
#' @keywords internal
pffr_pointwise_crit <- function(
  object,
  X,
  level,
  use_sandwich,
  cluster = NULL
) {
  setup <- pffr_crit_setup(object, use_sandwich, cluster = cluster)
  max_points <- getOption("refund.pffr_satterthwaite_max_points", 1e4)
  if (setup$mode == "satterthwaite" && nrow(X) > max_points) {
    message(
      "Using Gaussian critical values for ",
      nrow(X),
      " prediction points: Satterthwaite df are computed for at most ",
      max_points,
      " points (option refund.pffr_satterthwaite_max_points)."
    )
    setup <- list(mode = "z", core = NULL)
  }
  cv <- pffr_pointwise_crit_values(setup, level, as.matrix(X))
  list(mode = setup$mode, crit = cv$crit, df = cv$df)
}

# Undefined working-model moment df can be detected at several sites per
# coef.pffr() call (once per smooth term and once for the parametric
# coefficients). pffr_begin_undefined_df() opens a collection window for the
# duration of one call; inside it the sites only record the fact and
# pffr_end_undefined_df() emits a single message on exit. Outside a window the
# message is emitted immediately.
pffr_df_note_state <- new.env(parent = emptyenv())
pffr_df_note_state$active <- FALSE
pffr_df_note_state$seen <- FALSE

PFFR_UNDEFINED_DF_MSG <- paste0(
  "Satterthwaite df are undefined at some evaluation points (zero sampling ",
  "variance); Gaussian critical values are used there."
)

#' Note (or record) that a working-model moment df was undefined
#'
#' @returns `NULL`, invisibly. Called for the side effect.
#' @keywords internal
pffr_note_undefined_df <- function() {
  if (isTRUE(pffr_df_note_state$active)) {
    pffr_df_note_state$seen <- TRUE
    return(invisible(NULL))
  }
  message(PFFR_UNDEFINED_DF_MSG)
  invisible(NULL)
}

#' Open an undefined-df collection window
#'
#' Pair with [pffr_end_undefined_df()] via `on.exit()`. Nested windows keep
#' the outermost one in charge.
#'
#' @returns `TRUE` if this call opened the window, `FALSE` if one was already
#'   open.
#' @keywords internal
pffr_begin_undefined_df <- function() {
  if (isTRUE(pffr_df_note_state$active)) {
    return(FALSE)
  }
  pffr_df_note_state$active <- TRUE
  pffr_df_note_state$seen <- FALSE
  TRUE
}

#' Close an undefined-df collection window, noting at most once
#'
#' @param opened The value returned by the matching [pffr_begin_undefined_df()].
#' @returns `NULL`, invisibly. Called for the side effect.
#' @keywords internal
pffr_end_undefined_df <- function(opened) {
  if (!isTRUE(opened)) {
    return(invisible(NULL))
  }
  seen <- isTRUE(pffr_df_note_state$seen)
  pffr_df_note_state$active <- FALSE
  pffr_df_note_state$seen <- FALSE
  if (seen) message(PFFR_UNDEFINED_DF_MSG)
  invisible(NULL)
}

#' Inform once per session that NCV-centred intervals are not for inference
#'
#' @param object A fitted pffr model.
#' @returns `invisible(NULL)`.
#' @keywords internal
pffr_inform_ncv_intervals <- function(object) {
  if (!identical(object$method, "NCV")) return(invisible(NULL))
  pffr_inform_once(
    "ncv_intervals",
    paste0(
      "Intervals around NCV estimates undercover: use the NCV fit for the ",
      "shape of coefficient surfaces, and a REML fit of the same model for ",
      "intervals (see ?pffr)."
    )
  )
}

#' Why the CL2 sandwich is unavailable for a fit, if it is
#'
#' @param m The fitted gam (or gamm) object.
#' @param algorithm Algorithm name.
#' @param cluster_id Cluster vector of the fit (one entry per fitted
#'   observation).
#' @returns `NULL` if CL2 is available, otherwise a short reason.
#' @keywords internal
pffr_cl2_unavailable <- function(m, algorithm, cluster_id) {
  if (algorithm %in% c("gamm4", "gamm")) {
    return(paste0(
      "the sandwich is not available for algorithm = \"",
      algorithm,
      "\""
    ))
  }
  if (pffr_score_kind(m$family) == "custom") {
    return(paste0(
      "no cluster-robust covariance is available for family '",
      m$family$family,
      "'"
    ))
  }
  if (length(unique(cluster_id)) < 2L) {
    return("the CL2 sandwich needs at least two curves or clusters")
  }
  NULL
}

#' Does a fit have a binary response?
#'
#' @param m A fitted gam object.
#' @returns A single logical.
#' @keywords internal
pffr_is_binary_fit <- function(m) {
  fam <- tolower(as.character(m$family$family))
  fam %in%
    c("binomial", "quasibinomial") &&
    all(as.vector(m$y) %in% c(0, 1)) &&
    all(as.vector(m$prior.weights %||% 1) == 1)
}

#' Store the CL2 sandwich covariance on a fitted pffr model
#'
#' Computes the curve-clustered CL2 covariance (exact leverage adjustment,
#' Bayesian form) from the fit's model-based bread and stores it in
#' `$pffr$Vsandwich` together with `$pffr$sandwich_info`. The model-based
#' `$Vp`/`$Vc`/`$Ve` are left untouched, so recomputing a sandwich later never
#' double-applies the correction. Emits the small-\eqn{G} warning and the
#' one-time message for binary responses.
#'
#' @param m Fitted model (gam/bam with `$pffr` metadata).
#' @returns Model with the robust covariance stored in `$pffr$Vsandwich`.
#' @keywords internal
apply_sandwich_correction <- function(m) {
  if (!is.environment(m$pffr$Vsandwich_cache))
    m$pffr$Vsandwich_cache <- new.env(parent = emptyenv())
  core <- pffr_influence(m)
  Vsw <- gam_sandwich_cluster_cl2(
    pffr_model_based_gam(m),
    cluster_id = NULL,
    influence = core
  )
  G <- core$G
  m$pffr$Vsandwich <- Vsw
  m$pffr$sandwich_info <- list(
    type = "cl2",
    cluster_var = m$pffr$cluster,
    G = G,
    cl2_adjustment = "exact",
    n_adjusted = attr(Vsw, "n_adjusted") %||% 0L,
    max_leverage = attr(Vsw, "max_leverage") %||% NA_real_,
    min_block_eig = attr(Vsw, "min_block_eig") %||% NA_real_,
    max_block_kappa = attr(Vsw, "max_block_kappa") %||% NA_real_,
    # Hat-invariant monitors (study LB, claim P-LB5): make a degenerate fit
    # visible from the fit object.
    max_obs_leverage = attr(Vsw, "max_obs_leverage") %||% NA_real_,
    min_obs_leverage = attr(Vsw, "min_obs_leverage") %||% NA_real_,
    min_hat_eig = attr(Vsw, "min_hat_eig") %||% NA_real_,
    min_block_eig_rel = attr(Vsw, "min_block_eig_rel") %||% NA_real_,
    hat_invariant_violation = attr(Vsw, "hat_invariant_violation"),
    version = as.character(utils::packageVersion("refund")),
    inference_core_version = core$version,
    cluster_rank = core$diagnostics$rank,
    storage_format = PFFR_COV_STORAGE_FORMAT
  )
  # Small-cluster-count guard: warn below G = 40, stating what the simulation
  # evidence covers. Fires once, at fit time -- pffr_vcov() (coef/predict/plot)
  # never calls this function again.
  if (G < 40) {
    warning(warningCondition(
      sprintf(
        paste0(
          "Only G = %d clusters: CL2 intervals with Satterthwaite critical ",
          "values were evaluated down to G = 20 (near-nominal for ",
          "coefficient surfaces and fitted means; intercepts of binary ",
          "models can undercover). Treat intervals as approximate."
        ),
        G
      ),
      class = "pffr_small_G_warning"
    ))
  }
  if (pffr_is_binary_fit(m)) {
    pffr_inform_once(
      "binary_cl2",
      paste0(
        "Binary response: CL2 intervals for the functional intercept and for ",
        "coefficient functions of scalar covariates can undercover, also with ",
        "Satterthwaite critical values. Treat them as approximate."
      )
    )
  }
  m
}


# =============================================================================
# pffr() orchestration helpers
# =============================================================================

#' Prepare data, formula, and call for pffr()
#'
#' Internal helper that performs input validation, formula parsing and
#' transformation, construction of the mgcv data object, and assembly of the
#' mgcv call. The returned list contains everything required to fit and
#' post-process a pffr model.
#'
#' @param call The matched call from pffr().
#' @param formula The original pffr formula.
#' @param yind The y-index argument (may be `NULL` if missing).
#' @param yind_missing Logical, whether `yind` was missing in the original call.
#' @param yind_expr The yind expression from the original call (for naming).
#' @param data The data argument (list/data.frame).
#' @param ydata Sparse response data (or `NULL`).
#' @param algorithm User-specified algorithm (may be NA).
#' @param method The requested mgcv method (character).
#' @param tensortype Tensor product type (symbol).
#' @param bs_yindex Basis specification for y-index.
#' @param bs_int Basis specification for functional intercept.
#' @param dots The list of additional arguments (...).
#' @returns A list with preparation outputs, including `new_call` and
#'   `pffr_data`.
#' @keywords internal
pffr_prepare <- function(
  call,
  formula,
  yind,
  yind_missing,
  yind_expr,
  data,
  ydata,
  algorithm,
  method,
  tensortype,
  bs_yindex,
  bs_int,
  dots
) {
  validated <- pffr_validate_dots(call, algorithm, dots, check_ar = TRUE)
  dots <- validated$dots
  use_ar <- validated$use_ar
  gaulss <- validated$gaulss

  pffr_validate_ydata(ydata)

  parsed <- parse_pffr_model_formula(formula, data, ydata)
  tf <- parsed$tf
  term_strings <- parsed$term_strings
  terms <- parsed$terms
  frml_env <- parsed$frml_env
  where_specials <- parsed$where_specials
  response_name <- parsed$response_name

  formula_env <- new.env()
  eval_env <- data

  dims <- pffr_get_dimensions(ydata, data, response_name, frml_env, eval_env)
  is_sparse <- dims$is_sparse
  nobs <- dims$nobs
  nyindex <- dims$nyindex
  ntotal <- dims$ntotal

  if (is_sparse) {
    yind_info <- pffr_setup_yind_sparse(ydata)
    yind <- yind_info$yind
    yind_name <- yind_info$yind_name
    nyindex <- yind_info$nyindex
  } else {
    yind_info <- pffr_setup_yind_dense(
      yind = if (yind_missing) NULL else yind,
      yind_missing = yind_missing,
      nyindex = nyindex,
      where_specials = where_specials,
      terms = terms,
      eval_env = eval_env,
      frml_env = frml_env,
      data = data,
      yind_expr = yind_expr
    )
    yind <- yind_info$yind
    yind_name <- yind_info$yind_name
  }

  method_missing <- !("method" %in% names(call))
  algorithm <- pffr_configure_algorithm(
    algorithm = algorithm,
    ntotal = ntotal,
    call = call,
    where_specials = where_specials,
    use_ar = use_ar
  )

  if (use_ar) {
    if (method_missing) {
      call$method <- "fREML"
    } else if (!identical(method, "fREML")) {
      stop(
        "Autocorrelated errors via `rho` require `method = \"fREML\"` (see ?mgcv::bam)."
      )
    }
  } else if (as.character(algorithm) == "bam" && method_missing) {
    call$method <- "fREML"
  }

  if (as.character(algorithm) == "bam" && !("chunk.size" %in% names(call))) {
    call$chunk.size <- 10000
  }

  resp_info <- pffr_setup_response(
    is_sparse = is_sparse,
    ydata = ydata,
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

  new_term_strings <- attr(tf, "term.labels")

  if (parsed$has_intercept) {
    int_result <- transform_intercept_term(
      yindex_vec_name,
      bs_int,
      yind_name
    )
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

  if (length(c(where_specials$ff, where_specials$sff))) {
    ff_terms <- lapply(
      terms[c(where_specials$ff, where_specials$sff)],
      \(x) eval(x, envir = eval_env, enclos = frml_env)
    )
    new_term_strings[c(where_specials$ff, where_specials$sff)] <- sapply(
      ff_terms,
      \(x) safeDeparse(x$call)
    )
    pffr_process_ff_terms(
      ff_terms = ff_terms,
      yind_vec = yind_vec,
      obs_indices = obs_indices,
      formula_env = formula_env
    )
  } else {
    ff_terms <- NULL
  }

  if (length(where_specials$ffpc)) {
    ffpc_terms <- lapply(
      terms[where_specials$ffpc],
      \(x) eval(x, envir = eval_env, enclos = frml_env)
    )
    new_term_strings[where_specials$ffpc] <- pffr_process_ffpc_terms(
      ffpc_terms = ffpc_terms,
      obs_indices = obs_indices,
      yindex_vec_name = yindex_vec_name,
      formula_env = formula_env
    )
    ffpc_terms <- lapply(ffpc_terms, \(x) x[names(x) != "data"])
  } else {
    ffpc_terms <- NULL
  }

  if (length(where_specials$pcre)) {
    pcre_terms <- lapply(
      terms[where_specials$pcre],
      \(x) eval(x, envir = eval_env, enclos = frml_env)
    )
    new_term_strings[where_specials$pcre] <- pffr_process_pcre_terms(
      pcre_terms = pcre_terms,
      is_sparse = is_sparse,
      yind = yind,
      yind_vec = yind_vec,
      nyindex = nyindex,
      nobs = nobs,
      obs_indices = obs_indices,
      formula_env = formula_env
    )
  } else {
    pcre_terms <- NULL
  }

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
            bs_yindex,
            tensortype,
            algorithm
          )
      )
  }

  if (length(where_specials$par)) {
    new_term_strings[where_specials$par] <- sapply(
      terms[where_specials$par],
      \(x) transform_par_term(x, yindex_vec_name, bs_yindex)
    )
  }

  if (!is.null(data)) {
    var_env <- list2env(data, envir = new.env(parent = frml_env))
  } else {
    var_env <- frml_env
  }
  pffr_expand_variables(
    terms = terms,
    where_specials = where_specials,
    is_sparse = is_sparse,
    nyindex = nyindex,
    nobs = nobs,
    obs_indices = obs_indices,
    eval_env = var_env,
    formula_env = formula_env
  )

  new_formula <- build_mgcv_formula(
    response_name = response_name,
    intercept_string = if (add_f_int) int_string else NULL,
    term_strings = new_term_strings,
    has_intercept = add_f_int,
    formula_env = formula_env
  )

  if (gaulss) {
    if (is.null(dots$varformula)) {
      dots$varformula <- formula(paste(
        "~",
        safeDeparse(as.call(c(
          as.name("s"),
          x = as.symbol(yindex_vec_name),
          bs_int
        )))
      ))
    }
    environment(dots$varformula) <- formula_env
    new_formula <- list(new_formula, dots$varformula)
  }

  if (use_ar) {
    pffr_build_ar_start(response_name, formula_env, obs_indices)
  }

  pffr_data <- build_mgcv_data(formula_env)
  new_call <- pffr_build_call(
    call = call,
    algorithm = algorithm,
    new_formula = new_formula,
    pffr_data = pffr_data,
    dots = dots,
    use_ar = use_ar,
    nobs = nobs,
    nyindex = nyindex,
    obs_indices = obs_indices
  )
  new_call$data <- pffr_data
  if (use_ar) {
    new_call$AR.start <- pffr_data$AR.start
  }

  list(
    call = call,
    algorithm = algorithm,
    use_ar = use_ar,
    gaulss = gaulss,
    is_sparse = is_sparse,
    nobs = nobs,
    nyindex = nyindex,
    ntotal = ntotal,
    yind = yind,
    yind_name = yind_name,
    response_name = response_name,
    terms = terms,
    where_specials = where_specials,
    term_strings = term_strings,
    new_term_strings = new_term_strings,
    add_f_int = add_f_int,
    int_string = int_string,
    ff_terms = ff_terms,
    ffpc_terms = ffpc_terms,
    pcre_terms = pcre_terms,
    missing_indices = missing_indices,
    ydata = ydata,
    yindex_vec_name = yindex_vec_name,
    formula_env = formula_env,
    pffr_data = pffr_data,
    new_call = new_call,
    dots = dots,
    new_formula = new_formula
  )
}
