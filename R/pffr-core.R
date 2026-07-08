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
  newcall$dof_correction <- newcall$edf_type <- NULL
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
  sandwich,
  dof_correction = "none",
  edf_type = "trace"
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
    dof_correction = dof_correction,
    edf_type = edf_type,
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
#'   `NULL` (default) clusters by curve. Only supported for dense responses.
#' @returns Integer vector of length equal to the number of fitted rows.
#' @keywords internal
build_cluster_id <- function(pffr_meta, cluster = NULL) {
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
    cluster_id <- pffr_meta$ydata$.obs
  } else {
    cluster_id <- rep(seq_len(pffr_meta$nobs), each = pffr_meta$nyindex)
  }
  if (!is.null(pffr_meta$missing_indices)) {
    cluster_id <- cluster_id[-pffr_meta$missing_indices]
  }
  cluster_id
}

#' Compute per-observation scores for gaulss family
#'
#' gaulss uses `tau = 1/sigma` with logb link and defines `family$sandwich`,
#' which prevents the generic GLM score path. Scores are computed analytically
#' from the log-likelihood and transformed to linear predictor space.
#'
#' @param b Fitted GAM object.
#' @param X Model matrix with `lpi` attribute.
#' @returns Score matrix (n_obs x p).
#' @keywords internal
compute_gaulss_scores <- function(b, X) {
  lpi <- attr(X, "lpi")
  n_obs <- length(b$y)

  mu <- b$fitted.values[seq_len(n_obs)]
  tau <- b$fitted.values[n_obs + seq_len(n_obs)]
  eta1 <- b$linear.predictors[seq_len(n_obs)]
  eta2 <- b$linear.predictors[n_obs + seq_len(n_obs)]

  r <- b$y - mu

  # Scores w.r.t. natural parameters, transformed to LP space via chain rule
  w1 <- (tau^2 * r) * b$family$linfo[[1]]$mu.eta(eta1)
  w2 <- (1 / tau - tau * r^2) * b$family$linfo[[2]]$mu.eta(eta2)

  pw <- b$prior.weights
  if (!is.null(pw) && any(pw != 1)) {
    w1 <- pw * w1
    w2 <- pw * w2
  }

  # Use += to handle possible overlapping lpi indices (shared coefficients)
  S <- matrix(0, nrow = n_obs, ncol = ncol(X))
  S[, lpi[[1]]] <- S[, lpi[[1]]] + w1 * X[, lpi[[1]], drop = FALSE]
  S[, lpi[[2]]] <- S[, lpi[[2]]] + w2 * X[, lpi[[2]], drop = FALSE]
  S
}

#' Symmetric matrix inverse square root with eigenvalue floor
#'
#' @param M Symmetric matrix.
#' @param tol Eigenvalue floor for numerical stability.
#' @returns Matrix inverse square root of `M`.
#' @keywords internal
sym_inv_sqrt <- function(M, tol = 1e-8) {
  M <- 0.5 * (M + t(M))
  ee <- eigen(M, symmetric = TRUE)
  vals <- pmax(ee$values, tol)
  ee$vectors %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(ee$vectors)
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
#' (prior-weighted) score weights from [compute_gaulss_scores()]. Both blocks
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

  # Same score components used by compute_gaulss_scores()
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

#' Assemble cluster-robust sandwich from score matrix
#'
#' Given a per-observation score matrix, aggregates by cluster and forms
#' \eqn{V_{CL} = f \cdot c \cdot V_p M_{CL} V_p + B_2} with HC1 correction
#' \eqn{c = G / (G - 1)} and an optional small-sample dof factor \eqn{f}
#' (default 1; see [compute_dof_factor()]).
#'
#' @param scores Per-observation score matrix (n_obs x p).
#' @param cluster_id Cluster membership vector.
#' @param Vp Bayesian posterior covariance (p x p).
#' @param B2 Bias correction matrix (p x p, or scalar 0).
#' @param dof_factor Scalar multiplier on the meat (default 1). Used to apply
#'   the optional CR1 small-sample correction \eqn{(N-1)/(N-\mathrm{EDF})}.
#' @returns A p x p covariance matrix.
#' @keywords internal
#' Number of clusters, requiring at least two
#'
#' The cluster sandwich's G/(G-1) small-sample factor is undefined for a
#' single cluster.
#'
#' @param cluster_id Cluster membership vector.
#' @returns The number of distinct clusters.
#' @keywords internal
n_clusters_checked <- function(cluster_id) {
  G <- length(unique(cluster_id))
  if (G < 2) {
    stop("Need at least two clusters for cluster sandwich.", call. = FALSE)
  }
  G
}

assemble_cluster_sandwich <- function(
  scores,
  cluster_id,
  Vp,
  B2,
  dof_factor = 1,
  b2 = TRUE,
  center_scores = FALSE
) {
  G <- n_clusters_checked(cluster_id)
  U <- rowsum(scores, cluster_id)
  if (isTRUE(center_scores)) {
    # U_g^c = U_g - (sum_g U_g) / G; sum_g U_g = S theta-hat (the working
    # log-likelihood penalty gradient). Centering removes the rank-one penalty
    # direction from the meat (X6). Exact: crossprod(U - Ubar) =
    # crossprod(U) - (colSums U)(colSums U)' / G.
    Ubar <- colSums(U) / G
    U <- sweep(U, 2L, Ubar, "-")
  }
  meat <- crossprod(U)
  hc1 <- G / (G - 1)
  core <- dof_factor * hc1 * Vp %*% meat %*% Vp
  if (isTRUE(b2)) core + B2 else core
}

#' Compute the optional CR1 small-sample dof factor
#'
#' Returns the scalar \eqn{(N-1)/(N-\mathrm{EDF})} when `dof_correction = "edf"`,
#' or `1` when `dof_correction = "none"`. Here \eqn{N} is the number of fitted
#' observations (`length(cluster_id)`) and the effective degrees of freedom is
#' selected by `edf_type`: `"trace"` uses `sum(b$edf)` (= trace of the penalized
#' hat), `"edf2"` uses `sum(b$edf2)` (mgcv's bias-corrected EDF), and `"basis"`
#' uses the basis dimension (`length(b$coefficients)`).
#'
#' This is the textbook CR1 finite-sample factor (Stata's cluster-robust default
#' multiplies the meat by an analogous \eqn{(N-1)/(N-k)} term). It is applied
#' ONLY in the CR1 path ([gam_sandwich_cluster()]); the CL2 path
#' ([gam_sandwich_cluster_cl2()]) already performs a per-cluster leverage
#' correction targeting the same downward bias, so combining the two would
#' double-correct.
#'
#' @param b Fitted GAM object.
#' @param cluster_id Cluster membership vector (length = number of fitted obs).
#' @param dof_correction `"none"` (factor 1) or `"edf"`.
#' @param edf_type Which EDF to use: `"trace"`, `"edf2"`, or `"basis"`.
#' @returns A finite positive scalar.
#' @keywords internal
compute_dof_factor <- function(
  b,
  cluster_id,
  dof_correction = c("none", "edf"),
  edf_type = c("trace", "edf2", "basis")
) {
  dof_correction <- match.arg(dof_correction)
  if (dof_correction == "none") {
    return(1)
  }
  edf_type <- match.arg(edf_type)

  N <- length(cluster_id)
  if (!is.finite(N) || N < 2) {
    stop(
      "dof_correction = \"edf\" requires at least N = 2 fitted observations.",
      call. = FALSE
    )
  }
  G <- n_clusters_checked(cluster_id)

  edf <- switch(
    edf_type,
    trace = if (!is.null(b$edf)) sum(b$edf) else NA_real_,
    edf2 = if (!is.null(b$edf2)) sum(b$edf2) else NA_real_,
    basis = length(b$coefficients)
  )
  if (!is.finite(edf) || edf < 0 || edf >= N) {
    stop(
      "dof_correction = \"edf\": effective degrees of freedom (edf_type = \"",
      edf_type,
      "\") evaluated to ",
      format(edf),
      ", which must be finite and in [0, N) with N = ",
      N,
      ". (Is 'edf2' available for this fit?)",
      call. = FALSE
    )
  }
  (N - 1) / (N - edf)
}

#' Cluster-robust sandwich covariance estimator
#'
#' Computes a cluster-robust (CR1) sandwich covariance matrix for a fitted
#' GAM, clustering by curve. Handles both heteroskedasticity and within-curve
#' autocorrelation, requiring only that curves are independent.
#'
#' @param b Fitted GAM object (must not have class `"pffr"`).
#' @param cluster_id Integer vector of length `nrow(model.matrix(b))` mapping
#'   each vectorized observation to its curve.
#' @param freq If `TRUE`, use frequentist sandwich (`B2 = 0`).
#'   If `FALSE` (default), use Bayesian sandwich (`B2 = Vp - Ve`).
#' @param dof_correction Optional small-sample correction multiplying the meat:
#'   `"none"` (default, factor 1 = current behavior) or `"edf"`, which applies
#'   the textbook CR1 factor \eqn{(N-1)/(N-\mathrm{EDF})}. This is OFF by default
#'   and is applied only here (CR1), never in the CL2 path, which already
#'   corrects per-cluster leverage; see [compute_dof_factor()].
#' @param edf_type Which effective degrees of freedom the `"edf"` correction
#'   uses: `"trace"` (default, `sum(b$edf)`), `"edf2"`, or `"basis"`.
#' @param b2 Internal ablation switch (default `TRUE` = current behavior). When
#'   `FALSE`, drop the additive Bayesian smoothing-bias term \eqn{B_2 = V_p -
#'   V_e}, returning \eqn{c\,V_p (\sum_g U_g U_g^\top) V_p} only. Not a
#'   user-facing option; used by the B2-ablation experiment (X5).
#' @param center_scores Internal ablation switch (default `FALSE` = current
#'   behavior). When `TRUE`, center the per-cluster score sums
#'   \eqn{U_g^c = U_g - (\sum_g U_g)/G} before forming the meat (X6). Since
#'   \eqn{\sum_g U_g = S\hat\theta} (the working-log-likelihood penalty
#'   gradient), this removes the rank-one penalty direction from the meat.
#' @returns A p x p covariance matrix.
#' @keywords internal
gam_sandwich_cluster <- function(
  b,
  cluster_id,
  freq = FALSE,
  dof_correction = c("none", "edf"),
  edf_type = c("trace", "edf2", "basis"),
  b2 = TRUE,
  center_scores = FALSE
) {
  dof_correction <- match.arg(dof_correction)
  edf_type <- match.arg(edf_type)
  dof_factor <- compute_dof_factor(b, cluster_id, dof_correction, edf_type)

  B2 <- if (freq) 0 else b$Vp - b$Ve
  X <- model.matrix(b)

  if (b$family$family == "gaulss") {
    scores <- compute_gaulss_scores(b, X)
    return(assemble_cluster_sandwich(
      scores,
      cluster_id,
      b$Vp,
      B2,
      dof_factor = dof_factor,
      b2 = b2,
      center_scores = center_scores
    ))
  }

  # Families that define family$sandwich (e.g. multinom) use custom
  # score computation — cluster aggregation not yet implemented for these.
  if (!is.null(b$family$sandwich)) {
    warning(
      "Cluster-robust sandwich not yet implemented for family '",
      b$family$family,
      "'. Falling back to observation-level HC sandwich via mgcv::vcov.gam().",
      call. = FALSE
    )
    return(mgcv::vcov.gam(b, sandwich = TRUE, freq = freq))
  }

  # Standard GLM case: per-observation scores via general score weight
  # w = pw * (d mu/d eta) * (y - mu) / (phi * V(mu)); prior weights enter the
  # estimating equation (the CL2 and gaulss paths already include them).
  mu <- b$fitted.values
  pw <- b$prior.weights %||% 1
  scores <- pw *
    b$family$mu.eta(b$linear.predictors) *
    (b$y - mu) /
    (b$sig2 * b$family$variance(mu)) *
    X
  assemble_cluster_sandwich(
    scores,
    cluster_id,
    b$Vp,
    B2,
    dof_factor = dof_factor,
    b2 = b2,
    center_scores = center_scores
  )
}

#' Cluster-robust CL2 sandwich covariance estimator
#'
#' Computes a cluster-robust covariance matrix with a Bell-McCaffrey style
#' leverage adjustment (`CL2` in this package), clustering by curve.
#'
#' @param b Fitted GAM object (must not have class `"pffr"`).
#' @param cluster_id Integer vector mapping each vectorized observation to a
#'   curve.
#' @param freq If `TRUE`, use frequentist sandwich (`B2 = 0`).
#'   If `FALSE` (default), use Bayesian sandwich (`B2 = Vp - Ve`).
#' @param tol Eigenvalue floor for numerical stability.
#' @param leverage_cap Cap for cluster leverage eigenvalues (< 1).
#' @param b2 Internal ablation switch (default `TRUE` = current behavior). When
#'   `FALSE`, drop the additive Bayesian smoothing-bias term \eqn{B_2 = V_p -
#'   V_e} (X5).
#' @param center_scores Internal ablation switch (default `FALSE` = current
#'   behavior). When `TRUE`, center the leverage-adjusted per-cluster
#'   contributions \eqn{U_g^c = U_g - (\sum_g U_g)/G} before forming the meat
#'   (X6).
#' @returns A p x p covariance matrix with attributes `n_capped_clusters` and
#'   `max_leverage` for the CL2 leverage diagnostic.
#' @keywords internal
gam_sandwich_cluster_cl2 <- function(
  b,
  cluster_id,
  freq = FALSE,
  tol = 1e-8,
  leverage_cap = 0.999,
  b2 = TRUE,
  center_scores = FALSE
) {
  if (!is.finite(leverage_cap) || leverage_cap <= 0 || leverage_cap >= 1) {
    stop("`leverage_cap` must be in (0, 1).", call. = FALSE)
  }

  fam <- tolower(as.character(b$family$family))

  # Families with custom family$sandwich (e.g. multinom) use custom
  # score computation — CL2 cluster leverage correction is not implemented yet.
  if (fam != "gaulss" && !is.null(b$family$sandwich)) {
    warning(
      "CL2 sandwich not yet implemented for family '",
      b$family$family,
      "'. Falling back to observation-level HC sandwich via mgcv::vcov.gam().",
      call. = FALSE
    )
    return(mgcv::vcov.gam(b, sandwich = TRUE, freq = freq))
  }

  work <- if (fam == "gaulss") {
    build_cl2_working_gaulss(b, cluster_id)
  } else {
    build_cl2_working_standard(b, cluster_id)
  }

  Xw <- work$Xw
  z <- work$z
  cluster_id_work <- work$cluster_id

  G <- n_clusters_checked(cluster_id_work)
  groups <- unique(cluster_id_work)

  Vp <- b$Vp
  B2 <- if (freq) 0 else b$Vp - b$Ve
  p <- ncol(Xw)
  meat <- matrix(0, nrow = p, ncol = p)
  Usum <- numeric(p)
  n_capped_clusters <- 0L
  max_leverage <- NA_real_

  for (g in groups) {
    idx <- which(cluster_id_work == g)
    Xwg <- Xw[idx, , drop = FALSE]
    zg <- z[idx]

    Hgg <- Xwg %*% Vp %*% t(Xwg)
    Hgg <- 0.5 * (Hgg + t(Hgg))

    ee_H <- eigen(Hgg, symmetric = TRUE)
    cluster_max_leverage <- max(ee_H$values, na.rm = TRUE)
    if (is.finite(cluster_max_leverage)) {
      max_leverage <- if (is.na(max_leverage)) {
        cluster_max_leverage
      } else {
        max(max_leverage, cluster_max_leverage)
      }
    }
    if (any(ee_H$values > leverage_cap, na.rm = TRUE)) {
      n_capped_clusters <- n_capped_clusters + 1L
      ee_H$values <- pmin(ee_H$values, leverage_cap)
      Hgg <- ee_H$vectors %*%
        diag(ee_H$values, nrow = length(ee_H$values)) %*%
        t(ee_H$vectors)
    }

    Mg <- diag(length(idx)) - Hgg
    Ag <- sym_inv_sqrt(Mg, tol = tol)
    Ug <- crossprod(Xwg, Ag %*% zg)
    meat <- meat + Ug %*% t(Ug)
    Usum <- Usum + as.vector(Ug)
  }

  if (isTRUE(center_scores)) {
    # Exact centering: sum_g (U_g - Ubar)(U_g - Ubar)' = meat - Usum Usum'/G
    # with Ubar = Usum / G (X6).
    meat <- meat - tcrossprod(Usum) / G
  }

  hc1 <- G / (G - 1)
  V <- hc1 * Vp %*% meat %*% Vp
  if (isTRUE(b2)) V <- V + B2
  V <- 0.5 * (V + t(V))
  attr(V, "n_capped_clusters") <- n_capped_clusters
  attr(V, "max_leverage") <- max_leverage
  if (n_capped_clusters > 0) {
    max_leverage_label <- if (is.finite(max_leverage)) {
      sprintf("%.3f", max_leverage)
    } else {
      "NA"
    }
    warning(
      sprintf(
        paste0(
          "CL2 leverage adjustment hit the leverage cap %.3f in %d of %d ",
          "clusters (max pre-cap eigenvalue %s). CL2 can be unreliable ",
          "with small G or saturated per-cluster leverage; consider ",
          "sandwich = \"cluster\" as the safer choice."
        ),
        leverage_cap,
        n_capped_clusters,
        G,
        max_leverage_label
      ),
      call. = FALSE
    )
  }
  V
}

#' Resolve a pffr fit's covariance matrices to a canonical representation
#'
#' The current storage contract (format 2) keeps `Vp`/`Vc`/`Ve` model-based
#' ALWAYS and stores the robust covariance in `object$pffr$Vsandwich`
#' (+ `Vsandwich_freq`). Older fits (format 1) instead overwrote `Vp`/`Vc`/`Ve`
#' with the robust matrices and stashed the model-based originals in
#' `object$pffr$model_cov`. This helper resolves either layout into a single
#' shape, emitting a session-scoped one-time back-compat warning for old-format
#' objects. All covariance recomputation must use the model-based matrices as
#' the (penalized) bread, or the sandwich is applied on top of itself.
#'
#' @param object A fitted pffr model.
#' @returns A list with `model` (list of model-based `Vp`/`Vc`/`Ve`),
#'   `fit_type` (the fit-time sandwich type), `Vsandwich`/`Vsandwich_freq` (the
#'   fit-time robust matrices or `NULL`), and `format`
#'   (`"new"`/`"old"`/`"ancient"`).
#' @keywords internal
pffr_canonicalize_cov <- function(object) {
  meta <- object$pffr
  fit_type <- normalize_sandwich_type(meta$sandwich)

  # Format 2 (current): model-based Vp/Vc/Ve, robust in $pffr$Vsandwich.
  # The cov_format stamp is set at metadata-build time, so this branch is also
  # taken during fitting, before apply_sandwich_correction() has stored
  # Vsandwich (Vp/Vc/Ve are model-based there too).
  if (
    isTRUE((meta$cov_format %||% 0L) >= 2L) ||
      !is.null(meta$sandwich_info) ||
      !is.null(meta$Vsandwich)
  ) {
    return(list(
      model = list(Vp = object$Vp, Vc = object$Vc, Ve = object$Ve),
      fit_type = meta$sandwich_info$type %||% fit_type,
      Vsandwich = meta$Vsandwich,
      Vsandwich_freq = meta$Vsandwich_freq,
      format = "new"
    ))
  }

  # Plain fit: no sandwich applied, Vp/Vc/Ve are model-based.
  if (identical(fit_type, "none")) {
    return(list(
      model = list(Vp = object$Vp, Vc = object$Vc, Ve = object$Ve),
      fit_type = "none",
      Vsandwich = NULL,
      Vsandwich_freq = NULL,
      format = "new"
    ))
  }

  # Format 1 (old): Vp/Vc/Ve overwritten with robust, model-based in model_cov.
  mc <- meta$model_cov
  if (!is.null(mc)) {
    pffr_warn_once(
      "oldformat_modelcov",
      paste0(
        "This pffr fit was created by an older refund version: its $Vp/$Vc/$Ve ",
        "hold the sandwich-adjusted (robust) covariance and the model-based ",
        "matrices are stashed in $pffr$model_cov. Reading it in ",
        "backward-compatible mode. Call pffr_upgrade_fit() to convert it to ",
        "the current storage format ($Vp model-based, robust in ",
        "$pffr$Vsandwich)."
      )
    )
    return(list(
      model = list(Vp = mc$Vp, Vc = mc$Vc, Ve = mc$Ve),
      fit_type = fit_type,
      Vsandwich = object$Vp,
      Vsandwich_freq = object$Ve,
      format = "old"
    ))
  }

  # Ancient: robust matrices, model-based irrecoverably lost.
  pffr_warn_once(
    "oldformat_nostash",
    paste0(
      "This pffr fit was created with sandwich = \"",
      fit_type,
      "\" by a very old refund version and does not carry the model-based ",
      "covariance matrices. Model-based uncertainty cannot be recovered and ",
      "any recomputed sandwich would double-apply the correction; refit the ",
      "model with the current refund version."
    )
  )
  list(
    model = list(Vp = object$Vp, Vc = object$Vc, Ve = object$Ve),
    fit_type = fit_type,
    Vsandwich = object$Vp,
    Vsandwich_freq = object$Ve,
    format = "ancient"
  )
}

#' Return a fit's underlying gam with model-based covariance matrices
#'
#' Strips the `"pffr"` class and guarantees `Vp`/`Vc`/`Ve` are the model-based
#' matrices, so downstream mgcv machinery and the sandwich estimators use the
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

#' Back-compat shim: restore model-based covariance on a fit
#'
#' Retained for backward compatibility. Sets `Vp`/`Vc`/`Ve` to the model-based
#' matrices (a no-op for current-format fits, whose `Vp`/`Vc`/`Ve` are already
#' model-based; a restore-from-stash for old-format fits). See
#' [pffr_canonicalize_cov()].
#'
#' @param object A fitted pffr model (or its stripped gam version with the
#'   `$pffr` metadata still attached).
#' @returns The object with model-based `Vp`/`Vc`/`Ve`.
#' @keywords internal
restore_model_cov <- function(object) {
  canon <- pffr_canonicalize_cov(object)
  object$Vp <- canon$model$Vp
  object$Vc <- canon$model$Vc
  object$Ve <- canon$model$Ve
  object
}

#' Compute a robust covariance matrix from a model-based gam
#'
#' Single dispatch point used by both fit-time correction and on-demand
#' recomputation. `b` must be a (stripped) gam whose `Vp`/`Vc`/`Ve` are the
#' model-based (penalized) matrices used as the sandwich bread.
#'
#' @param b Stripped gam object with model-based covariance matrices.
#' @param type One of `"cluster"`, `"cl2"`, `"hc"`, `"none"`.
#' @param cluster_id Integer cluster vector (for `"cluster"`/`"cl2"`), else
#'   `NULL`.
#' @param freq If `TRUE`, frequentist sandwich (`B2 = 0`).
#' @param dof_correction,edf_type CR1 small-sample dof options (`"cluster"`
#'   only).
#' @param b2 Internal ablation switch (default `TRUE`); drop the additive
#'   \eqn{B_2} term when `FALSE` (X5). Applies to `"cluster"`/`"cl2"` only.
#' @param center_scores Internal ablation switch (default `FALSE`); center the
#'   per-cluster score sums before the meat when `TRUE` (X6). Applies to
#'   `"cluster"`/`"cl2"` only.
#' @returns A covariance matrix (with CL2 leverage attributes for `type =
#'   "cl2"`).
#' @keywords internal
pffr_compute_sandwich <- function(
  b,
  type,
  cluster_id,
  freq = FALSE,
  dof_correction = "none",
  edf_type = "trace",
  b2 = TRUE,
  center_scores = FALSE
) {
  switch(
    type,
    cluster = gam_sandwich_cluster(
      b,
      cluster_id,
      freq = freq,
      dof_correction = dof_correction,
      edf_type = edf_type,
      b2 = b2,
      center_scores = center_scores
    ),
    cl2 = gam_sandwich_cluster_cl2(
      b,
      cluster_id,
      freq = freq,
      b2 = b2,
      center_scores = center_scores
    ),
    hc = mgcv::vcov.gam(b, sandwich = TRUE, freq = freq),
    none = if (freq) b$Ve else (b$Vc %||% b$Vp),
    stop("Unknown sandwich type: ", type, call. = FALSE)
  )
}

#' Resolve the covariance matrix for a pffr fit (single accessor)
#'
#' The one internal accessor through which every refund consumer
#' (`coef.pffr`, `plot.pffr`, `predict.pffr`, simultaneous bands) obtains a
#' covariance matrix, so that `$Vp`/`$Vc`/`$Ve` stay model-based and robust
#' matrices are never double-applied.
#'
#' @param object A fitted pffr model.
#' @param sandwich `NULL` (default) uses the fit-time choice; otherwise one of
#'   `"none"`/`"cluster"`/`"cl2"`/`"hc"` (or a legacy logical). `"none"` returns
#'   the genuinely model-based covariance; any other value that differs from the
#'   fit-time choice (or supplies a custom `cluster`) recomputes from the
#'   model-based bread.
#' @param freq If `TRUE`, return the frequentist covariance.
#' @param cluster Optional per-curve grouping forcing recomputation at a custom
#'   cluster level.
#' @param dof_correction,edf_type CR1 dof options; `NULL` inherits the fit's.
#' @param b2 Internal ablation switch (default `TRUE` = current behavior). When
#'   `FALSE`, drop the additive Bayesian smoothing-bias term \eqn{B_2} from the
#'   cluster/CL2 sandwich (X5). Non-default values always force a fresh
#'   recomputation and are never served from or written to the cache.
#' @param center_scores Internal ablation switch (default `FALSE` = current
#'   behavior). When `TRUE`, center the per-cluster score sums before forming
#'   the cluster/CL2 meat (X6). Same cache semantics as `b2`.
#' @returns A covariance matrix.
#' @keywords internal
pffr_vcov <- function(
  object,
  sandwich = NULL,
  freq = FALSE,
  cluster = NULL,
  dof_correction = NULL,
  edf_type = NULL,
  b2 = TRUE,
  center_scores = FALSE
) {
  # Ablation variants (X5/X6) bypass the fit-time and recompute caches entirely,
  # so they never overwrite or shadow the standard cached matrices.
  ablation <- !isTRUE(b2) || isTRUE(center_scores)
  canon <- pffr_canonicalize_cov(object)
  requested <- if (is.null(sandwich)) {
    canon$fit_type
  } else {
    normalize_sandwich_type(sandwich)
  }
  requested <- match.arg(requested, c("none", "cluster", "cl2", "hc"))

  if (identical(requested, "none")) {
    mb <- canon$model
    return(if (freq) mb$Ve else (mb$Vc %||% mb$Vp))
  }

  dof_correction <- dof_correction %||% (object$pffr$dof_correction %||% "none")
  edf_type <- edf_type %||% (object$pffr$edf_type %||% "trace")

  # Serve the cached fit-time robust matrix when the request matches it exactly.
  opts_match <- if (requested == "cluster") {
    identical(dof_correction, object$pffr$dof_correction %||% "none") &&
      (dof_correction == "none" ||
        identical(edf_type, object$pffr$edf_type %||% "trace"))
  } else {
    TRUE
  }
  if (
    !ablation &&
      is.null(cluster) &&
      identical(requested, canon$fit_type) &&
      opts_match &&
      !is.null(canon$Vsandwich)
  ) {
    return(
      if (freq) {
        canon$Vsandwich_freq %||% canon$Vsandwich
      } else {
        canon$Vsandwich
      }
    )
  }

  # Recompute from the model-based bread (safe by construction). Cache the
  # result in fit$pffr$Vsandwich_cache[[type]] -- an environment, so the cache
  # actually persists across coef()/predict() calls on the same stored object
  # (a plain list slot could not, under copy-on-modify). Keys extend the type
  # with the freq/dof options; custom `cluster` requests are never cached.
  cache <- object$pffr$Vsandwich_cache
  key <- if (!ablation && is.null(cluster)) {
    if (requested == "cluster" && (freq || dof_correction != "none")) {
      paste(requested, freq, dof_correction, edf_type, sep = "|")
    } else if (freq) {
      paste(requested, "freq", sep = "|")
    } else {
      requested
    }
  } else {
    NULL
  }
  if (!is.null(cache) && !is.null(key) && !is.null(cache[[key]])) {
    return(cache[[key]])
  }

  b <- object
  b$Vp <- canon$model$Vp
  b$Vc <- canon$model$Vc
  b$Ve <- canon$model$Ve
  class(b) <- setdiff(class(b), "pffr")
  cluster_id <- if (requested %in% c("cluster", "cl2")) {
    build_cluster_id(object$pffr, cluster = cluster)
  } else {
    NULL
  }
  V <- pffr_compute_sandwich(
    b,
    requested,
    cluster_id,
    freq = freq,
    dof_correction = dof_correction,
    edf_type = edf_type,
    b2 = b2,
    center_scores = center_scores
  )
  if (!is.null(cache) && !is.null(key)) {
    cache[[key]] <- V
  }
  V
}

#' Satterthwaite degrees of freedom for cluster-robust pointwise intervals
#'
#' Working-iid Satterthwaite degrees of freedom for a set of scalar contrasts
#' (the rows of `Xp`) --- the pointwise half of the Bell--McCaffrey procedure.
#' For a contrast \eqn{a} with Fisher-whitened per-cluster design
#' \eqn{\tilde X_g = \sqrt{W_g}\,X_g}, model-based penalized bread \eqn{V_p} and
#' CL2 leverage adjustment \eqn{A_g = (I - H_{gg})^{-1/2}} (identity for the CR1
#' path),
#' \deqn{q_g = A_g\,\tilde X_g\,(V_p a), \qquad
#'   \nu(a) = \frac{\left(\sum_g \lVert q_g\rVert^2\right)^2}
#'                 {\sum_g \lVert q_g\rVert^4},
#'   \qquad \mathrm{crit} = t_{1-\alpha/2,\,\nu}.}
#' Rationale: the robust variance of \eqn{a^\top\hat\theta} is
#' \eqn{c\sum_g (q_g^\top z_g)^2}; under the working model the per-cluster terms
#' are independent \eqn{\lVert q_g\rVert^2\chi^2_1}-type variables, and matching
#' the first two moments of their sum gives \eqn{\nu}. This working-iid shortcut
#' drops the same cross-cluster residual terms that the shipped
#' \eqn{(I-H_{gg})^{-1/2}} CL2 shortcut drops (paper Appendix C); it therefore
#' returns \eqn{\approx G} for a perfectly balanced design where the *exact*
#' Bell--McCaffrey df is \eqn{G-1}. The exact-BM df is future work (task X15).
#'
#' Vectorized over the rows of `Xp`: \eqn{M = V_p X_p^\top} (`p x n_points`) is
#' formed once; each cluster contributes \eqn{Q_g = A_g\,\tilde X_g\,M}
#' (`D_g x n_points`) and the per-column squared norms
#' \eqn{\lVert q_g\rVert^2 = \mathrm{colSums}(Q_g^2)} accumulate into
#' \eqn{s_2 = \sum_g \lVert q_g\rVert^2} and
#' \eqn{s_4 = \sum_g \lVert q_g\rVert^4}; then \eqn{\nu = s_2^2 / s_4}. Cost
#' \eqn{O(\sum_g D_g\, p\, n_{points})}.
#'
#' @param Xw Fisher-whitened per-observation design (`n_work x p`), from
#'   [build_cl2_working_standard()] / [build_cl2_working_gaulss()].
#' @param cluster_id Work-level cluster membership (length `n_work`).
#' @param Vp Model-based penalized bread (`p x p`).
#' @param Xp Contrast matrix, one row per evaluation point (`n_points x p`, full
#'   coefficient space).
#' @param use_cl2 Apply the CL2 leverage adjustment `A_g`? (`FALSE` = CR1 path,
#'   `A_g = I`.)
#' @param leverage_cap,tol CL2 leverage cap / eigenvalue floor (match the
#'   shipped CL2 sandwich in [gam_sandwich_cluster_cl2()]).
#' @returns A list with `df` (length `n_points`; `NA` at zero-variance
#'   contrasts, otherwise clamped to `[1, G]` up to rounding) and `G`.
#' @keywords internal
satterthwaite_df_kernel <- function(
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

#' Per-cluster whitening context for Satterthwaite degrees of freedom
#'
#' Builds --- once per `coef()` call --- the shared pieces the per-point
#' Satterthwaite df needs, so [satterthwaite_df_kernel()] can be applied to each
#' term's contrast matrix without rebuilding the whitened design. Uses the same
#' Fisher-whitened two-block / standard construction as the CL2 sandwich (so the
#' per-cluster hat trace equals the model EDF).
#'
#' @param object A fitted pffr model.
#' @param sandwich_type Resolved sandwich path; a whitening context is only
#'   built for `"cluster"` / `"cl2"`.
#' @param cluster Optional custom per-curve grouping (as in [pffr_vcov()]).
#' @param leverage_cap,tol CL2 leverage cap / eigenvalue floor.
#' @returns A list with `ok` (`FALSE` when the sandwich path is not
#'   cluster/CL2, or the family has no whitened score factorization here), and
#'   when `ok`: `Xw`, `cluster_id`, `Vp`, `use_cl2`, `G`, `leverage_cap`, `tol`.
#' @keywords internal
pffr_df_context <- function(
  object,
  sandwich_type,
  cluster = NULL,
  leverage_cap = 0.999,
  tol = 1e-8
) {
  type <- normalize_sandwich_type(sandwich_type)
  if (!type %in% c("cluster", "cl2")) {
    return(list(ok = FALSE, type = type))
  }
  b <- pffr_model_based_gam(object)
  fam <- tolower(as.character(b$family$family))
  # Families with a custom family$sandwich other than gaulss have no whitened
  # score factorization here (same restriction as the CL2 sandwich path).
  if (fam != "gaulss" && !is.null(b$family$sandwich)) {
    return(list(ok = FALSE, type = type))
  }
  cluster_id_curve <- build_cluster_id(object$pffr, cluster = cluster)
  work <- if (fam == "gaulss") {
    build_cl2_working_gaulss(b, cluster_id_curve)
  } else {
    build_cl2_working_standard(b, cluster_id_curve)
  }
  list(
    ok = TRUE,
    type = type,
    Xw = work$Xw,
    cluster_id = work$cluster_id,
    Vp = b$Vp,
    use_cl2 = identical(type, "cl2"),
    G = length(unique(work$cluster_id)),
    leverage_cap = leverage_cap,
    tol = tol
  )
}

#' Per-point Satterthwaite df for a set of contrasts, from a df context
#'
#' Thin wrapper around [satterthwaite_df_kernel()] that returns just the df
#' vector, or all-`NA` when the context carries no whitened design (`ok =
#' FALSE`), so callers can transparently fall back to the Gaussian reference.
#'
#' @param ctx A [pffr_df_context()] result.
#' @param Xp Contrast matrix (`n_points x p`, full coefficient space).
#' @returns Numeric vector of per-point df (length `nrow(Xp)`).
#' @keywords internal
pffr_df_from_context <- function(ctx, Xp) {
  if (!isTRUE(ctx$ok)) {
    return(rep(NA_real_, nrow(Xp)))
  }
  satterthwaite_df_kernel(
    Xw = ctx$Xw,
    cluster_id = ctx$cluster_id,
    Vp = ctx$Vp,
    Xp = Xp,
    use_cl2 = ctx$use_cl2,
    leverage_cap = ctx$leverage_cap,
    tol = ctx$tol
  )$df
}

#' Resolve the pointwise critical-value reference for [coef.pffr()]
#'
#' Maps the user's `crit` (`"auto"`/`"z"`/`"tG1"`/`"satterthwaite"`) to a
#' concrete reference. `"auto"` selects the per-point Satterthwaite reference
#' when the pointwise SEs come from a cluster-robust sandwich
#' (`"cluster"`/`"cl2"`) and the number of independent curves/clusters is
#' moderate (`G < 150`), and the Gaussian reference otherwise. An explicit
#' `"satterthwaite"` on a non-cluster covariance has no cluster leverage
#' structure to match and degrades to `"z"` with a warning; `"tG1"` (the
#' \eqn{t_{G-1}} reference, pointwise counterpart of the simultaneous
#' `ci_ref = "t"`) uses the curve/cluster count regardless of the covariance.
#'
#' @param crit One of `"auto"`, `"z"`, `"tG1"`, `"satterthwaite"`.
#' @param sandwich_type The resolved covariance type used for the SEs.
#' @param G Number of independent curves / clusters.
#' @returns One of `"z"`, `"tG1"`, `"satterthwaite"`.
#' @keywords internal
resolve_crit_reference <- function(crit, sandwich_type, G) {
  crit <- match.arg(crit, c("auto", "z", "tG1", "satterthwaite"))
  type <- normalize_sandwich_type(sandwich_type)
  is_cluster <- type %in% c("cluster", "cl2")
  if (crit == "auto") {
    return(if (is_cluster && is.finite(G) && G < 150) "satterthwaite" else "z")
  }
  if (crit == "satterthwaite" && !is_cluster) {
    warning(
      "crit = \"satterthwaite\" requires a cluster-robust covariance ",
      "(sandwich = \"cluster\" or \"cl2\"); using crit = \"z\" instead.",
      call. = FALSE
    )
    return("z")
  }
  crit
}

#' Penalty-direction and B2-magnitude diagnostics for a pffr sandwich
#'
#' Fit-level diagnostics that quantify how much of the cluster-robust sandwich
#' is driven by the penalty (uncentered-meat) direction and by the additive
#' \eqn{B_2} smoothing-bias term. Used by the X5/X6 ablation experiment; needs
#' the raw per-cluster scores, so it lives in the package rather than being
#' derivable from the covariance accessor alone.
#'
#' Definitions (CR1 per-cluster score sums \eqn{U_g}, total score
#' \eqn{S\hat\theta = \sum_g U_g}, meat \eqn{M = \sum_g U_g U_g^\top}, penalized
#' bread \eqn{V_p}, \eqn{c = G/(G-1)}, \eqn{B_2 = V_p - V_e}):
#' \itemize{
#'   \item `pen_share` \eqn{= \lVert S\hat\theta\rVert^2 / \sum_g \lVert
#'     U_g\rVert^2} --- the fraction of the total score energy that lies in the
#'     (rank-one) penalty direction the uncentered meat injects (X6).
#'   \item `fro_ratio` \eqn{= \lVert B_2\rVert_F / \lVert c\,V_p M V_p\rVert_F}
#'     --- the Frobenius-norm size of the additive \eqn{B_2} term relative to
#'     the score-sandwich core (X5), for the CR1 meat.
#' }
#' The centered-meat components (`Sbeta`, `meat`, `Vp`, `hc1`, `G`) are returned
#' so callers can verify the exact rank-one centering identity
#' \eqn{M_c = M - S\hat\theta\, S\hat\theta^\top / G}.
#'
#' @param object A fitted pffr model (any sandwich type; the model-based bread
#'   is used regardless).
#' @returns A list with `pen_share`, `fro_ratio`, `G`, `hc1`, `Sbeta` (total
#'   score), `meat` (CR1 meat \eqn{\sum_g U_g U_g^\top}), `Vp`, and `B2`.
#' @keywords internal
pffr_sandwich_shares <- function(object) {
  b <- pffr_model_based_gam(object)
  cluster_id <- build_cluster_id(object$pffr)
  G <- n_clusters_checked(cluster_id)
  X <- model.matrix(b)

  scores <- if (identical(tolower(as.character(b$family$family)), "gaulss")) {
    compute_gaulss_scores(b, X)
  } else {
    mu <- b$fitted.values
    pw <- b$prior.weights %||% 1
    pw *
      b$family$mu.eta(b$linear.predictors) *
      (b$y - mu) /
      (b$sig2 * b$family$variance(mu)) *
      X
  }

  U <- rowsum(scores, cluster_id)
  Sbeta <- colSums(scores)
  # sum_g ||U_g||^2 = sum of all squared entries of the G x p matrix U
  pen_share <- sum(Sbeta^2) / sum(U^2)

  Vp <- b$Vp
  B2 <- Vp - b$Ve
  hc1 <- G / (G - 1)
  meat <- crossprod(U)
  core <- hc1 * Vp %*% meat %*% Vp
  fro_ratio <- norm(B2, "F") / norm(core, "F")

  list(
    pen_share = pen_share,
    fro_ratio = fro_ratio,
    G = G,
    hc1 = hc1,
    Sbeta = Sbeta,
    meat = meat,
    Vp = Vp,
    B2 = B2
  )
}

#' Leave-one-cluster-out jackknife of the fitted linear predictor (core)
#'
#' Internal engine for [pffr_jackknife_se()]. Computes, per evaluation point,
#' the exact leave-one-cluster-out (LOCO) jackknife variance of the linear
#' predictor \eqn{\hat\eta(x) = X_p(x)^\top\hat\theta} via a Sherman--Morrison--
#' Woodbury (SMW) downdate of the penalized (weighted) normal equations at fixed
#' smoothing parameters and converged working weights.
#'
#' For each cluster \eqn{g} the deleted-cluster coefficient is
#' \deqn{\hat\theta_{(-g)} = \hat\theta -
#'   A^{-1} \tilde X_g^\top (I - H_{gg})^{-1} \tilde r_g,\quad
#'   A^{-1} = V_p/\sigma^2,\; H_{gg} = \tilde X_g A^{-1} \tilde X_g^\top,}
#' with \eqn{\tilde X_g = \sqrt{W_g}\,X_g} the Fisher-whitened training design of
#' cluster \eqn{g} (\eqn{W_i = w_i (\mathrm{d}\mu/\mathrm{d}\eta)_i^2 /
#' \mathrm{Var}(\mu_i)}) and \eqn{\tilde r_g} the whitened working residual
#' (\eqn{\tilde r_i = \mathrm{sign}(\mathrm{d}\mu/\mathrm{d}\eta_i)
#' \sqrt{w_i/\mathrm{Var}(\mu_i)}\,(y_i - \mu_i)}). For a Gaussian-identity fit
#' \eqn{\tilde X_g = X_g}, \eqn{\tilde r_g = y_g - X_g\hat\theta} and the downdate
#' is the exact deleted-cluster penalized solve; for other families it is the
#' one-step approximation at the converged \eqn{(W, \lambda)}. `A^{-1}` uses the
#' MODEL-BASED penalized bread (`pffr_canonicalize_cov()$model$Vp`), never the
#' robust `$pffr$Vsandwich` slot, so the jackknife bread stays genuinely
#' model-based even on a sandwich fit.
#'
#' The mean-centered CV3 jackknife variance at evaluation point \eqn{j} is
#' \eqn{V_{jack}(j) = \frac{G-1}{G}\sum_g (\eta_{(-g),j} -
#' \bar\eta_{\cdot,j})^2}. Eigenvalues of \eqn{I - H_{gg}} are floored at
#' `eig_floor` before inversion; the number of floored clusters and the maximum
#' per-cluster leverage are returned for diagnostics.
#'
#' @param object A fitted [pffr()] model (single linear-predictor family).
#' @param newdata Optional prediction data (as in [predict.pffr()]); `NULL`
#'   evaluates at the fitted observation points.
#' @param cluster Optional per-curve grouping (one entry per curve) forcing a
#'   coarser leave-one-cluster-out level, as in [pffr_vcov()]. `NULL` clusters
#'   by curve.
#' @param eig_floor Eigenvalue floor for \eqn{(I - H_{gg})^{-1}} (default
#'   `1e-8`, matching the X3 prototype's `X3_JACK_EIG_FLOOR`).
#' @param smw_check If `TRUE`, additionally verify the SMW downdate against a
#'   direct solve of the deleted-cluster (whitened) penalized normal equations
#'   for the first cluster, returning the max abs coefficient discrepancy in
#'   `smw_max_abs_err`.
#' @returns A list with `var_link` (per evaluation point LOCO jackknife variance
#'   of \eqn{\hat\eta}), `eta_hat` (the fitted linear predictor at the evaluation
#'   points), `G`, `max_eig_Hgg` (per cluster), `n_floored`, and
#'   `smw_max_abs_err`.
#' @keywords internal
pffr_jackknife_core <- function(
  object,
  newdata = NULL,
  cluster = NULL,
  eig_floor = 1e-8,
  smw_check = FALSE
) {
  if (!inherits(object, "pffr")) {
    stop("`object` must be a fitted pffr model.", call. = FALSE)
  }
  fam <- object$family
  if (is.null(fam$mu.eta) || is.null(fam$variance)) {
    stop(
      "pffr_jackknife_se() supports single linear-predictor families with a ",
      "standard variance()/mu.eta() (e.g. gaussian, poisson, binomial, Gamma, ",
      "scat); family '",
      as.character(fam$family)[1],
      "' (e.g. gaulss / multivariate) is not supported.",
      call. = FALSE
    )
  }

  # Cluster ids at the fitted (missing-removed) training rows. Error if there is
  # no usable cluster structure (need >= 2 independent curves/clusters).
  cluster_id <- build_cluster_id(object$pffr, cluster = cluster)
  G <- length(unique(cluster_id))
  if (!is.finite(G) || G < 2L) {
    stop(
      "The leave-one-cluster-out jackknife needs at least two independent ",
      "curves/clusters; this fit resolves to G = ",
      G,
      ". Supply a `cluster` grouping with >= 2 levels.",
      call. = FALSE
    )
  }

  # MODEL-BASED penalized bread (S1: $Vp is inviolate/model-based; use the
  # canonical accessor so old-format fits are handled too). NEVER the sandwich.
  canon <- pffr_canonicalize_cov(object)
  Vp <- canon$model$Vp
  sig2 <- object$sig2
  if (is.null(sig2) || !is.finite(sig2) || sig2 <= 0) sig2 <- 1
  A_inv <- Vp / sig2 # penalized bread (X'WX + S_lambda)^{-1}
  theta <- object$coefficients

  # Training design at the FITTED rows (align with residuals + cluster_id).
  X_full <- predict(object, type = "lpmatrix", reformat = FALSE)
  mi <- object$pffr$missing_indices
  X_train <- if (!is.null(mi)) X_full[-mi, , drop = FALSE] else X_full
  y <- as.vector(object$y)
  mu <- as.vector(object$fitted.values)
  eta <- as.vector(object$linear.predictors)
  if (
    nrow(X_train) != length(y) ||
      length(cluster_id) != length(y) ||
      length(mu) != length(y) ||
      length(eta) != length(y)
  ) {
    stop(
      "pffr_jackknife_core: internal length mismatch between the training ",
      "design, residuals and cluster ids.",
      call. = FALSE
    )
  }

  # Fisher whitening of the training rows (undispersioned; sigma^2 lives in the
  # bread A_inv = Vp/sig2). For Gaussian-identity this is the identity map, so
  # Xt = X_train and rt = y - X_train theta exactly (prototype's unwhitened path).
  pw <- object$prior.weights
  if (is.null(pw)) pw <- rep(1, length(y))
  pw <- as.vector(pw)
  mu_eta <- as.vector(fam$mu.eta(eta))
  var_mu <- as.vector(fam$variance(mu))
  w_star <- pw * mu_eta^2 / var_mu
  s <- sqrt(w_star)
  s[!is.finite(s)] <- 0
  Xt <- X_train * s
  Xt[!is.finite(Xt)] <- 0
  rt <- sign(mu_eta) * sqrt(pw / var_mu) * (y - mu)
  rt[!is.finite(rt)] <- 0

  # Evaluation design: fitted points (default) or newdata.
  X_eval <- if (is.null(newdata)) {
    X_train
  } else {
    predict(object, newdata = newdata, type = "lpmatrix", reformat = FALSE)
  }
  eta_hat_eval <- as.vector(X_eval %*% theta)
  n_eval <- nrow(X_eval)

  groups <- unique(cluster_id)
  # D[g, ] = X_eval %*% delta_g (the LOO shift of eta at each evaluation point).
  Dmat <- matrix(0, nrow = G, ncol = n_eval)
  max_eig_Hgg <- numeric(G)
  n_floored <- 0L
  smw_max_abs_err <- NA_real_
  A <- if (isTRUE(smw_check)) solve(A_inv) else NULL

  for (gi in seq_along(groups)) {
    idx <- which(cluster_id == groups[gi])
    Xtg <- Xt[idx, , drop = FALSE]
    rtg <- rt[idx]
    Hgg <- Xtg %*% A_inv %*% t(Xtg)
    Hgg <- 0.5 * (Hgg + t(Hgg))

    ee <- eigen(diag(length(idx)) - Hgg, symmetric = TRUE)
    max_eig_Hgg[gi] <- 1 - min(ee$values)
    vals <- ee$values
    if (any(vals < eig_floor)) {
      n_floored <- n_floored + 1L
      vals <- pmax(vals, eig_floor)
    }
    # (I - H_gg)^{-1} rt_g via the (floored) eigendecomposition.
    u <- ee$vectors %*% (crossprod(ee$vectors, rtg) / vals)
    delta <- A_inv %*% crossprod(Xtg, u)
    Dmat[gi, ] <- as.vector(X_eval %*% delta)

    if (isTRUE(smw_check) && gi == 1L) {
      # Direct solve of the deleted-cluster whitened penalized normal equations:
      # (A - Xtg'Xtg) theta_(-g) = A theta - Xtg'Xtg theta - Xtg' rt_g, with
      # A = (Vp/sig2)^{-1} = Xt'Xt + S_lambda. For Gaussian-identity this equals
      # the prototype's (A - Xg'Xg) theta_(-g) = A theta - Xg' y_g.
      rhs <- as.vector(A %*% theta) -
        as.vector(crossprod(Xtg, Xtg %*% theta)) -
        as.vector(crossprod(Xtg, rtg))
      theta_direct <- solve(A - crossprod(Xtg), rhs)
      theta_smw <- theta - as.vector(delta)
      smw_max_abs_err <- max(abs(theta_direct - theta_smw))
    }
  }

  eta_bar <- colMeans(Dmat)
  dev <- sweep(Dmat, 2L, eta_bar)
  var_link <- (G - 1) / G * colSums(dev^2)

  list(
    var_link = var_link,
    eta_hat = eta_hat_eval,
    G = G,
    max_eig_Hgg = max_eig_Hgg,
    n_floored = n_floored,
    smw_max_abs_err = smw_max_abs_err
  )
}

#' Leave-one-cluster-out jackknife standard errors for the fitted mean
#'
#' Per-evaluation-point standard errors (and Wald intervals) for the fitted
#' linear predictor / response mean of a [pffr()] fit, from an exact
#' leave-one-cluster-out (LOCO) jackknife of the penalized normal equations. It
#' is the recommended small-\eqn{G} path for fitted-mean / response-scale
#' (\eqn{E(Y)}) intervals on cluster-robust fits, where the plug-in
#' cluster/CL2 sandwich under-propagates the aggregated variance because the
#' estimated cross-term blocks of the rank-\eqn{\le G} meat inject noise that
#' shrinks the aggregated \eqn{E(Y)} SE (see \sQuote{References}).
#'
#' The jackknife recomputes the deleted-cluster coefficients by an exact
#' Sherman--Morrison--Woodbury downdate of the penalized (weighted) normal
#' equations at fixed smoothing parameters \eqn{\lambda} and converged working
#' weights, using the MODEL-BASED penalized bread \eqn{A^{-1} = V_p/\sigma^2}
#' (the fit's inviolate `$Vp`, never the robust `$pffr$Vsandwich`). See
#' [pffr_jackknife_core()] for the algorithm.
#'
#' @section Calibration caveat (read this):
#' This interval reaches only \eqn{\approx 0.86}--\eqn{0.90} pointwise coverage
#' --- \emph{not} nominal --- at the hardest simulated cells (small \eqn{G},
#' e.g. \eqn{G = 20}, with strong within-curve AR(1) dependence): there it
#' over-inflates some replicates and compensates, so it is a workable interval,
#' not a calibrated pivot. It is \emph{exact} only for Gaussian-identity
#' responses at fixed smoothing parameters; for other families it is a
#' \emph{one-step approximation} at the converged working weights and fixed
#' \eqn{\lambda} (measured max coefficient gap versus a direct fixed-\eqn{sp}
#' deleted-cluster IRLS refit on the validation Poisson fit: about
#' \eqn{2\times 10^{-3}}, i.e. \eqn{\approx}0.2\% of the coefficient scale; see
#' the package tests). The default reference is \eqn{t_{G-1}}.
#'
#' @param object A fitted [pffr()] model with a single linear predictor
#'   (e.g. `family = gaussian()`, `poisson()`, `binomial()`, `Gamma()`,
#'   `scat()`). `gaulss` / multivariate families are not supported.
#' @param newdata Optional prediction data, in the format supplied to [pffr()]
#'   (as in [predict.pffr()]). `NULL` (default) evaluates at the fitted
#'   observation points.
#' @param alpha Significance level for the Wald intervals; default `0.05` for
#'   95\% intervals.
#' @param crit Critical-value reference: `"tG1"` (default, the \eqn{t_{G-1}}
#'   reference recommended at small \eqn{G}) or `"z"` (Gaussian).
#' @param se_scale `"link"` (default) returns the SE of the linear predictor;
#'   `"response"` applies the delta method via `family$mu.eta()` and returns the
#'   fitted mean, SE and interval on the response scale.
#' @param eig_floor Eigenvalue floor for the \eqn{(I - H_{gg})^{-1}} inversion
#'   in the SMW downdate (default `1e-8`).
#' @param cluster Optional per-curve grouping (one entry per curve) forcing a
#'   coarser leave-one-cluster-out level; `NULL` (default) clusters by curve.
#' @returns A data frame with one row per evaluation point (in `predict`
#'   lpmatrix order --- curve-major, index-fastest) and columns `fit`, `se`,
#'   `lower`, `upper`, plus attributes `G`, `n_floored`, `max_eig_Hgg` (per
#'   cluster), `crit`, `crit_value`, `df` (`G - 1` for `"tG1"`, else `NA`),
#'   `se_scale` and `alpha`.
#' @references
#' The \eqn{E(Y)} under-propagation mechanism (cross-term rank deficit of the
#' rank-\eqn{\le G} meat) and this jackknife's simulated coverage are documented
#' in the paper's \eqn{E(Y)} section (\code{sec-eymean}) and in
#' \code{notes/X3-phase2-findings.md} of the accompanying study repository:
#' the jackknife reaches 0.856 (\eqn{z}) / 0.864 (\eqn{t_{G-1}}) at AR(1),
#' \eqn{G = 20} versus 0.231 for the plug-in cluster sandwich.
#' @seealso [pffr_vcov()], [predict.pffr()] (with `se_method = "jackknife"`),
#'   [pffr_coefboot()].
#' @export
#' @author Fabian Scheipl
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- pffr_simulate(Y ~ ff(X1), n = 20, nxgrid = 15, nygrid = 15)
#' m <- pffr(Y ~ ff(X1), yind = attr(d, "yindex"), data = d,
#'           sandwich = "cluster")
#' jk <- pffr_jackknife_se(m)
#' head(jk)
#' attr(jk, "G")
#' }
pffr_jackknife_se <- function(
  object,
  newdata = NULL,
  alpha = 0.05,
  crit = c("tG1", "z"),
  se_scale = c("link", "response"),
  eig_floor = 1e-8,
  cluster = NULL
) {
  crit <- match.arg(crit)
  se_scale <- match.arg(se_scale)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single value in (0, 1).", call. = FALSE)
  }

  core <- pffr_jackknife_core(
    object,
    newdata = newdata,
    cluster = cluster,
    eig_floor = eig_floor
  )
  G <- core$G
  se_link <- sqrt(pmax(core$var_link, 0))
  eta_hat <- core$eta_hat

  crit_df <- if (crit == "tG1") G - 1 else NA_real_
  crit_value <- if (crit == "tG1") {
    qt(1 - alpha / 2, df = crit_df)
  } else {
    qnorm(1 - alpha / 2)
  }

  if (se_scale == "response") {
    mu_eta_eval <- as.vector(object$family$mu.eta(eta_hat))
    se <- se_link * abs(mu_eta_eval)
    fit <- as.vector(object$family$linkinv(eta_hat))
  } else {
    se <- se_link
    fit <- eta_hat
  }

  out <- data.frame(
    fit = fit,
    se = se,
    lower = fit - crit_value * se,
    upper = fit + crit_value * se
  )
  attr(out, "G") <- G
  attr(out, "n_floored") <- core$n_floored
  attr(out, "max_eig_Hgg") <- core$max_eig_Hgg
  attr(out, "crit") <- crit
  attr(out, "crit_value") <- crit_value
  attr(out, "df") <- crit_df
  attr(out, "se_scale") <- se_scale
  attr(out, "alpha") <- alpha
  out
}

#' Upgrade an old-format pffr fit to the current covariance storage contract
#'
#' Older refund versions (storage format 1) overwrote a sandwich-corrected
#' fit's `$Vp`/`$Vc`/`$Ve` with the robust covariance and stashed the
#' model-based matrices in `$pffr$model_cov`. This converts such a fit to the
#' current contract: `$Vp`/`$Vc`/`$Ve` become model-based again and the robust
#' covariance moves to `$pffr$Vsandwich`. Current-format fits and
#' `sandwich = "none"` fits are returned unchanged (with a message). Fits that
#' lack the model-based stash cannot be upgraded (they must be refitted).
#'
#' @param object A fitted pffr model.
#' @returns The fit in the current storage format.
#' @export
#' @seealso [pffr()]
pffr_upgrade_fit <- function(object) {
  if (!inherits(object, "pffr")) {
    stop("`object` must be a fitted pffr model.", call. = FALSE)
  }
  meta <- object$pffr
  if (
    isTRUE((meta$cov_format %||% 0L) >= 2L) ||
      !is.null(meta$sandwich_info) ||
      !is.null(meta$Vsandwich)
  ) {
    message("pffr fit is already in the current storage format; nothing to do.")
    return(object)
  }
  fit_type <- normalize_sandwich_type(meta$sandwich)
  if (identical(fit_type, "none")) {
    message("pffr fit uses sandwich = \"none\"; nothing to upgrade.")
    return(object)
  }
  mc <- meta$model_cov
  if (is.null(mc)) {
    warning(
      "Cannot upgrade: this fit does not carry the model-based covariance ",
      "(created by a very old refund version). Refit with the current version.",
      call. = FALSE
    )
    return(object)
  }
  robust <- object$Vp
  robust_freq <- object$Ve
  object$Vp <- mc$Vp
  object$Vc <- mc$Vc
  object$Ve <- mc$Ve
  object$pffr$Vsandwich <- robust
  object$pffr$Vsandwich_freq <- robust_freq
  object$pffr$sandwich_info <- list(
    type = fit_type,
    cluster_var = NULL,
    G = NA_integer_,
    n_capped = meta$cl2_n_capped %||% 0L,
    max_leverage = NA_real_,
    dof_correction = meta$dof_correction %||% "none",
    edf_type = meta$edf_type %||% "trace",
    version = as.character(utils::packageVersion("refund")),
    storage_format = PFFR_COV_STORAGE_FORMAT
  )
  object$pffr$Vsandwich_cache <- new.env(parent = emptyenv())
  object$pffr$model_cov <- NULL
  object$pffr$cov_format <- PFFR_COV_STORAGE_FORMAT
  message(
    "Upgraded pffr fit to the current storage format: $Vp/$Vc/$Ve are now ",
    "model-based; the robust covariance is in $pffr$Vsandwich."
  )
  object
}

#' Store a robust sandwich covariance on a fitted pffr model
#'
#' Computes the requested robust covariance (observation-level HC via
#' [mgcv::vcov.gam()], or cluster-robust CR1 via [gam_sandwich_cluster()] / CL2
#' via [gam_sandwich_cluster_cl2()]) from the fit's model-based bread and stores
#' it in `$pffr$Vsandwich` (+ `$pffr$Vsandwich_freq`) together with
#' `$pffr$sandwich_info`. The model-based `$Vp`/`$Vc`/`$Ve` are left untouched,
#' so recomputing a sandwich later never double-applies the correction.
#'
#' @param m Fitted model.
#' @param algorithm Algorithm symbol.
#' @param type `"cluster"` (default) for CR1 cluster-robust, `"cl2"` for
#'   leverage-adjusted cluster-robust CL2, or `"hc"` for observation-level HC.
#' @param dof_correction Optional CR1 small-sample dof correction
#'   (`"none"`/`"edf"`); applied only when `type = "cluster"`. Ignored (with a
#'   warning) for `"cl2"`, which already corrects per-cluster leverage.
#' @param edf_type Which EDF the `"edf"` correction uses (`"trace"`/`"edf2"`/
#'   `"basis"`).
#' @returns Model with the robust covariance stored in `$pffr$Vsandwich`.
#' @keywords internal
apply_sandwich_correction <- function(
  m,
  algorithm,
  type = "cluster",
  dof_correction = "none",
  edf_type = "trace"
) {
  gam_obj <- if (as.character(algorithm) %in% c("gamm4", "gamm")) m$gam else m

  if (type == "cl2" && !identical(dof_correction, "none")) {
    warning(
      "dof_correction = \"",
      dof_correction,
      "\" is ignored for sandwich = \"cl2\": the CL2 leverage adjustment ",
      "already targets the same small-sample bias.",
      call. = FALSE
    )
  }

  # $Vp/$Vc/$Ve stay model-based ALWAYS. The sandwich estimators read them as
  # the penalized bread; we store the resulting robust matrices separately.
  bread <- pffr_model_based_gam(gam_obj)
  cluster_id <- if (type %in% c("cluster", "cl2")) {
    build_cluster_id(gam_obj$pffr)
  } else {
    NULL
  }

  Vsw <- pffr_compute_sandwich(
    bread,
    type,
    cluster_id,
    freq = FALSE,
    dof_correction = dof_correction,
    edf_type = edf_type
  )
  Vsw_freq <- pffr_compute_sandwich(
    bread,
    type,
    cluster_id,
    freq = TRUE,
    dof_correction = dof_correction,
    edf_type = edf_type
  )

  n_capped <- attr(Vsw, "n_capped_clusters") %||% 0L
  max_lev <- attr(Vsw, "max_leverage") %||% NA_real_

  gam_obj$pffr$Vsandwich <- Vsw
  gam_obj$pffr$Vsandwich_freq <- Vsw_freq
  gam_obj$pffr$sandwich_info <- list(
    type = type,
    cluster_var = NULL,
    G = if (!is.null(cluster_id)) length(unique(cluster_id)) else NA_integer_,
    n_capped = n_capped,
    max_leverage = max_lev,
    dof_correction = if (type == "cluster") dof_correction else "none",
    edf_type = edf_type,
    version = as.character(utils::packageVersion("refund")),
    storage_format = PFFR_COV_STORAGE_FORMAT
  )
  # Keep the legacy CL2 leverage-cap diagnostic slot populated.
  gam_obj$pffr$cl2_n_capped <- if (type == "cl2") n_capped else NULL
  # Fresh cache for on-demand recomputation of other sandwich types
  # (fit$pffr$Vsandwich_cache[[type]]).
  gam_obj$pffr$Vsandwich_cache <- new.env(parent = emptyenv())

  if (as.character(algorithm) %in% c("gamm4", "gamm")) {
    m$gam <- gam_obj
  } else {
    m <- gam_obj
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
#' @param sandwich Character sandwich type (`"none"`, `"cluster"`, `"cl2"`,
#'   or `"hc"`).
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
  sandwich,
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
    sandwich = sandwich,
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
