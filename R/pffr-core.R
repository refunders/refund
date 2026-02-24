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
#' @param check_ar Whether to check AR-related arguments (TRUE for pffr, FALSE for pffr_gls).
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
    sandwich = sandwich
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
# Step 8: Formula Term Processing (shared by pffr and pffr_gls)
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
#' @returns Integer vector of length equal to the number of fitted rows.
#' @keywords internal
build_cluster_id <- function(pffr_meta) {
  if (isTRUE(pffr_meta$is_sparse)) {
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
#' Uses a pseudo-observation expansion with one pseudo-row per linear predictor
#' component (location and scale).
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

  # Same score components used by compute_gaulss_scores()
  w1 <- (tau^2 * r) * b$family$linfo[[1]]$mu.eta(eta1)
  w2 <- (1 / tau - tau * r^2) * b$family$linfo[[2]]$mu.eta(eta2)

  pw <- b$prior.weights
  if (!is.null(pw) && any(pw != 1)) {
    w1 <- pw * w1
    w2 <- pw * w2
  }

  make_block <- function(w, lp_idx) {
    abs_sqrt <- sqrt(abs(w))
    abs_sqrt[!is.finite(abs_sqrt)] <- 0

    z <- sign(w) * abs_sqrt
    z[!is.finite(z)] <- 0

    Xw_block <- matrix(0, nrow = n_obs, ncol = ncol(X))
    if (!is.null(lp_idx) && length(lp_idx) > 0) {
      Xw_block[, lp_idx] <- X[, lp_idx, drop = FALSE] * abs_sqrt
      Xw_block[!is.finite(Xw_block)] <- 0
    }
    list(Xw = Xw_block, z = z)
  }

  block1 <- make_block(w1, lpi[[1]])
  block2 <- make_block(w2, lpi[[2]])

  list(
    Xw = rbind(block1$Xw, block2$Xw),
    z = c(block1$z, block2$z),
    cluster_id = rep(cluster_id, times = 2L)
  )
}

#' Assemble cluster-robust sandwich from score matrix
#'
#' Given a per-observation score matrix, aggregates by cluster and forms
#' \eqn{V_{CL} = c \cdot V_p M_{CL} V_p + B_2} with HC1 correction.
#'
#' @param scores Per-observation score matrix (n_obs x p).
#' @param cluster_id Cluster membership vector.
#' @param Vp Bayesian posterior covariance (p x p).
#' @param B2 Bias correction matrix (p x p, or scalar 0).
#' @returns A p x p covariance matrix.
#' @keywords internal
assemble_cluster_sandwich <- function(scores, cluster_id, Vp, B2) {
  U <- rowsum(scores, cluster_id)
  meat <- crossprod(U)
  hc1 <- length(unique(cluster_id)) / (length(unique(cluster_id)) - 1)
  hc1 * Vp %*% meat %*% Vp + B2
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
#' @returns A p x p covariance matrix.
#' @keywords internal
gam_sandwich_cluster <- function(b, cluster_id, freq = FALSE) {
  B2 <- if (freq) 0 else b$Vp - b$Ve
  X <- model.matrix(b)

  if (b$family$family == "gaulss") {
    scores <- compute_gaulss_scores(b, X)
    return(assemble_cluster_sandwich(scores, cluster_id, b$Vp, B2))
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
  # w = (d mu/d eta) * (y - mu) / (phi * V(mu))
  mu <- b$fitted.values
  scores <- b$family$mu.eta(b$linear.predictors) *
    (b$y - mu) /
    (b$sig2 * b$family$variance(mu)) *
    X
  assemble_cluster_sandwich(scores, cluster_id, b$Vp, B2)
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
#' @returns A p x p covariance matrix.
#' @keywords internal
gam_sandwich_cluster_cl2 <- function(
  b,
  cluster_id,
  freq = FALSE,
  tol = 1e-8,
  leverage_cap = 0.999
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

  groups <- unique(cluster_id_work)
  G <- length(groups)
  if (G < 2) {
    stop("Need at least two clusters for cluster sandwich.", call. = FALSE)
  }

  Vp <- b$Vp
  B2 <- if (freq) 0 else b$Vp - b$Ve
  p <- ncol(Xw)
  meat <- matrix(0, nrow = p, ncol = p)
  n_capped_clusters <- 0L

  for (g in groups) {
    idx <- which(cluster_id_work == g)
    Xwg <- Xw[idx, , drop = FALSE]
    zg <- z[idx]

    Hgg <- Xwg %*% Vp %*% t(Xwg)
    Hgg <- 0.5 * (Hgg + t(Hgg))

    ee_H <- eigen(Hgg, symmetric = TRUE)
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
  }

  hc1 <- G / (G - 1)
  V <- hc1 * Vp %*% meat %*% Vp + B2
  V <- 0.5 * (V + t(V))
  attr(V, "n_capped_clusters") <- n_capped_clusters
  V
}

#' Apply sandwich correction to model covariance matrices
#'
#' Applies either observation-level HC sandwich (via [mgcv::vcov.gam()]) or
#' cluster-robust sandwich (CR1 via [gam_sandwich_cluster()] or CL2 via
#' [gam_sandwich_cluster_cl2()]) to a fitted pffr model.
#'
#' @param m Fitted model.
#' @param algorithm Algorithm symbol.
#' @param type `"cluster"` (default) for CR1 cluster-robust, `"cl2"` for
#'   leverage-adjusted cluster-robust CL2, or `"hc"` for observation-level HC.
#' @returns Model with corrected covariance matrices.
#' @keywords internal
apply_sandwich_correction <- function(m, algorithm, type = "cluster") {
  gam_obj <- if (as.character(algorithm) %in% c("gamm4", "gamm")) m$gam else m

  gam_obj_stripped <- gam_obj
  class(gam_obj_stripped) <- setdiff(class(gam_obj_stripped), "pffr")

  if (type %in% c("cluster", "cl2")) {
    cluster_id <- build_cluster_id(gam_obj$pffr)
    cov_fun <- if (type == "cl2") {
      gam_sandwich_cluster_cl2
    } else {
      gam_sandwich_cluster
    }
    gam_obj$Vp <- gam_obj$Vc <- cov_fun(
      gam_obj_stripped,
      cluster_id,
      freq = FALSE
    )
    gam_obj$Ve <- cov_fun(
      gam_obj_stripped,
      cluster_id,
      freq = TRUE
    )
  } else if (type == "hc") {
    gam_obj$Vp <- gam_obj$Vc <- mgcv::vcov.gam(
      gam_obj_stripped,
      sandwich = TRUE,
      freq = FALSE
    )
    gam_obj$Ve <- mgcv::vcov.gam(
      gam_obj_stripped,
      sandwich = TRUE,
      freq = TRUE
    )
  } else {
    stop("Unknown sandwich type: ", type, call. = FALSE)
  }

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
