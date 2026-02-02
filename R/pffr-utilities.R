# Utility functions for pffr:

safeDeparse <- function(expr) {
  # turn an expression into a _single_ string, regardless of the expression's length
  ret <- paste(deparse(expr), collapse = "")
  #rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}

# Null coalescing operator (returns y if x is NULL)
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Return call with all possible arguments
#'
#' Return a call in which all of the arguments which were supplied or have presets are specified by their full names and their supplied or default values.
#'
#' @param definition a function. See \code{\link[base]{match.call}}.
#' @param call an unevaluated call to the function specified by definition. See \code{\link[base]{match.call}}.
#' @param expand.dots logical. Should arguments matching ... in the call be included or left as a ... argument? See \code{\link[base]{match.call}}.
#' @return An object of mode "\code{\link[base]{call}}".
#' @author Fabian Scheipl
#' @seealso \code{\link[base]{match.call}}
expand.call <- function(
  definition = NULL,
  call = sys.call(sys.parent(1)),
  expand.dots = TRUE
) {
  call <- match.call(definition, call, expand.dots)
  #given args:
  ans <- as.list(call)
  this_function <- safeDeparse(ans[[1]])
  # rm namespace operators so that function finding works (mlr-org/mlr/#1559)
  if (grepl(":", this_function)) {
    this_function <- gsub("(^[^:]+:+)(.+$)", "\\2", this_function)
    if (!exists(this_function)) {
      warning(
        "expand.call couldn't find ",
        this_function,
        " and returned unchanged call:\n <",
        call,
        ">"
      )
      return(call)
    }
  }
  frmls <- formals(this_function)
  #remove formal args with no presets:
  frmls <- frmls[!sapply(frmls, is.symbol)]
  add <- which(!(names(frmls) %in% names(ans)))
  return(as.call(c(ans, frmls[add])))
}


list2df <- function(l) {
  # make a list into a dataframe -- matrices are left as matrices!
  nrows <- sapply(l, \(x) nrow(as.matrix(x)))
  stopifnot(length(unique(nrows)) == 1)
  ret <- data.frame(rep(NA, nrows[1]))
  for (i in 1:length(l)) ret[[i]] <- l[[i]]
  names(ret) <- names(l)
  return(ret)
}


#--------------------------------------
# Matrix Window Operations for ff() Terms with Limits
#--------------------------------------

#' Extract row-wise windows from a matrix
#'
#' For each row of a matrix, extracts elements at indices specified by the
#' corresponding row of `windows`. Used for `ff()` terms with integration limits
#' to reduce matrix size by keeping only the relevant integration region.
#'
#' @param X Matrix to extract windows from.
#' @param windows Integer matrix where each row specifies column indices to
#'   extract from the corresponding row of `X`.
#' @returns Matrix with same number of rows as `X` and same number of columns
#'   as columns in `windows`.
#' @keywords internal
shift_and_shorten_matrix <- function(X, windows) {
  t(vapply(
    seq_len(nrow(X)),
    function(i) X[i, windows[i, ]],
    FUN.VALUE = numeric(ncol(windows))
  ))
}


#' Compute integration window bounds from a boolean mask
#'
#' For `ff()` terms with limits (e.g., `limits = "s<t"`), computes the column
#' range containing TRUE values for each row of the mask.
#'
#' @param use_mask Logical matrix indicating which (s,t) pairs are within the
#'   integration limits.
#' @returns Integer matrix with 3 columns: start index, end index, and width
#'   for each row.
#' @keywords internal
compute_integration_windows <- function(use_mask) {
  windows <- t(apply(use_mask, 1, function(row) {
    indices <- which(row)
    # Edge case: no integration region (all FALSE)
    if (length(indices) == 0) return(c(1L, 1L))
    range(indices)
  }))
  # Add width column
  cbind(windows, windows[, 2] - windows[, 1] + 1L)
}


#' Expand windows to uniform width
#'
#' Expands each window to have the same width (the maximum observed width),
#' extending in whichever direction has room. This enables efficient matrix
#' operations on uniformly-sized windows.
#'
#' @param windows Integer matrix from [compute_integration_windows()] with
#'   columns: start, end, width.
#' @param max_col Maximum valid column index (typically `ncol(smat)`).
#' @returns Integer matrix where each row contains the column indices for the
#'   expanded window.
#' @keywords internal
expand_windows_to_maxwidth <- function(windows, max_col) {
  max_width <- max(windows[, 3])
  t(apply(windows, 1, function(window) {
    width <- window[3]
    # Extend toward end if possible, otherwise toward start
    if ((window[2] + max_width - width) <= max_col) {
      window[1]:(window[2] + max_width - width)
    } else {
      (window[1] + width - max_width):window[2]
    }
  }))
}


#--------------------------------------
# Shared Utilities for ff/sff Terms
#--------------------------------------

#' Compute integration weights for ff/sff terms
#'
#' Computes numerical integration weights using Simpson's rule, trapezoidal
#' rule, or Riemann sums. Used by [ff()] and [sff()] for constructing
#' functional regression terms.
#'
#' @param xind Matrix of x-index values (n x nxgrid). Each row contains the
#'   evaluation points for one observation.
#' @param integration Character: "simpson", "trapezoidal", or "riemann".
#' @returns Matrix of integration weights (n x nxgrid).
#' @keywords internal
compute_integration_weights <- function(xind, integration = "simpson") {
  n <- nrow(xind)
  nxgrid <- ncol(xind)

  switch(
    integration,
    simpson = {
      # Simpson's rule: int^b_a f(t) dt ≈ (b-a)/(3n) * [f(a) + 4f(t_1) + 2f(t_2) + ...]
      ((xind[, nxgrid] - xind[, 1]) / nxgrid) /
        3 *
        matrix(
          c(1, rep(c(4, 2), length.out = nxgrid - 2), 1),
          nrow = n,
          ncol = nxgrid,
          byrow = TRUE
        )
    },
    trapezoidal = {
      # Trapezoidal rule
      diffs <- t(apply(xind, 1, diff))
      0.5 *
        cbind(
          diffs[, 1],
          t(apply(diffs, 1, filter, filter = c(1, 1)))[, -(nxgrid - 1)],
          diffs[, nxgrid - 1]
        )
    },
    riemann = {
      # Simple Riemann sums
      diffs <- t(apply(xind, 1, diff))
      cbind(rep(mean(diffs), n), diffs)
    }
  )
}


#' Validate and expand xind to matrix form
#'
#' Validates that xind has the correct dimensions and is in increasing order,
#' then expands a vector or single-row matrix to a full n x nxgrid matrix.
#' Used by [ff()] and [sff()].
#'
#' @param xind Vector or matrix of x-index values.
#' @param n Number of observations.
#' @param nxgrid Number of grid points (columns expected).
#' @param arg_name Name of the argument for error messages.
#' @returns Matrix of x-index values (n x nxgrid).
#' @keywords internal
validate_and_expand_xind <- function(xind, n, nxgrid, arg_name = "xind") {
  # Convert to matrix if needed

  if (is.null(dim(xind))) {
    xind <- t(as.matrix(xind))
  }

  # Check column count
  if (ncol(xind) != nxgrid) {
    stop(
      "`",
      arg_name,
      "` must have ",
      nxgrid,
      " columns, not ",
      ncol(xind),
      "."
    )
  }

  # Expand single-row matrix to full matrix
  if (nrow(xind) == 1) {
    xind <- matrix(as.vector(xind), nrow = n, ncol = nxgrid, byrow = TRUE)
  } else if (nrow(xind) != n) {
    stop(
      "`",
      arg_name,
      "` must be a vector, single-row matrix, or have ",
      n,
      " rows."
    )
  }

  # Check ordering
  if (!identical(order(xind[1, ]), seq_len(nxgrid))) {
    stop("`", arg_name, "` values must be in increasing order.")
  }

  xind
}


#' Convert limits argument to a function
#'
#' Parses the `limits` argument for [ff()] and [sff()] terms. Accepts
#' `NULL` (no limits), a function, or string shortcuts "s<t" / "s<=t".
#'
#' @param limits Either `NULL`, a function(s, t) returning logical, or
#'   a string "s<t" or "s<=t".
#' @returns A function(s, t) returning logical, or `NULL`.
#' @keywords internal
build_limits_function <- function(limits) {
  if (is.null(limits) || is.function(limits)) {
    return(limits)
  }

  limits_fun <- switch(
    limits,
    "s<t" = \(s, t) s < t,
    "s<=t" = \(s, t) (s < t) | (s == t),
    NULL
  )

  if (is.null(limits_fun)) {
    stop("Unknown `limits` value: ", limits, ".")
  }

  limits_fun
}


#' Create short labels for pffr terms at model fit time
#'
#' Creates a named vector mapping mgcv smooth labels to human-readable pffr labels.
#' This is called during `pffr()` fitting and stored in `object$pffr$shortlabels`.
#'
#' @param label_map Named list mapping original term names to mgcv smooth labels.
#' @param m_smooth List of smooth objects from the fitted model.
#' @param yind_name Name of the y-index variable.
#' @param where_specials List indicating which terms are which type (par, ff, ffpc, etc.).
#' @param family The model family object (used to detect location-scale families).
#' @returns Named character vector: names are mgcv smooth labels, values are short labels.
#' @keywords internal
create_shortlabels <- function(
  label_map,
  m_smooth,
  yind_name,
  where_specials,
  family = NULL
) {
  # Build result: one short label per smooth in m_smooth

  shortlabels <- character()

  for (term_name in names(label_map)) {
    mgcv_labels <- label_map[[term_name]]

    # Skip if no labels (shouldn't happen)
    if (length(mgcv_labels) == 0 || all(is.na(mgcv_labels))) next

    # Check term type based on name patterns and where_specials
    term_idx <- match(term_name, names(label_map))
    is_par <- term_idx %in% where_specials$par
    is_ffpc <- grepl("^ffpc", term_name)
    is_ff <- term_idx %in% where_specials$ff
    is_sff <- term_idx %in% where_specials$sff

    for (i in seq_along(mgcv_labels)) {
      mgcv_label <- mgcv_labels[i]
      if (is.na(mgcv_label)) next

      sm <- NULL
      if (!is.null(names(m_smooth)) && mgcv_label %in% names(m_smooth)) {
        sm <- m_smooth[[mgcv_label]]
      } else {
        sm_idx <- which(vapply(
          m_smooth,
          \(x) identical(x$label, mgcv_label),
          logical(1)
        ))
        if (length(sm_idx) > 0) {
          sm <- m_smooth[[sm_idx[1]]]
        }
      }

      # Generate short label based on term type
      if (is_par && length(mgcv_labels) > 1) {
        # Factor varying coefficient: extract factor level from smooth
        if (!is.null(sm) && !is.null(sm$by.level)) {
          short <- paste0(sm$by, sm$by.level, "(", yind_name, ")")
        } else {
          short <- paste0(term_name, "(", yind_name, ")")
        }
      } else if (is_par) {
        # Scalar varying coefficient
        short <- paste0(term_name, "(", yind_name, ")")
      } else if (is_ffpc) {
        # ffpc: extract PC number and variable name
        # mgcv_label looks like "s(yindex.vec):X.PC1"
        pc_match <- regmatches(mgcv_label, regexpr("PC[0-9]+", mgcv_label))
        if (length(pc_match) > 0) {
          pc_num <- gsub("PC", "", pc_match)
          # Extract variable from term_name (ffpc(X, ...))
          var_match <- regmatches(term_name, regexpr("\\(([^,)]+)", term_name))
          var_name <- gsub("[\\(\\)]", "", var_match)
          short <- paste0("ffpc(", var_name, ")", pc_num)
        } else {
          short <- simplify_term_label(term_name, yind_name)
        }
      } else if (is_ff || is_sff) {
        # ff/sff: use simplified version of original term
        short <- simplify_term_label(term_name, yind_name)
      } else {
        # Other terms (s, te, t2, c-wrapped, intercept)
        short <- simplify_term_label(term_name, yind_name)
      }

      shortlabels[mgcv_label] <- short
    }
  }

  # Handle location-scale family smooths (e.g., gaulss)
  # These have labels like "s.1(...)" for scale parameter and aren't in label_map
  # Only do this for families with multiple linear predictors
  is_location_scale <- !is.null(family) &&
    !is.null(family$nlp) &&
    family$nlp > 1

  if (is_location_scale) {
    all_smooth_labels <- vapply(m_smooth, \(x) x$label, character(1))
    unlabeled <- setdiff(all_smooth_labels, names(shortlabels))

    for (mgcv_label in unlabeled) {
      # Check for scale parameter pattern: smooth.N(...) where N > 0
      # Only match known mgcv smooth types: s, te, ti, t2
      # e.g., "s.1(t_grid.vec)", "te.1(x,t)", "s.1(t_grid.vec):xlin"
      scale_match <- regexpr("^(s|te|ti|t2)\\.([0-9]+)\\(", mgcv_label)
      if (scale_match > 0) {
        # Extract the base term by removing the .N suffix
        base_label <- sub("^(s|te|ti|t2)\\.([0-9]+)", "\\1", mgcv_label)

        # Check if this is just the intercept (yind_name only) or has covariates
        has_by <- grepl(":", mgcv_label)

        if (has_by) {
          # Complex case: scale depends on covariate
          # Extract the covariate part after the colon
          by_part <- sub(".*:", "", mgcv_label)
          shortlabels[mgcv_label] <- paste0(
            "log(SD): ",
            by_part,
            "(",
            yind_name,
            ")"
          )
        } else {
          # Simple case: intercept-only scale
          shortlabels[mgcv_label] <- paste0("log(SD)(", yind_name, ")")
        }
      }
    }
  }

  shortlabels
}

#' Simplify a term label for display
#'
#' @param term_name Original term name string.
#' @param yind_name Name of the y-index variable.
#' @returns Simplified label string.
#' @keywords internal
simplify_term_label <- function(term_name, yind_name) {
  # Handle special cases
  if (grepl("^Intercept\\(", term_name)) {
    return(term_name)
  }

  # Try to parse the term
  expr <- tryCatch(
    parse(text = term_name)[[1]],
    error = \(e) NULL
  )

  if (is.null(expr)) {
    # Simple variable name: add (yind_name)
    return(paste0(term_name, "(", yind_name, ")"))
  }

  # For function calls, simplify by removing extra arguments
  if (is.call(expr)) {
    fn_name <- as.character(expr[[1]])

    # Extract variable names
    vars <- all.vars(expr)
    # Remove yindex if present
    vars <- vars[!grepl("yindex|vec", vars)]

    if (length(vars) == 0) {
      return(term_name)
    }

    # Build simplified label
    if (fn_name %in% c("s", "te", "ti", "t2")) {
      if (length(vars) == 1) {
        return(paste0(fn_name, "(", vars[1], ",", yind_name, ")"))
      } else {
        return(paste0(
          fn_name,
          "(",
          paste(vars, collapse = ","),
          ",",
          yind_name,
          ")"
        ))
      }
    } else if (fn_name %in% c("ff", "sff", "ffpc")) {
      return(paste0(fn_name, "(", vars[1], ")"))
    } else {
      # Default: just use the function call with main variables
      if (length(vars) == 1) {
        return(paste0(vars[1], "(", yind_name, ")"))
      } else {
        return(paste0(fn_name, "(", paste(vars, collapse = ","), ")"))
      }
    }
  }

  # Fallback: return as-is with yind_name
  paste0(term_name, "(", yind_name, ")")
}


#' P-spline constructor with modified 'shrinkage' penalty
#'
#' Construct a B-spline basis with a modified difference penalty
#' of full rank (i.e., that also penalizes low-order polynomials).
#'
#' This penalty-basis combination is useful to avoid non-identifiability issues for \code{\link{ff}} terms.
#' See 'ts' or 'cs' in \code{\link[mgcv]{smooth.terms}}
#' for similar "shrinkage penalties" for thin plate and cubic regression splines.
#' The basic idea is to replace the k-th zero eigenvalue of the original penalty by
#' \eqn{s^k \nu_m}, where \eqn{s} is the shrinkage factor (defaults to 0.1)
#' and \eqn{\nu_m} is the smallest non-zero eigenvalue. See reference for the
#' original idea, implementation follows that in the 'ts' and 'cs' constructors
#' (see \code{\link[mgcv]{smooth.terms}}).
#'
#' @param object see \code{\link[mgcv]{smooth.construct}}. The shrinkage factor can be specified via \code{object$xt$shrink}
#' @param data see \code{\link[mgcv]{smooth.construct}}.
#' @param knots see \code{\link[mgcv]{smooth.construct}}.
#'
#' @author Fabian Scheipl;  adapted from 'ts' and 'cs' constructors by S.N. Wood.
#'
#' @references Marra, G., & Wood, S. N. (2011). Practical variable selection for generalized additive models.
#' \emph{Computational Statistics & Data Analysis}, 55(7), 2372-2387.
#' @method smooth.construct pss.smooth.spec
#' @export
#' @importFrom mgcv smooth.construct.ps.smooth.spec
smooth.construct.pss.smooth.spec <- function(object, data, knots) {
  shrink <- ifelse(is.null(object$xt$shrink), 0.1, object$xt$shrink)
  stopifnot(shrink > 0, shrink < 1)

  object <- smooth.construct.ps.smooth.spec(object, data, knots)
  nk <- object$bs.dim
  difforder <- object$m[2]
  ## add shrinkage term to penalty:
  ## Modify the penalty by increasing the penalty on the
  ## unpenalized space from zero...
  es <- eigen(object$S[[1]], symmetric = TRUE)
  ## now add a penalty on the penalty null space
  es$values[(nk - difforder + 1):nk] <- es$values[nk - difforder] * shrink
  ## ... so penalty on null space is still less than that on range space.
  object$S[[1]] <- es$vectors %*% (as.numeric(es$values) * t(es$vectors))
  object$rank <- nk
  object$null.space.dim <- 0
  class(object) <- "pss.smooth"
  object
}

#' @importFrom mgcv Predict.matrix.pspline.smooth
#' @export
Predict.matrix.pss.smooth <- function(object, data) {
  Predict.matrix.pspline.smooth(object, data)
}

#' @importFrom mgcv smooth.construct.ps.smooth.spec
#' @export
smooth.construct.ps_c.smooth.spec <- function(object, data, knots) {
  object <- smooth.construct.ps.smooth.spec(object, data, knots)
  object$C <- object$xt$C1
  object
}

#' @importFrom mgcv smooth.construct.ps.smooth.spec
#' @export
smooth.construct.fame.smooth.spec <- function(object, data, knots) {
  object$bs <- "ps"
  object <- smooth.construct.ps.smooth.spec(object, data, knots)
  # anti-kernel penalty:
  object$S[[1]] <- object$xt$penK

  object$rank <- qr(object$S[[1]])$rank
  object$null.space.dim <- object$bs.dim - object$rank

  return(object)
}

# compute dim of overlap of the span of two orthonormal matrices
trace_lv <- function(A, B, tol = 1e-10) {
  ## A, B orthnormal!!
  #Rolf Larsson, Mattias Villani (2001)
  #"A distance measure between cointegration spaces"

  if (NCOL(A) == 0 | NCOL(B) == 0) {
    return(0)
  }

  if (NROW(A) != NROW(B) | NCOL(A) > NROW(A) | NCOL(B) > NROW(B)) {
    return(NA)
  }

  trace <- if (NCOL(B) <= NCOL(A)) {
    sum(diag(t(B) %*% A %*% t(A) %*% B))
  } else {
    sum(diag(t(A) %*% B %*% t(B) %*% A))
  }
  trace
}


#--------------------------------------=
# Phase 1: pffrSim Redesign - Core Infrastructure (Section 1.1)
#--------------------------------------=

#' Parse a pffr-style formula
#'
#' Extract term types and variable names from a pffr formula using
#' `terms.formula()` with specials.
#'
#' @param formula A formula object with pffr-style terms.
#' @returns A list with components:
#'   \item{response}{Character string of the response variable name.}
#'   \item{terms}{A list of parsed terms, each with `type`, `varname`, and `call`.}
#' @keywords internal
parse_pffr_formula <- function(formula) {
  # Define specials for pffr terms
  specials <- c("ff", "sff", "ffpc", "pcre", "s", "te", "ti", "t2", "c")

  # Parse formula with specials
  tf <- terms.formula(formula, specials = specials)

  # Get response name
  response <- if (attr(tf, "response") == 1L) {
    all.vars(formula)[1L]
  } else {
    NULL
  }

  # Get term labels
  term_labels <- attr(tf, "term.labels")

  # Parse each term
  terms_list <- lapply(
    term_labels,
    \(term_str) parse_single_term(term_str, specials)
  )

  # Name by term label

  names(terms_list) <- term_labels

  list(
    response = response,
    terms = terms_list
  )
}

#' Parse a single term string
#'
#' @param term_str Character string of the term.
#' @param specials Character vector of special term types.
#' @returns A list with `type`, `varname`, and `call`.
#' @keywords internal
parse_single_term <- function(term_str, specials) {
  # Try to parse as expression

  expr <- tryCatch(parse(text = term_str)[[1]], error = \(e) NULL)

  if (is.null(expr)) {
    return(list(type = "linear", varname = term_str, call = NULL))
  }

  # Check if it's a function call

  if (is.call(expr)) {
    fn_name <- as.character(expr[[1]])

    # Handle c() wrapper
    if (fn_name == "c") {
      inner <- parse_single_term(safeDeparse(expr[[2]]), specials)
      inner$wrapped <- TRUE
      return(inner)
    }

    # Check for special terms
    if (fn_name %in% specials) {
      # Extract first argument as varname
      varname <- if (length(expr) > 1) {
        arg1 <- expr[[2]]
        if (is.name(arg1)) {
          as.character(arg1)
        } else {
          safeDeparse(arg1)
        }
      } else {
        term_str
      }

      return(list(
        type = fn_name,
        varname = varname,
        call = expr,
        wrapped = FALSE
      ))
    }
  }

  # Default: linear term
  list(type = "linear", varname = term_str, call = NULL, wrapped = FALSE)
}


# ------------------------------------------------------------------------------
# P-spline Prior Sampling for Random Truth Generation
# ------------------------------------------------------------------------------

#' Sample coefficients from 1D P-spline prior
#'
#' Samples coefficients from N(0, P^{-1}) where
#' P = (1/wiggliness) * K'K + I and K is a 2nd-order difference matrix.
#'
#' @param k Number of basis functions.
#' @param wiggliness Controls smoothness. Higher = more wiggly.
#' @returns Numeric vector of coefficients (length k).
#' @keywords internal
sample_pspline_coef_1d <- function(k, wiggliness) {
  K <- diff(diag(k), differences = 2)
  P <- (1 / wiggliness) * crossprod(K) + diag(k)
  z <- rnorm(k)
  backsolve(chol(P), z)
}

#' Sample coefficients from 2D tensor P-spline prior
#'
#' Samples from tensor product P-spline prior with additive penalties:
#' P = (1/wiggliness) * (P_s ⊗ I_t + I_s ⊗ P_t) + I
#'
#' @param k_s Basis dimension for first margin.
#' @param k_t Basis dimension for second margin.
#' @param wiggliness Controls smoothness.
#' @returns Matrix of coefficients (k_s x k_t).
#' @keywords internal
sample_pspline_coef_2d <- function(k_s, k_t, wiggliness) {
  K_s <- diff(diag(k_s), differences = 2)
  K_t <- diff(diag(k_t), differences = 2)

  P_s <- crossprod(K_s)
  P_t <- crossprod(K_t)

  I_s <- diag(k_s)
  I_t <- diag(k_t)
  P_tensor <- kronecker(P_s, I_t) + kronecker(I_s, P_t)

  k_total <- k_s * k_t
  P <- (1 / wiggliness) * P_tensor + diag(k_total)

  z <- rnorm(k_total)
  coef_vec <- backsolve(chol(P), z)
  matrix(coef_vec, nrow = k_s, ncol = k_t)
}

#' Create B-spline specification for evaluation at new points
#'
#' Stores the knots and degree of a B-spline basis so it can be re-evaluated
#' at new points while maintaining the same basis functions.
#'
#' @param x Numeric vector of points to build the basis on.
#' @param k Number of basis functions.
#' @returns List with knots, boundary_knots, and degree.
#' @keywords internal
make_bs_spec <- function(x, k) {
  B <- splines::bs(x, df = k, degree = 3, intercept = TRUE)
  list(
    knots = attr(B, "knots"),
    boundary_knots = attr(B, "Boundary.knots"),
    degree = attr(B, "degree"),
    k = k
  )
}

#' Evaluate B-spline at new points using stored specification
#'
#' @param x Numeric vector of evaluation points.
#' @param spec Specification from [make_bs_spec()].
#' @returns B-spline basis matrix.
#' @keywords internal
eval_bs <- function(x, spec) {
  splines::bs(
    x,
    df = spec$k,
    knots = spec$knots,
    Boundary.knots = spec$boundary_knots,
    degree = spec$degree,
    intercept = TRUE
  )
}

#' Generate 1D random function from P-spline prior
#'
#' Samples from the P-spline prior and evaluates on a fixed grid.
#' Returns a numeric vector scaled to unit SD.
#'
#' @param grid Numeric vector of evaluation points.
#' @param k Number of B-spline basis functions.
#' @param wiggliness Controls smoothness. Higher = more wiggly.
#' @returns Numeric vector, scaled to unit SD.
#' @keywords internal
gen_random_1d <- function(grid, k, wiggliness) {
  B <- splines::bs(grid, df = k, degree = 3, intercept = TRUE)
  coef <- sample_pspline_coef_1d(k, wiggliness)
  f <- as.vector(B %*% coef)
  as.vector(scale(f))
}

#' Generate 2D random surface from tensor P-spline prior
#'
#' Samples from tensor product P-spline prior and evaluates on fixed grids.
#' Returns a matrix scaled to unit SD.
#'
#' @param s_grid First dimension evaluation points.
#' @param t_grid Second dimension evaluation points.
#' @param k_s Basis dimension for s.
#' @param k_t Basis dimension for t.
#' @param wiggliness Controls smoothness.
#' @returns Matrix (length(s_grid) x length(t_grid)), scaled to unit SD.
#' @keywords internal
gen_random_2d <- function(s_grid, t_grid, k_s, k_t, wiggliness) {
  B_s <- splines::bs(s_grid, df = k_s, degree = 3, intercept = TRUE)
  B_t <- splines::bs(t_grid, df = k_t, degree = 3, intercept = TRUE)

  coef_mat <- sample_pspline_coef_2d(k_s, k_t, wiggliness)
  f_mat <- B_s %*% coef_mat %*% t(B_t)

  f_centered <- f_mat - mean(f_mat)
  f_centered / sd(f_centered)
}

#' Create 2D random surface function from P-spline prior
#'
#' Returns a function that can be evaluated at arbitrary (s, t) points.
#' The surface is scaled to unit SD on the reference grids.
#'
#' @param s_ref Reference grid for first dimension (for scaling).
#' @param t_ref Reference grid for second dimension (for scaling).
#' @param k_s Basis dimension for s.
#' @param k_t Basis dimension for t.
#' @param wiggliness Controls smoothness.
#' @returns Function(s, t) returning a matrix.
#' @keywords internal
make_random_2d_fn <- function(s_ref, t_ref, k_s, k_t, wiggliness) {
  # Store basis specs for evaluation at arbitrary points
  s_spec <- make_bs_spec(s_ref, k_s)
  t_spec <- make_bs_spec(t_ref, k_t)

  # Sample coefficients
  coef_mat <- sample_pspline_coef_2d(k_s, k_t, wiggliness)

  # Compute scaling factor on reference grid
  B_s_ref <- eval_bs(s_ref, s_spec)
  B_t_ref <- eval_bs(t_ref, t_spec)
  f_ref <- B_s_ref %*% coef_mat %*% t(B_t_ref)
  scale_center <- mean(f_ref)
  scale_sd <- sd(f_ref)

  # Return function that evaluates at arbitrary points
  function(s, t) {
    B_s <- eval_bs(s, s_spec)
    B_t <- eval_bs(t, t_spec)
    f_mat <- B_s %*% coef_mat %*% t(B_t)
    (f_mat - scale_center) / scale_sd
  }
}

# ------------------------------------------------------------------------------
# Random Truth Generator for Reproducible Random Effects
# ------------------------------------------------------------------------------

#' Create random truth generator with P-spline prior
#'
#' Generates reproducible random truth functions by sampling from P-spline priors.
#' Uses the same penalty structure as mgcv, so wiggliness maps to effective
#' degrees of freedom. Truth is in the B-spline span, avoiding approximation bias.
#'
#' @param yind Response grid points.
#' @param xind Functional covariate grid points.
#' @param z_ref Reference z values for smooth terms.
#' @param wiggliness Controls smoothness (default: 1). Higher = more wiggly.
#'   Typical range: 0.001 (very smooth) to 10 (very wiggly).
#' @param k_truth List with basis dimensions per term type. Defaults:
#'   ff_s = 8, ff_t = 8, smooth_z = 8, smooth_t = 8,
#'   linear = 8, intercept = 8, concurrent = 8.
#' @returns List with generator functions for each term type.
#' @keywords internal
make_random_truth_generator <- function(
  yind,
  xind,
  z_ref = NULL,
  wiggliness = 1,
  k_truth = list()
) {
  # Merge with defaults
  # Use k=8 for all terms so wiggliness parameter has consistent meaning
  k <- modifyList(
    list(
      ff_s = 8L,
      ff_t = 8L,
      smooth_z = 8L,
      smooth_t = 8L,
      linear = 8L,
      intercept = 8L,
      concurrent = 8L
    ),
    k_truth
  )

  z_ref <- z_ref %||% seq(-2, 2, length.out = 20)

  # Pre-generate random functions (deterministic order for reproducibility)
  # 1D functions are evaluated on fixed yind grid
  f_intercept <- gen_random_1d(yind, k = k$intercept, wiggliness = wiggliness)
  f_linear <- gen_random_1d(yind, k = k$linear, wiggliness = wiggliness)
  f_concurrent <- gen_random_1d(yind, k = k$concurrent, wiggliness = wiggliness)

  # 2D ff is evaluated on fixed xind x yind grid
  f_ff <- gen_random_2d(
    xind,
    yind,
    k_s = k$ff_s,
    k_t = k$ff_t,
    wiggliness = wiggliness
  )

  # Smooth needs to evaluate at arbitrary z values, so return a function
  f_smooth_fn <- make_random_2d_fn(
    z_ref,
    yind,
    k_s = k$smooth_z,
    k_t = k$smooth_t,
    wiggliness = wiggliness
  )

  list(
    intercept = function(t) f_intercept,
    ff = function(s, t) f_ff,
    linear = function(t) f_linear,
    concurrent = function(t) f_concurrent,
    smooth = f_smooth_fn
  )
}

# ------------------------------------------------------------------------------
# Preset Effect Libraries (Internal, not exported)
# ------------------------------------------------------------------------------

#' Function-on-function effect presets
#'
#' Named list of preset beta(s,t) coefficient surfaces for ff() terms.
#' Each function takes arguments (s, t) and returns a matrix.
#' @keywords internal
ff_presets <- list(
  cosine = function(s, t) {
    s * cos(pi * abs(s - t)) - 0.19
  },
  product = function(s, t) {
    s * t - 0.25
  },
  gaussian = function(s, t) {
    exp(-((s - 0.5)^2 + (t - 0.5)^2) / 0.2)
  },
  separable = function(s, t) {
    sin(2 * pi * s) * cos(2 * pi * t)
  },
  historical = function(s, t) {
    # Only effect where s < t (historical/concurrent model)
    ifelse(s < t, cos(pi * (t - s)), 0)
  }
)

#' Smooth effect presets for s() terms
#'
#' Named list of preset f(x,t) surfaces for smooth varying coefficient terms.
#' Each function takes arguments (x, t) and returns a matrix.
#' @keywords internal
smooth_presets <- list(
  beta = function(x, t) {
    # Rescale x to [0,1] for beta density
    x_scaled <- (x - min(x)) / (max(x) - min(x) + 1e-10)
    outer(x_scaled, t, function(xi, ti) dbeta(ti, 2 + 3 * xi, 2))
  },
  dnorm = function(x, t) {
    outer(x, t, function(xi, ti) dnorm(ti, mean = 0.5, sd = 0.2) * xi)
  },
  sine = function(x, t) {
    outer(cos(x), t - 0.5, `*`)
  },
  cosine = function(x, t) {
    outer(sin(x), cos(2 * pi * t), `*`)
  },
  polynomial = function(x, t) {
    outer(x, t, function(xi, ti) xi * ti * (1 - ti))
  },
  step = function(x, t) {
    outer(x, t, function(xi, ti) ifelse(ti > 0.5, xi, -xi))
  }
)

#' Constant-over-t effect presets for c() terms
#'
#' Named list of preset f(...) functions for constant effects.
#' Single-covariate presets take argument x; multi-covariate presets
#' (e.g., `gaussian_2d` for `c(te(x1, x2))`) take multiple arguments.
#' @keywords internal
const_presets <- list(
  constant = function(x) {
    rep(1.5, length(x))
  },
  gaussian_2d = function(x1, x2) {
    # For te(x1, x2) wrapped in c()
    -x1 * x2^2
  },
  linear = function(x) {
    2 * x
  }
)

#' Linear varying coefficient presets
#'
#' Named list of preset beta(t) functions for linear varying coefficient terms.
#' Each function takes argument t and returns a vector.
#' @keywords internal
linear_presets <- list(
  dnorm = function(t) {
    # Scaled normal density centered at 0.2
    as.vector(scale(-dnorm(4 * (t - 0.2))))
  },
  sine = function(t) {
    sin(2 * pi * t)
  },
  cosine = function(t) {
    cos(2 * pi * t)
  },
  constant = function(t) {
    rep(1, length(t))
  },
  linear = function(t) {
    t - 0.5
  }
)

#' Concurrent functional effect presets
#'
#' Named list of preset beta_c(t) functions for concurrent functional covariates.
#' These are varying coefficients where Xc(t) * beta_c(t) contributes to the response.
#' Each function takes argument t and returns a vector.
#' @keywords internal
concurrent_presets <- list(
  dnorm = function(t) {
    as.vector(scale(-dnorm(4 * (t - 0.2))))
  },
  sine = function(t) {
    sin(2 * pi * t)
  },
  cosine = function(t) {
    cos(2 * pi * t)
  },
  constant = function(t) {
    rep(1, length(t))
  },
  linear = function(t) {
    t - 0.5
  }
)

#' Intercept presets for mu(t)
#'
#' Named list of preset intercept functions mu(t).
#' Each function takes argument t and returns a vector.
#' @keywords internal
intercept_presets <- list(
  constant = function(t) {
    rep(1, length(t))
  },
  beta = function(t) {
    1 + dbeta(t, 2, 7)
  },
  sine = function(t) {
    sin(2 * pi * t)
  },
  zero = function(t) {
    rep(0, length(t))
  }
)


#' Resolve an effect specification to a coefficient function/matrix
#'
#' Accept a preset string, custom function, or list with parameters,
#' and return the evaluated coefficient function/matrix.
#'
#' @param spec Effect specification: a string (preset name), function, numeric,
#'   or list with `fun` and additional parameters. The special preset `"random"`
#'   uses the random_generator to create reproducible random effects.
#' @param term_type Character string: one of "ff", "smooth", "const", "intercept",
#'   "linear", or "concurrent".
#' @param xind Numeric vector of x evaluation points (for ff terms).
#' @param yind Numeric vector of y (response) evaluation points.
#' @param random_generator Optional random truth generator from
#'   [make_random_truth_generator()], required when spec = "random".
#' @returns Evaluated coefficient: matrix for ff/smooth terms, vector for others.
#' @keywords internal
resolve_effect <- function(
  spec,
  term_type,
  xind = NULL,
  yind = NULL,
  random_generator = NULL
) {
  # Get the appropriate preset library
  presets <- switch(
    term_type,
    ff = ff_presets,
    smooth = smooth_presets,
    const = const_presets,
    intercept = intercept_presets,
    linear = linear_presets,
    concurrent = concurrent_presets,
    NULL
  )

  # Track whether spec is from random generator (returns full matrix, not scalar)
  is_random <- FALSE

  # Handle different spec types
  if (is.character(spec)) {
    if (spec == "random") {
      # Use random generator for reproducible random effects
      if (is.null(random_generator)) {
        stop("'random' preset requires random_generator")
      }
      fn <- random_generator[[term_type]]
      if (is.null(fn)) {
        stop("Random generator does not support term type '", term_type, "'")
      }
      is_random <- TRUE
    } else {
      # Standard preset name
      if (is.null(presets) || !spec %in% names(presets)) {
        stop("Unknown preset '", spec, "' for term type '", term_type, "'")
      }
      fn <- presets[[spec]]
    }
  } else if (is.function(spec)) {
    fn <- spec
  } else if (is.numeric(spec)) {
    # Constant value - return as-is for const type, replicate for others
    if (term_type == "ff") {
      return(matrix(spec, nrow = length(xind), ncol = length(yind)))
    } else if (term_type %in% c("smooth", "linear", "concurrent")) {
      return(matrix(spec, nrow = 1, ncol = length(yind)))
    } else if (term_type == "const") {
      # For const terms, return numeric as-is (will be used by compute_const_effect)
      return(spec)
    } else {
      return(rep(spec, length(yind)))
    }
  } else if (is.list(spec) && "fun" %in% names(spec)) {
    fn <- spec$fun
  } else if (is.list(spec) && "type" %in% names(spec)) {
    # Handle list syntax for concurrent: list(type = "concurrent", effect = ...)
    # Return the spec itself for later processing
    return(spec)
  } else {
    stop("Invalid effect specification for term type '", term_type, "'")
  }

  # Evaluate the function
  if (term_type == "ff") {
    if (is_random) {
      # Random generator functions take vectors and return full matrix
      fn(xind, yind)
    } else {
      # Preset functions are designed for outer() (take scalars)
      outer(xind, yind, fn)
    }
  } else if (term_type == "smooth") {
    if (is_random) {
      # Random generator returns full matrix function
      fn
    } else {
      fn
    }
  } else if (term_type == "intercept") {
    fn(yind)
  } else if (term_type == "const") {
    fn
  } else if (term_type == "concurrent") {
    # Concurrent: return function for later evaluation at yind
    fn
  } else {
    fn
  }
}


# Centering helpers for pffr identifiability ----------------------------------

#' Center functional covariate column-wise
#'
#' Centers X so that colMeans(X) = 0, i.e., each evaluation point s has
#' mean 0 across observations. This is the recommended centering for ff terms
#' so the intercept can be interpreted as the global mean function.
#'
#' @param X Numeric matrix (n x n_xind) of functional covariate.
#' @returns Centered matrix with same dimensions.
#' @keywords internal
center_functional_covariate <- function(X) {
  sweep(X, 2, colMeans(X), "-")
}

#' Compute Simpson integration weights
#'
#' Computes Simpson's rule weights for numerical integration on an equidistant
#' grid. This matches pffr's default integration method in ff() terms.
#' Note: pffr uses (b-a)/n/3 (divides by grid length n), not (b-a)/(n-1)/3.
#'
#' @param xind Numeric vector of evaluation points.
#' @returns Numeric vector of Simpson weights.
#' @keywords internal
simpson_weights <- function(xind) {
  n <- length(xind)
  if (n < 2) return(1)
  # pffr's formula: ((xind[n] - xind[1]) / n) / 3 * [1, 4, 2, ..., 1]
  ((xind[n] - xind[1]) / n / 3) * c(1, rep(c(4, 2), length.out = n - 2), 1)
}

#' Center beta(s,t) surface for ff identifiability
#'
#' Centers beta so that weighted sum over s equals 0 for each t:
#' sum_s w(s) * beta(s, t) = 0 for all t.
#' This matches pffr's identifiability constraint for ff terms.
#' Uses Simpson weights to be consistent with pffr's default integration.
#'
#' @param beta_st Numeric matrix (n_xind x n_yind) coefficient surface.
#' @param xind Numeric vector of x evaluation points (for integration weights).
#' @returns Centered beta matrix with same dimensions.
#' @keywords internal
center_ff_beta <- function(beta_st, xind) {
  # Use Simpson integration weights (consistent with ff() defaults)
  w <- simpson_weights(xind)
  w_sum <- sum(w)
  # Weighted column means
  col_means <- colSums(beta_st * w) / w_sum
  sweep(beta_st, 2, col_means, "-")
}

#' Center effect matrix column-wise for smooth term identifiability
#'
#' Centers effect matrix so that colMeans(effect) = 0, i.e., for each t,
#' the effect averages to 0 across observations: sum_i f(z_i, t) = 0 for all t.
#' This matches pffr's identifiability constraint for smooth varying coefficient terms.
#'
#' @param effect Numeric matrix (n x n_yind) of effect values.
#' @returns Centered matrix with same dimensions.
#' @keywords internal
center_smooth_effect <- function(effect) {
  sweep(effect, 2, colMeans(effect), "-")
}

#' Center concurrent functional covariate column-wise
#'
#' Centers Xc so that colMeans(Xc) = 0 for each t. This ensures the intercept
#' can be interpreted as the mean response function.
#'
#' @param Xc Numeric matrix (n x n_yind) of concurrent functional covariate.
#' @returns Centered matrix with same dimensions.
#' @keywords internal
center_concurrent_covariate <- function(Xc) {
  sweep(Xc, 2, colMeans(Xc), "-")
}

#' Compute concurrent functional effect: Xc(t) * beta_c(t)
#'
#' Computes the pointwise product of a concurrent functional covariate
#' with a varying coefficient function.
#'
#' @param Xc Numeric matrix (n x n_yind) of concurrent covariate values.
#' @param beta_t Numeric vector (length n_yind) of concurrent coefficient.
#' @returns Numeric matrix (n x n_yind) of effect contributions.
#' @keywords internal
compute_concurrent_effect <- function(Xc, beta_t) {
  sweep(Xc, 2, beta_t, "*")
}

# Effect computation helpers --------------------------------------------------

#' Compute function-on-function effect
#'
#' Compute the integral effect: integral of X(s) * beta(s,t) ds
#' using Simpson integration weights consistent with `ff()` defaults.
#'
#' @param X Numeric matrix (n x length(xind)) of functional covariate values.
#' @param beta_st Numeric matrix (length(xind) x length(yind)) coefficient surface.
#' @param xind Numeric vector of x evaluation points.
#' @returns Numeric matrix (n x length(yind)) of effect contributions.
#' @keywords internal
compute_ff_effect <- function(X, beta_st, xind) {
  n <- nrow(X)

  # Simpson integration weights (consistent with ff() defaults)
  w <- simpson_weights(xind)
  L <- matrix(w, nrow = n, ncol = length(xind), byrow = TRUE)

  # Compute integral: (L * X) %*% beta_st
  (L * X) %*% beta_st
}

#' Compute linear varying coefficient effect
#'
#' Compute effect: x * beta(t) for scalar covariate with varying coefficient.
#'
#' @param x Numeric vector (length n) of scalar covariate values.
#' @param beta_t Numeric vector (length(yind)) of varying coefficient.
#' @param yind Numeric vector of y evaluation points (unused, for API consistency).
#' @returns Numeric matrix (n x length(yind)) of effect contributions.
#' @keywords internal
compute_linear_effect <- function(x, beta_t, yind = NULL) {
  n <- length(x)
  # Outer product: x %o% beta_t
  outer(x, beta_t)
}

#' Compute smooth varying coefficient effect
#'
#' Compute effect: f(x, t) for smooth term varying over response index.
#'
#' @param x Numeric vector (length n) of covariate values.
#' @param f_xt Function or matrix. If function, takes (x, t) and returns matrix.
#'   If matrix, should be (n x length(yind)).
#' @param yind Numeric vector of y evaluation points.
#' @returns Numeric matrix (n x length(yind)) of effect contributions.
#' @keywords internal
compute_smooth_effect <- function(x, f_xt, yind) {
  n <- length(x)
  if (is.function(f_xt)) {
    f_xt(x, yind)
  } else if (is.matrix(f_xt)) {
    # If matrix has 1 row (e.g., from numeric spec), replicate to n rows
    if (nrow(f_xt) == 1 && n > 1) {
      matrix(rep(f_xt, each = n), nrow = n, ncol = ncol(f_xt))
    } else {
      f_xt
    }
  } else {
    stop("f_xt must be a function or matrix")
  }
}

#' Compute constant-over-t effect
#'
#' Compute effect: f(x1, x2, ...) that is constant over the response index t.
#' Supports multi-argument functions for terms like `c(te(x1, x2))`.
#'
#' @param ... Numeric vectors of covariate values (all same length n).
#' @param f_x Function or numeric. If function, called with all covariates as
#'   separate arguments. If numeric, used directly (length 1 or n).
#' @param n_yind Integer, number of y evaluation points.
#' @returns Numeric matrix (n x n_yind) of effect contributions (constant across columns).
#' @keywords internal
compute_const_effect <- function(..., f_x, n_yind) {
  covariates <- list(...)

  if (is.function(f_x)) {
    vals <- do.call(f_x, covariates)
  } else if (is.numeric(f_x)) {
    n <- if (length(covariates) > 0) length(covariates[[1]]) else 1
    vals <- if (length(f_x) == 1) rep(f_x, n) else f_x
  } else {
    stop("f_x must be a function or numeric")
  }

  # Replicate across yind (constant over t)
  matrix(vals, nrow = length(vals), ncol = n_yind)
}

#' Compute functional intercept
#'
#' Compute the intercept term mu(t) replicated for all observations.
#'
#' @param mu_t Numeric vector (length(yind)) of intercept function values.
#' @param n Integer, number of observations.
#' @returns Numeric matrix (n x length(mu_t)) of intercept contributions.
#' @keywords internal
compute_intercept <- function(mu_t, n) {
  matrix(mu_t, nrow = n, ncol = length(mu_t), byrow = TRUE)
}


#--------------------------------------=
# Phase 1: Covariate Generation Functions (Section 1.3)
#--------------------------------------=

#' Generate functional covariate
#'
#' Generate smooth random functions via B-spline basis for use in simulations.
#'
#' @param n Integer, number of observations (curves to generate).
#' @param xind Numeric vector of evaluation points for the functional covariate.
#' @param type Character string specifying generation method. Currently only
#'   "bspline" is supported.
#' @param bs_dim Integer, dimension of B-spline basis (default 7).
#' @param seed Optional integer seed for reproducibility within this call.
#' @returns Numeric matrix (n x length(xind)) of functional covariate values.
#' @keywords internal
generate_functional_covariate <- function(
  n,
  xind,
  type = "bspline",
  bs_dim = NULL,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  if (type != "bspline") {
    stop("Currently only type = 'bspline' is supported")
  }

  if (is.null(bs_dim)) {
    # Use a larger basis by default so simulated functional covariates are not
    # overly low-rank relative to typical ff() basis sizes.
    bs_dim <- min(5L + as.integer(floor(length(xind) / 2)), 25L)
  }

  # B-spline basis construction (adapted from existing rf() in pffrSim)
  nk <- bs_dim - 2
  xu <- max(xind)
  xl <- min(xind)
  xr <- xu - xl
  xl <- xl - xr * 0.001
  xu <- xu + xr * 0.001
  dx <- (xu - xl) / (nk - 1)
  kn <- seq(xl - dx * 3, xu + dx * 3, length = nk + 4 + 2)

  # Design matrix
  X_basis <- splines::spline.des(kn, xind, 4, xind * 0)$design

  # Generate n curves
  coefs <- matrix(rnorm(n * bs_dim), nrow = n, ncol = bs_dim)
  X <- coefs %*% t(X_basis)

  X
}

#' Generate scalar covariate
#'
#' Generate scalar covariates of various types for use in simulations.
#'
#' @param n Integer, number of observations.
#' @param type Character string specifying distribution: "normal", "uniform",
#'   or "factor".
#' @param levels For type = "factor", integer number of factor levels (default 3).
#' @param seed Optional integer seed for reproducibility within this call.
#' @returns For "normal" or "uniform": numeric vector of length n.
#'   For "factor": factor of length n.
#' @keywords internal
generate_scalar_covariate <- function(
  n,
  type = "normal",
  levels = 3,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  switch(
    type,
    normal = rnorm(n),
    uniform = runif(n),
    factor = factor(sample(seq_len(levels), n, replace = TRUE)),
    stop("Unknown covariate type: ", type)
  )
}


#--------------------------------------=
# Phase 1.4-1.5: New pffrSim with formula interface and backward compatibility
#--------------------------------------=

#' Scenario-to-formula mapping for backward compatibility
#'
#' Maps legacy scenario strings to formula and effects specifications.
#' @param scenario Character vector of scenario names.
#' @param nxgrid Number of x grid points.
#' @param nygrid Number of y grid points.
#' @param limits Optional limits function for ff terms.
#' @returns List with `formula`, `effects`, and `intercept` components.
#' @keywords internal
scenario_to_formula <- function(scenario, nxgrid, nygrid, limits = NULL) {
  # Define individual scenario components
  scenario_defs <- list(
    int = list(terms = NULL, effects = list()),
    ff = list(terms = "ff(X1, xind = s)", effects = list(X1 = "cosine")),
    lin = list(terms = "xlin", effects = list(xlin = "dnorm")),
    smoo = list(terms = "s(xsmoo)", effects = list(xsmoo = "sine")),
    te = list(
      terms = "c(te(xte1, xte2))",
      effects = list(`te(xte1,xte2)` = "gaussian_2d")
    ),
    const = list(terms = "c(xconst)", effects = list(xconst = 2)),
    factor = list(terms = "xfactor", effects = list(xfactor = "factor_effect"))
  )

  # Handle "all" scenario

  if (length(scenario) == 1 && scenario == "all") {
    scenario <- c("ff", "lin", "smoo", "te", "const", "factor")
  }

  # Combine selected scenarios
  all_terms <- character(0)
  all_effects <- list()

  for (sc in scenario) {
    if (sc == "int") next
    if (!sc %in% names(scenario_defs)) {
      stop("Unknown scenario: ", sc)
    }
    def <- scenario_defs[[sc]]
    if (!is.null(def$terms)) {
      all_terms <- c(all_terms, def$terms)
    }
    all_effects <- c(all_effects, def$effects)
  }

  # Build formula
  if (length(all_terms) == 0) {
    formula <- Y ~ 1
  } else {
    formula_str <- paste("Y ~", paste(all_terms, collapse = " + "))
    formula <- as.formula(formula_str)
  }

  list(
    formula = formula,
    effects = all_effects,
    intercept = "beta"
  )
}


#' Generate data from formula specification
#'
#' Internal function that generates simulated data based on a pffr formula.
#'
#' @param formula A formula specifying the model structure.
#' @param n Number of observations.
#' @param yind Numeric vector of y (response) evaluation points.
#' @param xind Numeric vector of x (functional covariate) evaluation points.
#' @param data Optional data frame with pre-generated covariates.
#' @param effects Named list mapping term names to effect specifications.
#'   For concurrent terms, use list syntax: `list(Xc = list(type = "concurrent", effect = "sine"))`.
#' @param intercept Intercept specification (preset name, function, or numeric).
#' @param SNR Signal-to-noise ratio.
#' @param family A family object for response distribution.
#' @param propmissing Proportion of missing values in response.
#' @param limits Optional limits function for ff terms.
#' @param wiggliness Controls smoothness for "random" preset (default: 1).
#' @param k_truth Named list of basis dimensions for random truth generation.
#' @returns A data frame with simulated data and truth attributes.
#' @keywords internal
pffrSim_formula <- function(
  formula,
  n,
  yind,
  xind,
  data = NULL,
  effects = list(),
  intercept = "beta",
  SNR = 10,
  family = gaussian(),
  propmissing = 0,
  limits = NULL,
  wiggliness = 1,
  k_truth = list()
) {
  nygrid <- length(yind)
  nxgrid <- length(xind)

  # Parse the formula
  parsed <- parse_pffr_formula(formula)

  # Check if formula requests intercept (respects Y ~ 0 + ... or Y ~ -1 + ...)
  has_intercept <- attr(terms(formula), "intercept") == 1L

  # Check if any effect uses "random"
  uses_random <- function(spec) {
    if (is.character(spec) && spec == "random") return(TRUE)
    if (is.list(spec) && isTRUE(spec$effect == "random")) return(TRUE)
    FALSE
  }
  needs_random <- any(vapply(effects, uses_random, logical(1))) ||
    (is.character(intercept) && intercept == "random")

  # Initialize data
  if (is.null(data)) {
    data <- list()
  } else {
    data <- as.list(data)
  }

  normalize_key <- function(x) gsub("[[:space:]]+", "", x)
  effects_norm <- effects
  if (length(effects_norm)) {
    names(effects_norm) <- normalize_key(names(effects_norm))
  }

  # ---------------------------------------------------------------------------

  # Phase 1: Pre-generate covariates to determine z-range for smooth terms
  # This ensures B-spline bases cover the actual data range
  # ---------------------------------------------------------------------------
  z_values_all <- numeric(0)

  for (term_label in names(parsed$terms)) {
    term_info <- parsed$terms[[term_label]]
    term_type <- term_info$type
    varname <- term_info$varname
    wrapped <- isTRUE(term_info$wrapped)

    term_key <- normalize_key(term_label)
    call_key <- if (!is.null(term_info$call)) {
      normalize_key(safeDeparse(term_info$call))
    } else {
      NULL
    }

    effect_spec <- effects[[term_label]] %||% effects_norm[[term_key]]
    if (!is.null(term_info$call)) {
      effect_spec <- effect_spec %||%
        effects[[safeDeparse(term_info$call)]] %||%
        effects_norm[[call_key]]
    }
    effect_spec <- effect_spec %||%
      effects[[varname]] %||%
      get_default_effect(term_type, varname, wrapped)

    is_concurrent <- is.list(effect_spec) &&
      isTRUE(effect_spec$type == "concurrent")

    if (is_concurrent) {
      # Concurrent: functional covariate on yind grid
      if (is.null(data[[varname]])) {
        data[[varname]] <- I(generate_functional_covariate(n, yind))
      }
    } else if (term_type == "ff" || term_type == "sff") {
      # FF: functional covariate on xind grid
      if (is.null(data[[varname]])) {
        data[[varname]] <- I(generate_functional_covariate(n, xind))
      }
    } else if (term_type == "s" && !wrapped) {
      # Smooth: scalar covariate - collect z values for B-spline range
      if (is.null(data[[varname]])) {
        data[[varname]] <- I(generate_scalar_covariate(n, "normal"))
      }
      x_centered <- data[[varname]] - mean(data[[varname]])
      data[[varname]] <- I(x_centered)
      z_values_all <- c(z_values_all, as.vector(x_centered))
    } else if (term_type == "te" || term_type == "ti" || term_type == "t2") {
      # Tensor product terms
      if (!is.null(term_info$call)) {
        vars <- all.vars(term_info$call)
      } else {
        vars <- varname
      }
      for (v in vars) {
        if (is.null(data[[v]])) {
          data[[v]] <- I(generate_scalar_covariate(n, "normal"))
        }
      }
      if (!wrapped) {
        # Varying over t - these are smooth-like, collect z values
        for (v in vars) {
          x_centered <- data[[v]] - mean(data[[v]])
          data[[v]] <- I(x_centered)
          z_values_all <- c(z_values_all, as.vector(x_centered))
        }
      }
    } else if (wrapped || term_type == "c") {
      # Constant over t
      if (is.null(data[[varname]])) {
        data[[varname]] <- I(generate_scalar_covariate(n, "normal"))
      }
    } else if (term_type == "linear") {
      # Linear varying coefficient
      if (is.null(data[[varname]])) {
        if (grepl("factor", varname, ignore.case = TRUE)) {
          data[[varname]] <- generate_scalar_covariate(n, "factor", levels = 3)
        } else {
          data[[varname]] <- I(generate_scalar_covariate(n, "normal"))
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Phase 2: Create random generator with actual z-range (if needed)
  # ---------------------------------------------------------------------------
  random_generator <- if (needs_random) {
    # Compute z_ref from actual data range with 10% padding
    if (length(z_values_all) > 0) {
      z_range <- range(z_values_all)
      z_pad <- 0.1 * diff(z_range)
      z_ref <- seq(z_range[1] - z_pad, z_range[2] + z_pad, length.out = 25)
    } else {
      z_ref <- NULL # Use default in make_random_truth_generator
    }
    make_random_truth_generator(
      yind,
      xind,
      z_ref = z_ref,
      wiggliness = wiggliness,
      k_truth = k_truth
    )
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # Phase 3: Compute effects using pre-generated covariates
  # ---------------------------------------------------------------------------
  etaTerms <- list()
  beta_coefs <- list()

  # Compute intercept only if formula includes it
  if (has_intercept) {
    mu_t <- resolve_effect(
      intercept,
      "intercept",
      yind = yind,
      random_generator = random_generator
    )
    etaTerms$intercept <- compute_intercept(mu_t, n)
    beta_coefs$intercept <- mu_t
  }

  # Process each term
  for (term_label in names(parsed$terms)) {
    term_info <- parsed$terms[[term_label]]
    term_type <- term_info$type
    varname <- term_info$varname
    wrapped <- isTRUE(term_info$wrapped)

    # Get effect specification (default to preset based on term type)
    term_key <- normalize_key(term_label)
    call_key <- if (!is.null(term_info$call)) {
      normalize_key(safeDeparse(term_info$call))
    } else {
      NULL
    }

    effect_spec <- effects[[term_label]] %||% effects_norm[[term_key]]
    if (!is.null(term_info$call)) {
      effect_spec <- effect_spec %||%
        effects[[safeDeparse(term_info$call)]] %||%
        effects_norm[[call_key]]
    }
    effect_spec <- effect_spec %||%
      effects[[varname]] %||%
      get_default_effect(term_type, varname, wrapped)

    # Check if this is a concurrent term specified via list syntax
    is_concurrent <- is.list(effect_spec) &&
      isTRUE(effect_spec$type == "concurrent")

    if (is_concurrent) {
      # Concurrent functional covariate: Xc(t) * beta_c(t)
      # Covariate already generated in Phase 1
      # Center concurrent covariate column-wise for identifiability
      Xc_centered <- center_concurrent_covariate(as.matrix(data[[varname]]))
      data[[varname]] <- I(Xc_centered)

      # Get the effect specification
      conc_effect_spec <- effect_spec$effect %||% "sine"
      beta_t <- resolve_effect(
        conc_effect_spec,
        "concurrent",
        yind = yind,
        random_generator = random_generator
      )
      if (is.function(beta_t)) {
        beta_vec <- beta_t(yind)
      } else if (is.matrix(beta_t)) {
        beta_vec <- beta_t[1, ]
      } else {
        beta_vec <- beta_t
      }
      etaTerms[[term_label]] <- compute_concurrent_effect(Xc_centered, beta_vec)
      beta_coefs[[term_label]] <- beta_vec
    } else if (term_type == "ff" || term_type == "sff") {
      # Function-on-function term - covariate already generated in Phase 1
      # Center functional covariate column-wise for identifiability
      X_centered <- center_functional_covariate(as.matrix(data[[varname]]))
      data[[varname]] <- I(X_centered)

      beta_st <- resolve_effect(
        effect_spec,
        "ff",
        xind,
        yind,
        random_generator = random_generator
      )
      if (!is.null(limits)) {
        limit_mask <- outer(xind, yind, limits)
        beta_st <- beta_st * limit_mask
      }
      # Center beta surface for pffr identifiability constraint
      beta_st <- center_ff_beta(beta_st, xind)

      etaTerms[[term_label]] <- compute_ff_effect(X_centered, beta_st, xind)
      beta_coefs[[term_label]] <- beta_st
    } else if (term_type == "s" && !wrapped) {
      # Smooth varying coefficient term - covariate already centered in Phase 1
      x_centered <- data[[varname]]

      f_xt <- resolve_effect(
        effect_spec,
        "smooth",
        yind = yind,
        random_generator = random_generator
      )
      effect_raw <- compute_smooth_effect(x_centered, f_xt, yind)
      # Center effect column-wise for pffr identifiability constraint
      etaTerms[[term_label]] <- center_smooth_effect(effect_raw)
      beta_coefs[[term_label]] <- f_xt
    } else if (term_type == "te" || term_type == "ti" || term_type == "t2") {
      # Tensor product term - extract variable names from call
      if (!is.null(term_info$call)) {
        vars <- all.vars(term_info$call)
      } else {
        vars <- varname
      }

      if (wrapped) {
        # Constant over t: c(te(x1, x2))
        for (v in vars) {
          if (is.null(data[[v]])) {
            data[[v]] <- I(generate_scalar_covariate(n, "normal"))
          }
        }
        f_x <- resolve_effect(
          effect_spec,
          "const",
          yind = yind,
          random_generator = random_generator
        )
        # Get all covariate vectors
        cov_list <- lapply(vars, function(v) data[[v]])
        etaTerms[[term_label]] <- do.call(
          compute_const_effect,
          c(cov_list, list(f_x = f_x, n_yind = nygrid))
        )
        beta_coefs[[term_label]] <- f_x
      } else {
        # Varying over t - use all variables in the effect
        for (v in vars) {
          if (is.null(data[[v]])) {
            data[[v]] <- I(generate_scalar_covariate(n, "normal"))
          }
        }
        f_xt <- resolve_effect(
          effect_spec,
          "smooth",
          yind = yind,
          random_generator = random_generator
        )

        if (is.function(f_xt)) {
          # Check if function accepts multiple x variables
          n_args <- length(formals(f_xt))
          if (n_args >= length(vars) + 1) {
            # Multi-variable function: f(x1, x2, ..., t)
            cov_vals <- lapply(vars, \(v) data[[v]])
            etaTerms[[term_label]] <- do.call(f_xt, c(cov_vals, list(yind)))
          } else {
            # Standard 2-arg function: combine covariates via product
            combined_x <- Reduce(`*`, lapply(vars, \(v) data[[v]]))
            etaTerms[[term_label]] <- f_xt(combined_x, yind)
          }
        } else {
          # Matrix effect - use first variable (for backward compat)
          etaTerms[[term_label]] <- compute_smooth_effect(
            data[[vars[1]]],
            f_xt,
            yind
          )
        }
        beta_coefs[[term_label]] <- f_xt
      }
    } else if (wrapped || term_type == "c") {
      # Constant over t term: c(x) or c(s(x))
      if (is.null(data[[varname]])) {
        data[[varname]] <- I(generate_scalar_covariate(n, "normal"))
      }
      f_x <- resolve_effect(
        effect_spec,
        "const",
        yind = yind,
        random_generator = random_generator
      )
      if (is.function(f_x)) {
        etaTerms[[term_label]] <- compute_const_effect(
          data[[varname]],
          f_x = f_x,
          n_yind = nygrid
        )
      } else {
        etaTerms[[term_label]] <- compute_const_effect(
          data[[varname]],
          f_x = f_x,
          n_yind = nygrid
        )
      }
      beta_coefs[[term_label]] <- f_x
    } else if (term_type == "linear") {
      # Linear varying coefficient term
      if (is.null(data[[varname]])) {
        # Check if it might be a factor
        if (grepl("factor", varname, ignore.case = TRUE)) {
          data[[varname]] <- generate_scalar_covariate(n, "factor", levels = 3)
        } else {
          data[[varname]] <- I(generate_scalar_covariate(n, "normal"))
        }
      }

      if (is.factor(data[[varname]])) {
        # Factor term - varying effect by level (no centering for factors)
        fac <- data[[varname]]
        fac_effect <- resolve_effect(
          effect_spec,
          "linear",
          yind = yind,
          random_generator = random_generator
        )
        if (is.function(fac_effect)) {
          # Check if function expects 2 args (custom) or 1 arg (preset)
          n_args <- length(formals(fac_effect))
          if (n_args >= 2) {
            # Custom 2-arg function: f(factor_numeric, t)
            etaTerms[[term_label]] <- fac_effect(as.numeric(fac), yind)
          } else {
            # Linear preset (1-arg): use beta(t) * factor_level
            beta_vec <- fac_effect(yind)
            etaTerms[[term_label]] <- outer(as.numeric(fac), beta_vec)
          }
        } else {
          # Default factor effect
          etaTerms[[term_label]] <- 2 *
            as.numeric(fac) +
            sin(2 * outer(as.numeric(fac), yind))
        }
        beta_coefs[[term_label]] <- "factor"
      } else {
        # Numeric linear term with varying coefficient
        # Center scalar covariate for identifiability (recommended by pffr docs)
        x_centered <- data[[varname]] - mean(data[[varname]])
        data[[varname]] <- I(x_centered)

        beta_t <- resolve_effect(
          effect_spec,
          "linear",
          yind = yind,
          random_generator = random_generator
        )
        if (is.function(beta_t)) {
          beta_vec <- beta_t(yind)
        } else if (is.matrix(beta_t)) {
          beta_vec <- beta_t[1, ]
        } else {
          beta_vec <- beta_t
        }
        etaTerms[[term_label]] <- compute_linear_effect(x_centered, beta_vec)
        beta_coefs[[term_label]] <- beta_vec
      }
    }
  }

  # Compute eta (linear predictor)
  if (length(etaTerms) == 0) {
    # No terms (e.g., Y ~ 0): eta is zero
    eta <- matrix(0, nrow = n, ncol = nygrid)
  } else {
    eta <- Reduce(`+`, etaTerms)
  }

  # Generate noise and response
  if (family$family == "gaussian" || family$family == "gaulss") {
    # Gaussian (or gaulss location-scale): use SNR definition with normal errors
    eta_var <- var(as.vector(eta))
    # Guard against zero variance (e.g., no signal): use unit variance noise
    sigma2 <- if (eta_var > 0) eta_var / SNR else 1
    epsilon <- matrix(rnorm(n * nygrid, sd = sqrt(sigma2)), nrow = n)
    Y <- eta + epsilon
  } else if (family$family == "scaled t") {
    # Scaled t-distribution: use SNR definition with t-distributed errors
    # Default df = 4 for moderately heavy tails
    # scat() stores df in theta[1], but we use a sensible default
    df <- tryCatch(
      {
        theta <- family$getTheta(TRUE)
        if (length(theta) > 0 && is.finite(theta[1]) && theta[1] > 0) {
          theta[1]
        } else {
          4
        }
      },
      error = \(e) 4
    )
    eta_var <- var(as.vector(eta))
    sigma2 <- if (eta_var > 0) eta_var / SNR else 1
    # Scale t-distributed errors to have variance = sigma2
    # Var(t_df) = df/(df-2) for df > 2, so scale by sqrt(sigma2 * (df-2)/df)
    if (df > 2) {
      scale_factor <- sqrt(sigma2 * (df - 2) / df)
    } else {
      scale_factor <- sqrt(sigma2) # Fallback for df <= 2
    }
    epsilon <- matrix(rt(n * nygrid, df = df) * scale_factor, nrow = n)
    Y <- eta + epsilon
  } else {
    # Non-Gaussian: sample from the actual family distribution
    mu <- family$linkinv(eta)
    Y <- simulate_from_family(mu, family, n, nygrid)
    # For non-Gaussian, epsilon represents deviation from mean (less meaningful)
    epsilon <- Y - mu
  }

  # Use response name from formula (default to "Y" if not specified)
  resp_name <- parsed$response %||% "Y"
  data[[resp_name]] <- I(Y)

  # Create result structure
  result_data <- as.data.frame(data, row.names = seq_len(n))

  if (propmissing == 0) {
    result <- structure(
      result_data,
      xindex = xind,
      yindex = yind,
      truth = list(
        eta = eta,
        etaTerms = etaTerms,
        beta = beta_coefs,
        epsilon = epsilon
      )
    )
  } else {
    # Handle missing data
    n_total <- n * nygrid
    n_missing <- round(propmissing * n_total)
    missing <- rep(FALSE, n_total)
    missing[sample(n_total, n_missing)] <- TRUE

    ydata <- data.frame(
      .obs = rep(seq_len(n), each = nygrid)[!missing],
      .index = rep(yind, times = n)[!missing],
      .value = as.vector(t(result_data[[resp_name]]))[!missing]
    )

    result <- structure(
      list(data = result_data, ydata = ydata),
      xindex = xind,
      yindex = yind,
      truth = list(
        eta = eta,
        etaTerms = etaTerms,
        beta = beta_coefs,
        epsilon = epsilon
      )
    )
  }

  result
}


#' Simulate response from a GLM family
#'
#' Draw random samples from the distribution specified by a GLM family object.
#'
#' @param mu Numeric matrix (n x nygrid) of mean values on the response scale.
#' @param family A family object (e.g., binomial(), poisson()).
#' @param n Number of observations.
#' @param nygrid Number of response grid points.
#' @returns Numeric matrix (n x nygrid) of simulated response values.
#' @importFrom stats rbinom rpois rgamma rt
#' @keywords internal
simulate_from_family <- function(mu, family, n, nygrid) {
  # Supported families for simulation
  supported <- c(
    "binomial",
    "poisson",
    "Gamma",
    "quasibinomial",
    "quasipoisson"
  )

  if (!family$family %in% supported) {
    stop(
      "pffrSim does not support family '",
      family$family,
      "' for simulation.\n",
      "Supported families: ",
      paste(supported, collapse = ", "),
      ".\n",
      "Use family = gaussian() for continuous responses.",
      call. = FALSE
    )
  }

  mu_vec <- as.vector(mu)
  n_total <- n * nygrid

  Y_vec <- switch(
    family$family,
    binomial = {
      # Binomial with n=1 (Bernoulli)
      rbinom(n_total, size = 1, prob = mu_vec)
    },
    poisson = {
      rpois(n_total, lambda = mu_vec)
    },
    Gamma = {
      # Gamma: mean = mu, use shape = 1/dispersion, scale = mu * dispersion
      # Default dispersion = 1, so shape = 1, scale = mu
      rgamma(n_total, shape = 1, scale = mu_vec)
    },
    quasibinomial = {
      rbinom(n_total, size = 1, prob = mu_vec)
    },
    quasipoisson = {
      rpois(n_total, lambda = mu_vec)
    }
  )

  matrix(Y_vec, nrow = n, ncol = nygrid)
}


#' Get default effect specification for a term type
#'
#' @param term_type Character string of term type.
#' @param varname Character string of variable name.
#' @returns Default effect specification.
#' @keywords internal
get_default_effect <- function(term_type, varname, wrapped = FALSE) {
  # For te/ti/t2: use gaussian_2d when wrapped (const_presets),

  # sine when unwrapped (smooth_presets)
  if (term_type %in% c("te", "ti", "t2")) {
    return(if (wrapped) "gaussian_2d" else "sine")
  }
  switch(
    term_type,
    ff = "cosine",
    sff = "cosine",
    s = "sine",
    linear = "dnorm",
    "cosine"
  )
}
