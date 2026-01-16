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

#' Create short labels for pffr terms at model fit time
#'
#' Creates a named vector mapping mgcv smooth labels to human-readable pffr labels.
#' This is called during `pffr()` fitting and stored in `object$pffr$shortlabels`.
#'
#' @param labelmap Named list mapping original term names to mgcv smooth labels.
#' @param m.smooth List of smooth objects from the fitted model.
#' @param yindname Name of the y-index variable.
#' @param where.specials List indicating which terms are which type (par, ff, ffpc, etc.).
#' @param family The model family object (used to detect location-scale families).
#' @returns Named character vector: names are mgcv smooth labels, values are short labels.
#' @keywords internal
create_shortlabels <- function(
  labelmap,
  m.smooth,
  yindname,
  where.specials,
  family = NULL
) {
  # Build result: one short label per smooth in m.smooth

  shortlabels <- character()

  for (term_name in names(labelmap)) {
    mgcv_labels <- labelmap[[term_name]]

    # Skip if no labels (shouldn't happen)
    if (length(mgcv_labels) == 0 || all(is.na(mgcv_labels))) next

    # Check term type based on name patterns and where.specials
    term_idx <- match(term_name, names(labelmap))
    is_par <- term_idx %in% where.specials$par
    is_ffpc <- grepl("^ffpc", term_name)
    is_ff <- term_idx %in% where.specials$ff
    is_sff <- term_idx %in% where.specials$sff

    for (i in seq_along(mgcv_labels)) {
      mgcv_label <- mgcv_labels[i]
      if (is.na(mgcv_label)) next

      sm <- NULL
      if (!is.null(names(m.smooth)) && mgcv_label %in% names(m.smooth)) {
        sm <- m.smooth[[mgcv_label]]
      } else {
        sm_idx <- which(vapply(
          m.smooth,
          \(x) identical(x$label, mgcv_label),
          logical(1)
        ))
        if (length(sm_idx) > 0) {
          sm <- m.smooth[[sm_idx[1]]]
        }
      }

      # Generate short label based on term type
      if (is_par && length(mgcv_labels) > 1) {
        # Factor varying coefficient: extract factor level from smooth
        if (!is.null(sm) && !is.null(sm$by.level)) {
          short <- paste0(sm$by, sm$by.level, "(", yindname, ")")
        } else {
          short <- paste0(term_name, "(", yindname, ")")
        }
      } else if (is_par) {
        # Scalar varying coefficient
        short <- paste0(term_name, "(", yindname, ")")
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
          short <- simplify_term_label(term_name, yindname)
        }
      } else if (is_ff || is_sff) {
        # ff/sff: use simplified version of original term
        short <- simplify_term_label(term_name, yindname)
      } else {
        # Other terms (s, te, t2, c-wrapped, intercept)
        short <- simplify_term_label(term_name, yindname)
      }

      shortlabels[mgcv_label] <- short
    }
  }

  # Handle location-scale family smooths (e.g., gaulss)
  # These have labels like "s.1(...)" for scale parameter and aren't in labelmap
  # Only do this for families with multiple linear predictors
  is_location_scale <- !is.null(family) &&
    !is.null(family$nlp) &&
    family$nlp > 1

  if (is_location_scale) {
    all_smooth_labels <- vapply(m.smooth, \(x) x$label, character(1))
    unlabeled <- setdiff(all_smooth_labels, names(shortlabels))

    for (mgcv_label in unlabeled) {
      # Check for scale parameter pattern: smooth.N(...) where N > 0
      # Only match known mgcv smooth types: s, te, ti, t2
      # e.g., "s.1(t_grid.vec)", "te.1(x,t)", "s.1(t_grid.vec):xlin"
      scale_match <- regexpr("^(s|te|ti|t2)\\.([0-9]+)\\(", mgcv_label)
      if (scale_match > 0) {
        # Extract the base term by removing the .N suffix
        base_label <- sub("^(s|te|ti|t2)\\.([0-9]+)", "\\1", mgcv_label)

        # Check if this is just the intercept (yindname only) or has covariates
        has_by <- grepl(":", mgcv_label)

        if (has_by) {
          # Complex case: scale depends on covariate
          # Extract the covariate part after the colon
          by_part <- sub(".*:", "", mgcv_label)
          shortlabels[mgcv_label] <- paste0(
            "log(SD): ",
            by_part,
            "(",
            yindname,
            ")"
          )
        } else {
          # Simple case: intercept-only scale
          shortlabels[mgcv_label] <- paste0("log(SD)(", yindname, ")")
        }
      }
    }
  }

  shortlabels
}

#' Simplify a term label for display
#'
#' @param term_name Original term name string.
#' @param yindname Name of the y-index variable.
#' @returns Simplified label string.
#' @keywords internal
simplify_term_label <- function(term_name, yindname) {
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
    # Simple variable name: add (yindname)
    return(paste0(term_name, "(", yindname, ")"))
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
        return(paste0(fn_name, "(", vars[1], ",", yindname, ")"))
      } else {
        return(paste0(
          fn_name,
          "(",
          paste(vars, collapse = ","),
          ",",
          yindname,
          ")"
        ))
      }
    } else if (fn_name %in% c("ff", "sff", "ffpc")) {
      return(paste0(fn_name, "(", vars[1], ")"))
    } else {
      # Default: just use the function call with main variables
      if (length(vars) == 1) {
        return(paste0(vars[1], "(", yindname, ")"))
      } else {
        return(paste0(fn_name, "(", paste(vars, collapse = ","), ")"))
      }
    }
  }

  # Fallback: return as-is with yindname
  paste0(term_name, "(", yindname, ")")
}

#' Simulate example data for pffr
#'
#' Simulates example data for \code{\link{pffr}} from a variety of terms.
#' Supports both a formula-based interface (recommended) and a legacy
#' scenario-based interface for backward compatibility.
#'
#' **Formula Interface (Recommended):**
#' Specify a pffr-style formula and optional effect specifications:
#' \preformatted{
#' pffr_simulate(Y ~ ff(X1) + xlin + s(xsmoo), n = 100,
#'               effects = list(X1 = "cosine", xlin = "dnorm"))
#' }
#'
#' **Scenario Interface (Deprecated):**
#' Scenario "all" generates data from a complex multivariate model
#' \deqn{Y_i(t) = \mu(t) + \int X_{1i}(s)\beta_1(s,t)ds + xlin \beta_3(t) +
#' f(xte1, xte2) + f(xsmoo, t) + \beta_4 xconst + f(xfactor, t) + \epsilon_i(t)}
#'
#' Scenarios "int", "ff", "lin", "te", "smoo", "const", "factor" generate data
#' from simpler models containing only the respective terms.
#'
#' Sparse/irregular response trajectories can be generated by setting
#' \code{propmissing} > 0.
#'
#' @param formula A formula specifying the model structure (e.g.,
#'   \code{Y ~ ff(X1) + xlin}). If provided, the \code{scenario} argument is
#'   ignored.
#' @param scenario \strong{Deprecated.} Character string or vector specifying
#'   predefined scenarios. Use the \code{formula} argument instead.
#' @param n Number of observations.
#' @param nxgrid Number of evaluation points for functional covariates.
#'   Ignored if \code{xind} is provided.
#' @param nygrid Number of evaluation points for the functional response.
#'   Ignored if \code{yind} is provided.
#' @param yind Numeric vector of evaluation points for the response.
#'   Defaults to \code{seq(0, 1, length.out = nygrid)}.
#' @param xind Numeric vector of evaluation points for functional covariates.
#'   Defaults to \code{seq(0, 1, length.out = nxgrid)}.
#' @param data Optional data frame with pre-generated covariates.
#' @param effects Named list mapping term labels to effect specifications.
#'   Each entry can be a preset name (e.g., "cosine"), a function, or a numeric
#'   value. See Details.
#' @param intercept Intercept specification: preset name ("beta", "constant",
#'   "sine", "zero"), a function of \code{t}, or a numeric value.
#' @param SNR Signal-to-noise ratio: \code{var(eta) / var(epsilon)}.
#' @param family A family object for the response distribution. Defaults to
#'   \code{gaussian()}.
#' @param propmissing Proportion of missing data in the response (0 to 1).
#' @param limits A function defining integration limits for \code{ff()} terms,
#'   e.g., \code{function(s, t) s < t}.
#' @param seed Optional random seed for reproducibility.
#'
#' @returns A data frame (or list if \code{propmissing > 0}) with simulated
#'   data and attributes:
#'   \describe{
#'     \item{xindex}{Evaluation points for functional covariates}
#'     \item{yindex}{Evaluation points for the response}
#'     \item{truth}{List containing \code{eta} (linear predictor),
#'       \code{etaTerms} (individual term contributions), \code{beta}
#'       (coefficient functions), and \code{epsilon} (noise)}
#'   }
#'
#' @section Effect Presets:
#' For \code{ff()} terms: "cosine", "product", "gaussian", "separable", "historical"
#'
#' For \code{s()} terms: "beta", "dnorm", "sine", "cosine", "polynomial", "step"
#'
#' For \code{c()} terms: "constant", "gaussian_2d", "linear"
#'
#' For intercept: "constant", "beta", "sine", "zero"
#'
#' @export
#' @importFrom splines spline.des
#' @importFrom stats dnorm rnorm runif dbeta gaussian var as.formula
#'
#' @examples
#' # Formula interface
#' dat <- pffr_simulate(Y ~ ff(X1) + xlin + s(xsmoo), n = 50,
#'                      effects = list(X1 = "cosine", xlin = "dnorm", xsmoo = "sine"))
#'
#' # Legacy scenario interface (deprecated)
#' dat_legacy <- suppressWarnings(pffr_simulate(scenario = "ff", n = 50))
#'
#' # Access true coefficients
#' truth <- attr(dat, "truth")
#' str(truth$beta)
pffr_simulate <- function(
  formula = NULL,
  scenario = NULL,
  n = 100,
  nxgrid = 40,
  nygrid = 60,
  yind = NULL,
  xind = NULL,
  data = NULL,
  effects = list(),
  intercept = "beta",
  SNR = 10,
  family = gaussian(),
  propmissing = 0,
  limits = NULL,
  seed = NULL
) {
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Set up grid indices
  if (is.null(yind)) yind <- seq(0, 1, length.out = nygrid)
  if (is.null(xind)) xind <- seq(0, 1, length.out = nxgrid)

  # Determine which interface to use
  if (!is.null(formula)) {
    # New formula interface
    return(pffrSim_formula(
      formula = formula,
      n = n,
      yind = yind,
      xind = xind,
      data = data,
      effects = effects,
      intercept = intercept,
      SNR = SNR,
      family = family,
      propmissing = propmissing,
      limits = limits
    ))
  }

  # Legacy scenario interface
  if (is.null(scenario)) {
    scenario <- "all"
  }

  # Issue deprecation warning once per session (not on every call)
  if (is.null(getOption("refund.pffr_simulate.scenario.warned"))) {
    warning(
      "The 'scenario' argument in pffr_simulate() is deprecated.\n",
      "Use the formula interface instead, e.g.:\n",
      "  pffr_simulate(Y ~ ff(X1) + xlin, effects = list(X1 = 'cosine'))",
      call. = FALSE
    )
    options(refund.pffr_simulate.scenario.warned = TRUE)
  }

  # Use legacy implementation for backward compatibility
  pffrSim_legacy(
    scenario = scenario,
    n = n,
    nxgrid = length(xind),
    nygrid = length(yind),
    SNR = SNR,
    propmissing = propmissing,
    limits = limits
  )
}


#' Simulate example data for pffr (deprecated)
#'
#' @description
#' **Deprecated**
#'
#' `pffrSim()` was renamed to [pffr_simulate()] for consistency with the
#' package naming conventions.
#'
#' @inheritParams pffr_simulate
#' @export
#' @keywords internal
pffrSim <- function(
  formula = NULL,
  scenario = NULL,
  n = 100,

  nxgrid = 40,
  nygrid = 60,

  yind = NULL,
  xind = NULL,
  data = NULL,
  effects = list(),
  intercept = "beta",
  SNR = 10,
  family = gaussian(),
  propmissing = 0,
  limits = NULL,
  seed = NULL
) {
  .Deprecated("pffr_simulate")
  pffr_simulate(
    formula = formula,
    scenario = scenario,
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    yind = yind,
    xind = xind,
    data = data,
    effects = effects,
    intercept = intercept,
    SNR = SNR,
    family = family,
    propmissing = propmissing,
    limits = limits,
    seed = seed
  )
}


#' Legacy pffr_simulate implementation
#'
#' Internal function preserving the original pffr_simulate behavior for backward
#' compatibility. This is called when the deprecated \code{scenario} argument
#' is used.
#'
#' @inheritParams pffr_simulate
#' @keywords internal
pffrSim_legacy <- function(
  scenario = "all",
  n = 100,
  nxgrid = 40,
  nygrid = 60,
  SNR = 10,
  propmissing = 0,
  limits = NULL
) {
  mc <- match.call()

  ## generates random functions...
  rf <- function(x = seq(0, 1, length = 100), bs.dim = 7, center = FALSE) {
    nk <- bs.dim - 2
    xu <- max(x)
    xl <- min(x)
    xr <- xu - xl
    xl <- xl - xr * 0.001
    xu <- xu + xr * 0.001
    dx <- (xu - xl) / (nk - 1)
    kn <- seq(xl - dx * 3, xu + dx * 3, length = nk + 4 + 2)
    X <- splines::spline.des(kn, x, 4, x * 0)$design
    drop(X %*% rnorm(bs.dim))
  }

  test1 <- function(s, t) {
    s * cos(pi * abs(s - t)) - 0.19
  }

  s <- seq(0, 1, length = nxgrid)
  t <- seq(0, 1, length = nygrid)

  mu.t <- matrix(1 + dbeta(t, 2, 7), nrow = n, ncol = nygrid, byrow = TRUE)

  data <- list()
  etaTerms <- list()
  etaTerms$int <- mu.t

  # Functional covariates
  data$X1 <- I(t(replicate(n, rf(s))))

  L <- matrix(1 / nxgrid, ncol = nxgrid, nrow = n)
  LX1 <- L * data$X1
  beta1.st <- outer(s, t, test1)
  if (!is.null(limits)) {
    range <- outer(s, t, limits)
    beta1.st <- beta1.st * range
  }
  etaTerms$X1 <- LX1 %*% beta1.st

  # Scalar covariates
  data$xlin <- I(rnorm(n))
  beta.t <- matrix(
    scale(-dnorm(4 * (t - 0.2))),
    nrow = n,
    ncol = nygrid,
    byrow = TRUE
  )
  etaTerms$xlin <- data$xlin * beta.t

  data$xsmoo <- I(rnorm(n))
  etaTerms$xsmoo <- outer(drop(scale(cos(data$xsmoo))), (t - 0.5), "*")

  data$xfactor <- sample(gl(3, n / 3), replace = TRUE)
  etaTerms$xfactor <- 2 *
    as.numeric(data$xfactor) +
    sin(2 * outer(as.numeric(data$xfactor), t))
  if ("2factor" %in% scenario) {
    data$x2factor <- sample(gl(3, n / 3), replace = TRUE)
    etaTerms$x2factor <- 2 *
      as.numeric(data$x2factor) +
      cos(2 * outer(as.numeric(data$x2factor), t))
  }

  data$xte1 <- I(rnorm(n))
  data$xte2 <- I(rnorm(n))
  etaTerms$xte <- matrix(
    drop(scale(-data$xte1 * data$xte2^2)),
    ncol = nygrid,
    nrow = n
  )

  data$xconst <- I(rnorm(n))
  etaTerms$xconst <- matrix(2 * data$xconst, ncol = nygrid, nrow = n)

  if (length(scenario) == 1) {
    eta <- mu.t +
      switch(
        scenario,
        "int" = 0,
        "all" = Reduce("+", etaTerms),
        "ff" = etaTerms$X1,
        "lin" = etaTerms$xlin,
        "smoo" = etaTerms$xsmoo,
        "te" = etaTerms$xte,
        "const" = etaTerms$xconst,
        "factor" = etaTerms$xfactor
      )
  } else {
    stopifnot(all(
      scenario %in%
        c("int", "ff", "lin", "smoo", "te", "const", "factor", "2factor")
    ))
    eta <- 0 * mu.t
    if ("int" %in% scenario) eta <- eta + mu.t
    if ("ff" %in% scenario) eta <- eta + etaTerms$X1
    if ("lin" %in% scenario) eta <- eta + etaTerms$xlin
    if ("smoo" %in% scenario) eta <- eta + etaTerms$xsmoo
    if ("te" %in% scenario) eta <- eta + etaTerms$xte
    if ("const" %in% scenario) eta <- eta + etaTerms$xconst
    if ("factor" %in% scenario) eta <- eta + etaTerms$xfactor
    if ("2factor" %in% scenario) eta <- eta + etaTerms$x2factor
  }

  eps <- sd(as.vector(eta)) /
    sqrt(SNR) *
    matrix(scale(rnorm(n * nygrid)), nrow = n)
  data$Y <- I(eta + eps)

  if (propmissing == 0) {
    return(structure(
      as.data.frame(data, row.names = 1:n),
      xindex = s,
      yindex = t,
      truth = list(eta = eta, etaTerms = etaTerms),
      call = mc
    ))
  } else {
    missing <- sample(c(
      rep(TRUE, propmissing * n * nygrid),
      rep(FALSE, n * nygrid - propmissing * n * nygrid)
    ))
    data <- as.data.frame(data, row.names = 1:n)

    ydata <- data.frame(
      .obs = rep(1:n, each = nygrid)[!missing],
      .index = rep(t, times = n)[!missing],
      .value = as.vector(t(data$Y))[!missing]
    )

    return(structure(
      list(data = data, ydata = ydata),
      xindex = s,
      yindex = t,
      truth = list(eta = eta, etaTerms = etaTerms),
      call = mc
    ))
  }
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


# ==============================================================================
# Phase 1: pffrSim Redesign - Core Infrastructure (Section 1.1)
# ==============================================================================

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
#'   or list with `fun` and additional parameters.
#' @param term_type Character string: one of "ff", "smooth", "const", "intercept", "linear".
#' @param xind Numeric vector of x evaluation points (for ff terms).
#' @param yind Numeric vector of y (response) evaluation points.
#' @returns Evaluated coefficient: matrix for ff/smooth terms, vector for others.
#' @keywords internal
resolve_effect <- function(spec, term_type, xind = NULL, yind = NULL) {
  # Get the appropriate preset library

  presets <- switch(
    term_type,
    ff = ff_presets,
    smooth = smooth_presets,
    const = const_presets,
    intercept = intercept_presets,
    linear = linear_presets,
    NULL
  )

  # Handle different spec types
  if (is.character(spec)) {
    # Preset name
    if (is.null(presets) || !spec %in% names(presets)) {
      stop("Unknown preset '", spec, "' for term type '", term_type, "'")
    }
    fn <- presets[[spec]]
  } else if (is.function(spec)) {
    fn <- spec
  } else if (is.numeric(spec)) {
    # Constant value - return as-is for const type, replicate for others
    if (term_type == "ff") {
      return(matrix(spec, nrow = length(xind), ncol = length(yind)))
    } else if (term_type %in% c("smooth", "linear")) {
      return(matrix(spec, nrow = 1, ncol = length(yind)))
    } else if (term_type == "const") {
      # For const terms, return numeric as-is (will be used by compute_const_effect)
      return(spec)
    } else {
      return(rep(spec, length(yind)))
    }
  } else if (is.list(spec) && "fun" %in% names(spec)) {
    fn <- spec$fun
  } else {
    stop("Invalid effect specification for term type '", term_type, "'")
  }

  # Evaluate the function
  if (term_type == "ff") {
    outer(xind, yind, fn)
  } else if (term_type == "intercept") {
    fn(yind)
  } else if (term_type == "const") {
    fn
  } else {
    fn
  }
}


# ==============================================================================
# Phase 1: Effect Computation Functions (Section 1.2)
# ==============================================================================

#' Compute function-on-function effect
#'
#' Compute the integral effect: integral of X(s) * beta(s,t) ds
#' using rectangular integration weights consistent with `ff()` defaults.
#'
#' @param X Numeric matrix (n x length(xind)) of functional covariate values.
#' @param beta_st Numeric matrix (length(xind) x length(yind)) coefficient surface.
#' @param xind Numeric vector of x evaluation points.
#' @returns Numeric matrix (n x length(yind)) of effect contributions.
#' @keywords internal
compute_ff_effect <- function(X, beta_st, xind) {
  n <- nrow(X)
  n_xind <- length(xind)

  # Rectangular integration weights (consistent with ff() defaults)
  dx <- if (length(xind) > 1) diff(xind)[1] else 1
  L <- matrix(dx, nrow = n, ncol = n_xind)

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


# ==============================================================================
# Phase 1: Covariate Generation Functions (Section 1.3)
# ==============================================================================

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
  bs_dim = 7,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  if (type != "bspline") {
    stop("Currently only type = 'bspline' is supported")
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


# ==============================================================================
# Phase 1.4-1.5: New pffrSim with formula interface and backward compatibility
# ==============================================================================

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
#' @param intercept Intercept specification (preset name, function, or numeric).
#' @param SNR Signal-to-noise ratio.
#' @param family A family object for response distribution.
#' @param propmissing Proportion of missing values in response.
#' @param limits Optional limits function for ff terms.
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
  limits = NULL
) {
  nygrid <- length(yind)
  nxgrid <- length(xind)

  # Parse the formula
  parsed <- parse_pffr_formula(formula)

  # Check if formula requests intercept (respects Y ~ 0 + ... or Y ~ -1 + ...)
  has_intercept <- attr(terms(formula), "intercept") == 1L

  # Initialize data and eta terms
  if (is.null(data)) {
    data <- list()
  } else {
    data <- as.list(data)
  }
  etaTerms <- list()
  beta_coefs <- list()

  # Compute intercept only if formula includes it
  if (has_intercept) {
    mu_t <- resolve_effect(intercept, "intercept", yind = yind)
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
    effect_spec <- effects[[term_label]] %||%
      effects[[varname]] %||%
      get_default_effect(term_type, varname, wrapped)

    if (term_type == "ff" || term_type == "sff") {
      # Function-on-function term
      if (is.null(data[[varname]])) {
        data[[varname]] <- I(generate_functional_covariate(n, xind))
      }
      beta_st <- resolve_effect(effect_spec, "ff", xind, yind)
      if (!is.null(limits)) {
        limit_mask <- outer(xind, yind, limits)
        beta_st <- beta_st * limit_mask
      }
      etaTerms[[term_label]] <- compute_ff_effect(
        data[[varname]],
        beta_st,
        xind
      )
      beta_coefs[[term_label]] <- beta_st
    } else if (term_type == "s" && !wrapped) {
      # Smooth varying coefficient term
      if (is.null(data[[varname]])) {
        data[[varname]] <- I(generate_scalar_covariate(n, "normal"))
      }
      f_xt <- resolve_effect(effect_spec, "smooth", yind = yind)
      etaTerms[[term_label]] <- compute_smooth_effect(
        data[[varname]],
        f_xt,
        yind
      )
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
        f_x <- resolve_effect(effect_spec, "const", yind = yind)
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
        f_xt <- resolve_effect(effect_spec, "smooth", yind = yind)

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
      f_x <- resolve_effect(effect_spec, "const", yind = yind)
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
        # Factor term - varying effect by level
        fac <- data[[varname]]
        fac_effect <- resolve_effect(effect_spec, "linear", yind = yind)
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
        beta_t <- resolve_effect(effect_spec, "linear", yind = yind)
        if (is.function(beta_t)) {
          beta_vec <- beta_t(yind)
        } else if (is.matrix(beta_t)) {
          beta_vec <- beta_t[1, ]
        } else {
          beta_vec <- beta_t
        }
        etaTerms[[term_label]] <- compute_linear_effect(
          data[[varname]],
          beta_vec
        )
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
