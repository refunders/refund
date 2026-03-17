#' Modular Formula Parsing for pffr
#'
#' These functions provide a modular implementation of pffr formula parsing,
#' extracting logical units from the main pffr() function for better
#' maintainability and testing.
#'
#' @name pffr-formula
#' @keywords internal
NULL

# Constants for pffr special terms
PFFR_SPECIALS <- c("s", "te", "ti", "t2", "ff", "c", "sff", "ffpc", "pcre")

#' Parse a pffr model formula and classify terms
#'
#' Internal function for pffr() model fitting. Different from the simulation
#' formula parser used by pffrSim.
#'
#' @param formula The pffr formula.
#' @param data The data list/frame (unused, for API consistency).
#' @param ydata Optional sparse/irregular response data (unused).
#' @returns A list with:
#'   - `tf`: The terms object
#'   - `term_strings`: Character vector of term labels
#'   - `terms`: Named list of parsed term calls
#'   - `frml_env`: Formula environment
#'   - `where_specials`: Named list of term indices by type
#'   - `response_name`: Response variable name
#'   - `has_intercept`: Logical, whether formula has intercept
#' @keywords internal
parse_pffr_model_formula <- function(formula, data = NULL, ydata = NULL) {
  tf <- terms.formula(formula, specials = PFFR_SPECIALS)
  term_strings <- attr(tf, "term.labels")

  # Parse each term string into a call
  terms <- lapply(term_strings, \(trm) as.call(parse(text = trm))[[1]])
  names(terms) <- term_strings

  # Classify terms by type
  where_specials <- lapply(
    PFFR_SPECIALS,
    \(sp) attr(tf, "specials")[[sp]] - 1
  )
  names(where_specials) <- PFFR_SPECIALS

  # Parametric terms are those not in any special category
  all_special_indices <- unlist(attr(tf, "specials")) - 1
  where_specials$par <- if (length(term_strings)) {
    which(!seq_along(term_strings) %in% all_special_indices)
  } else {
    numeric(0)
  }

  list(
    tf = tf,
    term_strings = term_strings,
    terms = terms,
    frml_env = environment(formula),
    where_specials = where_specials,
    response_name = attr(tf, "variables")[2][[1]],
    has_intercept = attr(tf, "intercept") == 1
  )
}

#' Transform intercept term to functional intercept
#'
#' Creates an s() smooth term for the functional intercept.
#'
#' @param yind_vec_name Name of the y-index vector variable (character).
#' @param bs_int Spline basis specification for intercept (named list).
#' @param yind_name Name of y-index for labeling.
#' @returns List with `term_string` and `label`.
#' @keywords internal
transform_intercept_term <- function(yind_vec_name, bs_int, yind_name) {
  # Build: s(yind_vec_name, bs=..., k=..., m=...)
  intcall <- as.call(c(
    list(as.name("s"), as.name(yind_vec_name)),
    bs_int
  ))

  intstring <- safeDeparse(intcall)
  names(intstring) <- paste0("Intercept(", yind_name, ")")

  list(term_string = intstring, label = names(intstring))
}

#' Transform c() wrapped terms by removing the wrapper
#'
#' @param term_string Original term string with c() wrapper.
#' @returns Term string with c() removed.
#' @keywords internal
transform_c_term <- function(term_string) {
  sub("\\)$", "", sub("^c\\(", "", term_string))
}

#' Transform smooth terms (s, te, ti, t2) to include y-index
#'
#' Converts univariate/bivariate smooth terms into tensor products that
#' include the functional response index.
#'
#' @param term The term call object.
#' @param yind_vec_name Name of y-index vector (character).
#' @param bs_yindex Default basis spec for y-index (named list).
#' @param tensortype Type of tensor product (symbol: ti or t2).
#' @param algorithm Fitting algorithm name (character or symbol).
#' @returns Transformed term string.
#' @keywords internal
transform_smooth_term <- function(
  term,
  yind_vec_name,
  bs_yindex,
  tensortype,
  algorithm
) {
  xnew <- term
  term_type <- deparse(term[[1]])
  is_gamm4 <- as.character(algorithm) == "gamm4"

  # Convert term type for gamm4 compatibility
  xnew[[1]] <- get_smooth_type(term_type, tensortype, is_gamm4)

  # Set dimension parameter
  xnew$d <- get_smooth_dimensions(term, term_type)
  n_base_marginals <- length(xnew$d) - 1L

  # Add y-index as last variable
  xnew[[length(xnew) + 1]] <- as.name(yind_vec_name)

  # Get term-specific or default bs.yindex
  this_bs_yindex <- term$bs.yindex %||% bs_yindex
  if (!is.null(term$bs.yindex)) {
    this_bs_yindex <- eval(term$bs.yindex)
  }
  xnew <- xnew[names(xnew) != "bs.yindex"]

  # Set constraints for ti terms
  if (deparse(xnew[[1]]) == "ti") {
    xnew$mc <- c(rep(TRUE, length(xnew$d) - 1), FALSE)
  }

  # Set spline parameters
  xnew$bs <- get_smooth_bs(
    term,
    this_bs_yindex,
    n_base_marginals = n_base_marginals
  )
  m_spec <- get_smooth_m(term, this_bs_yindex, n_marginals = length(xnew$bs))
  if (!is.null(m_spec)) xnew$m <- m_spec
  xnew$k <- get_smooth_k(term, xnew$d, this_bs_yindex)

  # Preserve xt if present
  if ("xt" %in% names(term)) {
    xnew$xt <- term$xt
  }

  safeDeparse(xnew)
}

# Determine smooth term type based on original type and algorithm
# @param term_type Original term type ("s", "te", "ti", "t2")
# @param tensortype Default tensor type for s() terms
# @param is_gamm4 Whether using gamm4 algorithm
# @returns Symbol for the transformed term type
get_smooth_type <- function(term_type, tensortype, is_gamm4) {
  if (term_type %in% c("te", "ti") && is_gamm4) {
    return(quote(t2))
  }
  if (term_type == "s") {
    return(if (is_gamm4) quote(t2) else tensortype)
  }
  as.name(term_type)
}

# Determine dimension parameter for tensor smooth
# @param term The original term call
# @param term_type The term type string
# @returns Integer vector of dimensions for tensor product
get_smooth_dimensions <- function(term, term_type) {
  if (term_type == "s") {
    nvars <- if (!is.null(names(term))) {
      length(all.vars(term[names(term) == ""]))
    } else {
      length(all.vars(term))
    }
    c(nvars, 1)
  } else if ("d" %in% names(term)) {
    c(eval(term$d), 1)
  } else {
    rep(1, length(all.vars(term)) + 1)
  }
}

# Determine basis specification (bs) for tensor smooth
# @param term The original term call
# @param bs_yindex Basis spec for y-index marginal
# @returns Character vector of basis types for each marginal
get_smooth_bs <- function(term, bs_yindex, n_base_marginals) {
  term_bs <- if ("bs" %in% names(term)) eval(term$bs) else NULL
  yindex_bs <- bs_yindex$bs %||% "tp"

  if (!is.null(term_bs)) {
    if (length(term_bs) == 1 && n_base_marginals > 1) {
      term_bs <- rep(term_bs, n_base_marginals)
    }
    c(term_bs, yindex_bs)
  } else {
    c(rep("tp", n_base_marginals), yindex_bs)
  }
}

# Determine penalty order (m) for tensor smooth
# @param term The original term call
# @param bs_yindex Basis spec for y-index marginal
# @returns Penalty order specification
get_smooth_m <- function(term, bs_yindex, n_marginals) {
  yindex_m <- bs_yindex$m

  if ("m" %in% names(term)) {
    if ("m" %in% names(bs_yindex)) {
      warning("overriding bs.yindex for m in ", deparse(term))
    }
    term_m <- eval(term$m)
    term_m_len <- if (is.list(term_m)) length(term_m) else length(term_m)

    if (term_m_len == n_marginals) {
      return(term$m)
    }
    if (term_m_len == (n_marginals - 1)) {
      y_m <- yindex_m %||% NA
      if (is.list(term_m)) {
        return(c(term_m, list(y_m)))
      }
      if (length(y_m) == 1) {
        return(c(term_m, y_m))
      }
      return(c(as.list(term_m), list(y_m)))
    }
    return(term$m)
  }

  if (is.null(yindex_m)) {
    return(NULL)
  }
  if (n_marginals == 1) {
    return(yindex_m)
  }

  c(rep(list(NA), n_marginals - 1), list(yindex_m))
}

# Determine basis dimension (k) for tensor smooth
# @param term The original term call
# @param d Dimension vector from get_smooth_dimensions()
# @param bs_yindex Basis spec for y-index marginal
# @returns Integer vector of basis dimensions for each marginal
get_smooth_k <- function(term, d, bs_yindex) {
  yindex_k <- bs_yindex$k %||% 8
  k <- if ("k" %in% names(term)) {
    c(eval(term$k), yindex_k)
  } else {
    c(pmax(8, 5^head(d, -1)), yindex_k)
  }
  unlist(k)
}

#' Transform parametric terms to varying coefficient terms
#'
#' Converts scalar covariates to smooth varying coefficient terms.
#'
#' @param term The term (call or symbol).
#' @param yind_vec_name Name of y-index vector (character).
#' @param bs_yindex Basis spec for y-index (named list).
#' @returns Transformed term string.
#' @keywords internal
transform_par_term <- function(term, yind_vec_name, bs_yindex) {
  xnew <- as.call(c(quote(s), as.name(yind_vec_name), by = term, bs_yindex))
  safeDeparse(xnew)
}

#' Build the mgcv data frame from pffr components
#'
#' @param newfrml_env Environment containing the transformed variables.
#' @returns A data frame for mgcv fitting.
#' @keywords internal
build_mgcv_data <- function(newfrml_env) {
  list2df(as.list(newfrml_env))
}

#' Build the final mgcv formula from transformed terms
#'
#' @param response_name Response variable name (symbol or character).
#' @param intercept_string Intercept term string (or NULL if no intercept).
#' @param term_strings Character vector of transformed term strings.
#' @param has_intercept Logical, whether model has intercept.
#' @param formula_env Environment for formula.
#' @returns Formula object for mgcv.
#' @keywords internal
build_mgcv_formula <- function(
  response_name,
  intercept_string,
  term_strings,
  has_intercept,
  formula_env
) {
  frml_str <- if (has_intercept && !is.null(intercept_string)) {
    paste(response_name, "~", intercept_string)
  } else {
    paste(response_name, "~ 0")
  }

  if (length(term_strings) > 0) {
    frml_str <- paste(c(frml_str, term_strings), collapse = "+")
  }

  frml <- formula(frml_str)
  environment(frml) <- formula_env
  frml
}
