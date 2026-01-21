# =============================================================================
# Shared Test Fixtures for pffr Tests
# =============================================================================
#
# This file provides lazily-computed shared fixtures to reduce redundant
# data generation and model fitting across tests. Each fixture is computed
# once and cached for reuse within a test session.
#
# Usage: Call the getter function (e.g., get_basic_pffr_data()) to retrieve
# the fixture. The first call computes it; subsequent calls return the cache.

# Environment for caching fixtures
pffr_test_env <- new.env(parent = emptyenv())

# =============================================================================
# Basic Data Fixtures
# =============================================================================

#' Get basic pffr data with ff(X1) + xlin terms
#'
#' Generates a small dataset suitable for testing core pffr functionality.
#' Uses a fixed seed for reproducibility.
#'
#' @returns A pffr_simulate result with n=30, nxgrid=25, nygrid=30
get_basic_pffr_data <- function() {
  if (is.null(pffr_test_env$basic_data)) {
    set.seed(12345)
    pffr_test_env$basic_data <- pffr_simulate(
      Y ~ ff(X1, xind = s) + xlin,
      n = 30,
      nxgrid = 25,
      nygrid = 30,
      effects = list(X1 = "cosine", xlin = "dnorm"),
      SNR = 50
    )
  }
  pffr_test_env$basic_data
}

#' Get simple xlin-only data
#'
#' Generates a minimal dataset with only a linear varying coefficient term.
#' Useful for testing algorithm variants and basic functionality.
#'
#' @param n Number of observations (default: 30)
#' @param nygrid Number of y-grid points (default: 25)
#' @param SNR Signal-to-noise ratio (default: 50)
#' @returns A pffr_simulate result
get_xlin_data <- function(n = 30, nygrid = 25, SNR = 50) {
  cache_key <- paste0("xlin_data_", n, "_", nygrid, "_", SNR)
  if (is.null(pffr_test_env[[cache_key]])) {
    set.seed(23456)
    pffr_test_env[[cache_key]] <- pffr_simulate(
      Y ~ xlin,
      n = n,
      nygrid = nygrid,
      effects = list(xlin = "dnorm"),
      SNR = SNR
    )
  }
  pffr_test_env[[cache_key]]
}

#' Get multi-term data for comprehensive tests
#'
#' Generates data with multiple term types: ff, xlin, s(xsmoo), c(xconst).
#' Useful for testing term extraction, predict(type="terms"), etc.
#'
#' @returns A pffr_simulate result with n=35, nxgrid=25, nygrid=30
get_multiterm_data <- function() {
  if (is.null(pffr_test_env$multiterm_data)) {
    set.seed(34567)
    pffr_test_env$multiterm_data <- pffr_simulate(
      Y ~ ff(X1, xind = s) + xlin + s(xsmoo) + c(xconst),
      n = 35,
      nxgrid = 25,
      nygrid = 30,
      effects = list(
        X1 = "cosine",
        xlin = "dnorm",
        xsmoo = "sine",
        xconst = 2
      ),
      SNR = 50
    )
  }
  pffr_test_env$multiterm_data
}

#' Get historical ff (limits) data
#'
#' Generates data with ff term using limits="s<t" for historical effects.
#'
#' @returns A pffr_simulate result with historical ff term
get_historical_ff_data <- function() {
  if (is.null(pffr_test_env$historical_data)) {
    set.seed(45678)
    pffr_test_env$historical_data <- pffr_simulate(
      Y ~ ff(X1, xind = s),
      n = 30,
      nxgrid = 25,
      nygrid = 30,
      effects = list(X1 = "historical"),
      limits = function(s, t) s < t,
      SNR = 50
    )
  }
  pffr_test_env$historical_data
}

#' Get sparse/irregular data
#'
#' Generates data with missing observations for testing sparse functionality.
#'
#' @returns A pffr_simulate result with ydata component for sparse data
get_sparse_data <- function() {
  if (is.null(pffr_test_env$sparse_data)) {
    set.seed(56789)
    pffr_test_env$sparse_data <- pffr_simulate(
      Y ~ s(xsmoo),
      n = 30,
      nygrid = 25,
      effects = list(xsmoo = "sine"),
      propmissing = 0.5,
      SNR = 50
    )
  }
  pffr_test_env$sparse_data
}

#' Get data for AR tests
#'
#' Small dataset suitable for AR(1) correlation testing.
#'
#' @returns A pffr_simulate result with n=10, nygrid=12
get_ar_data <- function() {
  if (is.null(pffr_test_env$ar_data)) {
    set.seed(67890)
    # Use legacy scenario for intercept-only data
    pffr_test_env$ar_data <- suppressWarnings(
      pffr_simulate(scenario = "int", n = 10, nygrid = 12)
    )
  }
  pffr_test_env$ar_data
}

#' Get "all scenario" data (legacy compatibility)
#'
#' Generates data with all major term types for comprehensive testing.
#'
#' @returns A pffr_simulate result with ff, xlin, te, s, c terms
get_all_scenario_data <- function() {
  if (is.null(pffr_test_env$all_scenario_data)) {
    set.seed(78901)
    pffr_test_env$all_scenario_data <- suppressWarnings(
      pffr_simulate(scenario = "all", n = 30)
    )
  }
  pffr_test_env$all_scenario_data
}

# =============================================================================
# Fitted Model Fixtures
# =============================================================================

#' Get basic fitted pffr model
#'
#' Fits a pffr model with ff(X1) + xlin using the basic data fixture.
#' Cached for reuse across tests needing a fitted model.
#'
#' @returns A fitted pffr object
get_basic_pffr_model <- function() {
  if (is.null(pffr_test_env$basic_model)) {
    dat <- get_basic_pffr_data()
    s <- attr(dat, "xindex")
    t <- attr(dat, "yindex")
    pffr_test_env$basic_model <- pffr(
      Y ~ ff(X1, xind = s) + xlin,
      yind = t,
      data = dat
    )
  }
  pffr_test_env$basic_model
}

#' Get simple xlin-only fitted model
#'
#' Fits a minimal pffr model with only xlin term.
#' Cached for reuse across algorithm comparison tests.
#'
#' @returns A fitted pffr object
get_xlin_model <- function() {
  if (is.null(pffr_test_env$xlin_model)) {
    dat <- get_xlin_data()
    t <- attr(dat, "yindex")
    pffr_test_env$xlin_model <- pffr(
      Y ~ xlin,
      yind = t,
      data = dat
    )
  }
  pffr_test_env$xlin_model
}

#' Get multi-term fitted model
#'
#' Fits a model with multiple term types for comprehensive method tests.
#'
#' @returns A fitted pffr object with ff, xlin, s, and c terms
get_multiterm_model <- function() {
  if (is.null(pffr_test_env$multiterm_model)) {
    dat <- get_multiterm_data()
    s <- attr(dat, "xindex")
    t <- attr(dat, "yindex")
    pffr_test_env$multiterm_model <- pffr(
      Y ~ ff(X1, xind = s) + xlin + s(xsmoo) + c(xconst),
      yind = t,
      data = dat
    )
  }
  pffr_test_env$multiterm_model
}

#' Get sparse fitted model
#'
#' Fits a pffr model on sparse/irregular data.
#'
#' @returns A fitted pffr object for sparse data
get_sparse_model <- function() {
  if (is.null(pffr_test_env$sparse_model)) {
    dat <- get_sparse_data()
    t <- attr(dat, "yindex")
    pffr_test_env$sparse_model <- pffr(
      Y ~ s(xsmoo),
      data = dat$data,
      ydata = dat$ydata,
      yind = t
    )
  }
  pffr_test_env$sparse_model
}

#' Get "all scenario" fitted model
#'
#' Fits a model with all major term types (legacy compatibility).
#'
#' @returns A fitted pffr object with ff, xlin, te, s, c terms
get_all_scenario_model <- function() {
  if (is.null(pffr_test_env$all_scenario_model)) {
    dat <- get_all_scenario_data()
    s <- attr(dat, "xindex")
    t <- attr(dat, "yindex")
    pffr_test_env$all_scenario_model <- pffr(
      Y ~ ff(X1, xind = s) + xlin + c(te(xte1, xte2)) + s(xsmoo) + c(xconst),
      yind = t,
      data = dat
    )
  }
  pffr_test_env$all_scenario_model
}

#' Get gaulss fitted model
#'
#' Fits a pffr model with gaulss family for location-scale testing.
#'
#' @returns A fitted pffr object with gaulss family
get_gaulss_model <- function() {
  if (is.null(pffr_test_env$gaulss_model)) {
    dat <- get_xlin_data()
    t <- attr(dat, "yindex")
    pffr_test_env$gaulss_model <- pffr(
      Y ~ xlin,
      yind = t,
      data = dat,
      family = mgcv::gaulss()
    )
  }
  pffr_test_env$gaulss_model
}

#' Get pffrGLS fitted model
#'
#' Fits a pffrGLS model with AR(1) covariance structure.
#'
#' @returns A fitted pffr object from pffrGLS
get_gls_model <- function() {
  if (is.null(pffr_test_env$gls_model)) {
    dat <- get_xlin_data()
    t <- attr(dat, "yindex")
    nygrid <- length(t)
    rho <- 0.5
    hatSigma <- rho^abs(outer(1:nygrid, 1:nygrid, "-"))
    pffr_test_env$gls_model <- pffr_gls(
      Y ~ xlin,
      yind = t,
      data = dat,
      hatSigma = hatSigma
    )
  }
  pffr_test_env$gls_model
}

# =============================================================================
# Fixture Management
# =============================================================================

#' Clear all cached fixtures
#'
#' Removes all cached data and models from the fixture environment.
#' Useful for testing fixture behavior or freeing memory.
clear_pffr_fixtures <- function() {
  rm(list = ls(envir = pffr_test_env), envir = pffr_test_env)
  invisible(NULL)
}

#' List all cached fixtures
#'
#' @returns Character vector of cached fixture names
list_pffr_fixtures <- function() {
  ls(envir = pffr_test_env)
}
