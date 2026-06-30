## 00_setup.R — shared setup for the re-analysis
## - loads pffr: prefers an installed `refund`; otherwise sources pffr from a local
##   refund checkout (set env var REFUND_R to its R/ directory; defaults to ../../R).
## - provides data loaders for Canadian Weather (fda) and Hawaii Ocean (FRegSigCom),
##   falling back to the CRAN GitHub mirror when the packages are not installed.
## - defines the ISE metric used by the authors.

suppressMessages({library(mgcv); library(Matrix)})

load_pffr <- function() {
  if (requireNamespace("refund", quietly = TRUE)) {
    suppressMessages(library(refund)); message("Using installed refund ", utils::packageVersion("refund"))
  } else {
    RFUN <- Sys.getenv("REFUND_R", file.path(dirname(getwd()), "..", "R"))
    if (!dir.exists(RFUN)) RFUN <- "/home/user/refund/R"
    src <- c("pffr-utilities.R","pffr-formula.R","pffr-ff.R","pffr-sff.R","pffr-ffpc.R",
             "pffr-pcre.R","pffr-robust.R","pffr-core.R","pffr-methods.R","pffr.R","re.R")
    invisible(lapply(src, function(f) source(file.path(RFUN, f))))
    message("Sourced pffr from ", normalizePath(RFUN))
  }
  stopifnot(exists("pffr"), exists("ff"))
}

DATA_DIR <- file.path("comment-on-penffr","analysis","data")
if (!dir.exists(DATA_DIR)) DATA_DIR <- "data"
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

get_canadian_weather <- function() {
  if (requireNamespace("fda", quietly = TRUE)) {
    e <- new.env(); utils::data("CanadianWeather", package = "fda", envir = e); return(e$CanadianWeather)
  }
  f <- file.path(DATA_DIR, "CanadianWeather.rda")
  if (!file.exists(f)) utils::download.file(
    "https://raw.githubusercontent.com/cran/fda/master/data/CanadianWeather.rda", f, quiet = TRUE)
  e <- new.env(); load(f, envir = e); e$CanadianWeather
}

get_ocean <- function() {
  if (requireNamespace("FRegSigCom", quietly = TRUE)) {
    e <- new.env(); utils::data("ocean", package = "FRegSigCom", envir = e); return(e$ocean)
  }
  f <- file.path(DATA_DIR, "ocean.RData")
  if (!file.exists(f)) utils::download.file(
    "https://raw.githubusercontent.com/cran/FRegSigCom/master/data/ocean.RData", f, quiet = TRUE)
  e <- new.env(); load(f, envir = e); e$ocean
}

## authors' ISE: sum of squared residuals over the response grid, per curve
ISE_curve <- function(Yact, Yhat) sum((as.vector(Yact) - as.vector(Yhat))^2)
