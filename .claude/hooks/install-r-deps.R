#!/usr/bin/env Rscript
# Installs the refund package's dependencies (plus fastFMM) for Claude Code web
# sessions. CRAN and its mirrors are blocked by the web sandbox network policy,
# so we install precompiled binaries from Ubuntu's r-cran-* apt packages where
# available and fall back to building from the GitHub CRAN mirror (github.com/cran)
# for the handful that apt does not ship. Idempotent: already-installed packages
# are skipped.

options(timeout = 600)

PROJECT  <- Sys.getenv("CLAUDE_PROJECT_DIR", ".")
base_pkgs <- rownames(installed.packages(priority = "base"))
gh_dir   <- "/tmp/gh_pkgs"
dir.create(gh_dir, showWarnings = FALSE, recursive = TRUE)
in_progress <- character(0)

is_installed <- function(p) p %in% rownames(installed.packages())

apt_try <- function(p) {
  lc <- tolower(p)
  if (system(sprintf("apt-cache show r-cran-%s >/dev/null 2>&1", lc)) != 0) return(FALSE)
  message(sprintf("  apt: r-cran-%s", lc))
  system(sprintf("DEBIAN_FRONTEND=noninteractive apt-get install -y r-cran-%s >/dev/null 2>&1", lc))
  is_installed(p)
}

fetch_desc <- function(p) {
  for (br in c("master", "main")) {
    dest <- tempfile()
    ok <- tryCatch(
      download.file(sprintf("https://raw.githubusercontent.com/cran/%s/%s/DESCRIPTION", p, br),
                    dest, quiet = TRUE) == 0,
      error = function(e) FALSE)
    if (ok && file.info(dest)$size > 0) return(dest)
  }
  NA_character_
}

parse_deps <- function(descfile) {
  dcf <- read.dcf(descfile)
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(dcf))
  parts <- unlist(strsplit(paste(dcf[1, fields], collapse = ","), ","))
  parts <- trimws(gsub("\\(.*?\\)", "", parts))
  parts <- parts[parts != "" & parts != "R"]
  setdiff(unique(parts), base_pkgs)
}

gh_install <- function(p) {
  for (br in c("master", "main")) {
    dest <- file.path(gh_dir, paste0(p, ".tar.gz"))
    ok <- tryCatch(
      download.file(sprintf("https://github.com/cran/%s/archive/refs/heads/%s.tar.gz", p, br),
                    dest, quiet = TRUE) == 0,
      error = function(e) FALSE)
    if (!ok || file.info(dest)$size < 200) next
    untar(dest, exdir = gh_dir)
    srcdir <- file.path(gh_dir, paste0(p, "-", br))
    if (!dir.exists(srcdir)) next
    system(sprintf("R CMD INSTALL --no-test-load '%s' >/tmp/inst_%s.log 2>&1", srcdir, p))
    return(is_installed(p))
  }
  FALSE
}

ensure <- function(p) {
  if (p %in% base_pkgs || is_installed(p) || p %in% in_progress) return(invisible(TRUE))
  in_progress[[length(in_progress) + 1]] <<- p
  message("ensure: ", p)
  if (apt_try(p)) return(invisible(TRUE))
  desc <- fetch_desc(p)
  if (is.na(desc)) stop(sprintf("Cannot find DESCRIPTION for '%s' on GitHub CRAN mirror", p))
  for (d in parse_deps(desc)) ensure(d)
  if (!gh_install(p)) {
    log <- sprintf("/tmp/inst_%s.log", p)
    msg <- if (file.exists(log)) paste(tail(readLines(log, warn = FALSE), 15), collapse = "\n") else ""
    stop(sprintf("Failed to install '%s'\n%s", p, msg))
  }
  invisible(TRUE)
}

# ggdist (a transitive dep of fastFMM) needs ggplot2 >= 3.5.0, but apt only ships
# 3.4.4. Pin 3.5.2 from the GitHub mirror into /usr/local (higher libPath priority).
ensure_ggplot2 <- function() {
  ensure("ggplot2")
  if (utils::packageVersion("ggplot2") < "3.5.0") {
    message("ensure: ggplot2 3.5.2 (apt version too old for ggdist)")
    dest <- file.path(gh_dir, "ggplot2-3.5.2.tar.gz")
    download.file("https://github.com/cran/ggplot2/archive/refs/tags/3.5.2.tar.gz", dest, quiet = TRUE)
    untar(dest, exdir = gh_dir)
    system(sprintf("R CMD INSTALL --no-test-load '%s' >/tmp/inst_ggplot2.log 2>&1",
                   file.path(gh_dir, "ggplot2-3.5.2")))
  }
}

# Targets: refund's own declared dependencies (from the local DESCRIPTION) plus fastFMM.
dcf <- read.dcf(file.path(PROJECT, "DESCRIPTION"))
fields <- intersect(c("Depends", "Imports", "Suggests"), colnames(dcf))
targets <- trimws(gsub("\\(.*?\\)", "", unlist(strsplit(paste(dcf[1, fields], collapse = ","), ","))))
targets <- setdiff(unique(targets[targets != "" & targets != "R"]), base_pkgs)
targets <- unique(c(targets, "fastFMM"))

ensure_ggplot2()
for (t in targets) ensure(t)

# Install refund itself from the local source so library(refund) and the test
# suite work out of the box.
system(sprintf("R CMD INSTALL --no-test-load '%s' >/tmp/inst_refund.log 2>&1", PROJECT))

cat("\n--- R dependency install summary ---\n")
cat("fastFMM installed:", is_installed("fastFMM"), "\n")
cat("refund  installed:", is_installed("refund"), "\n")
missing <- targets[!vapply(targets, is_installed, logical(1))]
if (length(missing)) {
  cat("MISSING:", paste(missing, collapse = ", "), "\n")
  quit(status = 1)
}
cat("All dependencies present.\n")
