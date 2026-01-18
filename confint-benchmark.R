# confint-benchmark.R
# Benchmark CI coverage/width from coef.pffr() under different error structures.
# Methods: pffr, pffr_sandwich, pffr_gls_est, pffr_ar
# Families: gaussian, gaulss, scat

`%||%` <- function(x, y) if (is.null(x)) y else x

# Setup ----

load_refund_dev <- function(pkg_dir = ".") {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
    devtools::load_all(pkg_dir, quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(refund))
  }
}

require_pkgs <- function() {
  load_refund_dev(".")
  stopifnot(requireNamespace("mgcv", quietly = TRUE),
            requireNamespace("mvtnorm", quietly = TRUE))
}

# Truth functions ----

make_bs_spec <- function(x, df) {
  B <- splines::bs(x, df = df, intercept = TRUE)
  list(knots = attr(B, "knots"), bknots = attr(B, "Boundary.knots"),
       degree = attr(B, "degree"), df = df)
}

eval_bs <- function(x, spec) {
  splines::bs(x, df = spec$df, knots = spec$knots,
              Boundary.knots = spec$bknots, degree = spec$degree, intercept = TRUE)
}

make_truth_functions <- function(s_grid, t_grid, z_range, wiggliness,
                                  seed = NULL, k_te = 5L, k_uni = 7L) {
  if (!is.null(seed)) set.seed(seed)

  rw_coef <- function(k, sd) { c <- cumsum(rnorm(k, sd = sd)); c - mean(c) }
  rw_coef_mat <- function(k, sd) {
    m <- matrix(rnorm(k * k, sd = sd), k)
    m <- apply(m, 2, cumsum); m <- t(apply(m, 1, cumsum)); m - mean(m)
  }

  s_spec <- make_bs_spec(s_grid, k_te)
  t_te_spec <- make_bs_spec(t_grid, k_te)
  t_uni_spec <- make_bs_spec(t_grid, k_uni)
  z_spec <- make_bs_spec(z_range, k_te)

  cf_ff <- rw_coef_mat(k_te, wiggliness)
  cf_con <- rw_coef(k_uni, wiggliness)
  cf_lin <- rw_coef(k_uni, wiggliness)
  cf_mu <- rw_coef(k_uni, wiggliness)
  cf_smooth <- rw_coef_mat(k_te, wiggliness)

  list(
    beta_ff = function(s, t) drop(eval_bs(s, s_spec) %*% cf_ff %*% t(eval_bs(t, t_te_spec))),
    beta_concurrent = function(t) drop(eval_bs(t, t_uni_spec) %*% cf_con),
    beta_linear = function(t) drop(eval_bs(t, t_uni_spec) %*% cf_lin),
    mu_t = function(t) drop(eval_bs(t, t_uni_spec) %*% cf_mu),
    f_smooth = function(z, t) drop(eval_bs(z, z_spec) %*% cf_smooth %*% t(eval_bs(t, t_te_spec)))
  )
}

# Integration weights ----

ff_weights <- function(xind) {
  nx <- length(xind)
  if (nx < 2) return(1)
  ((xind[nx] - xind[1]) / (nx - 1)) / 3 * c(1, rep(c(4, 2), length = nx - 2), 1)
}

center_beta_ff <- function(beta_mat, s_grid) {
  w <- ff_weights(s_grid)
  sweep(beta_mat, 2, colSums(beta_mat * w) / sum(w), "-")
}

# Error structure ----

make_error_cov <- function(t_grid, corr_type, hetero_type) {
  ny <- length(t_grid)

  corr <- switch(corr_type,
    iid = diag(ny),
    ar1 = 0.4^abs(outer(1:ny, 1:ny, "-")),
    gauss = exp(-abs(outer(t_grid, t_grid, "-"))^2 / (2 * 0.15^2)),
    fourier = {
      B <- cbind(sin(2*pi*t_grid/0.33), cos(2*pi*t_grid/0.33),
                 sin(2*pi*t_grid/0.17), cos(2*pi*t_grid/0.17))
      B <- scale(B, center = TRUE, scale = FALSE)
      cov2cor(tcrossprod(B) + diag(0.25, ny))
    }
  )

  sigma_t <- switch(hetero_type,
    none = rep(1, ny),
    bump = 1 + exp(-0.5 * ((t_grid - 0.7) / 0.10)^2)
  )

  cov_mat <- diag(sigma_t) %*% corr %*% diag(sigma_t)
  cov_mat + diag(1e-8, ny)  # ensure PD
}

sample_errors <- function(n, cov_mat, error_dist) {
  ny <- nrow(cov_mat)
  if (error_dist == "gaussian") {
    mvtnorm::rmvnorm(n, rep(0, ny), cov_mat)
  } else {
    mvtnorm::rmvt(n, cov_mat * 4/6, df = 6, delta = rep(0, ny))
  }
}

# Data simulation ----

simulate_dataset <- function(n, nxgrid, nygrid, snr, wiggliness, error_dist,
                              corr_type, hetero_type, seed, terms) {
  set.seed(seed)

  s_grid <- seq(0, 1, length.out = nxgrid)
  t_grid <- seq(0, 1, length.out = nygrid)
  z_range <- c(-2, 2)

  truth_funs <- make_truth_functions(s_grid, t_grid, z_range, wiggliness, seed + 1e5L)

  # Generate covariates (centered)
  gen_X <- function() { m <- matrix(rnorm(n * nxgrid), n); sweep(m, 2, colMeans(m)) }
  gen_z <- function() { z <- rnorm(n); z - mean(z) }

  dat <- data.frame(id = 1:n)
  eta <- matrix(truth_funs$mu_t(t_grid), n, nygrid, byrow = TRUE)  # intercept
  beta <- list()

  if ("ff" %in% terms) {
    dat$X1 <- I(gen_X())
    beta$ff <- center_beta_ff(truth_funs$beta_ff(s_grid, t_grid), s_grid)
    w <- ff_weights(s_grid)
    eta <- eta + (dat$X1 * rep(w, each = n)) %*% beta$ff
  }
  if ("concurrent" %in% terms) {
    dat$Xc <- I(gen_X())
    beta$concurrent <- truth_funs$beta_concurrent(t_grid)
    eta <- eta + sweep(dat$Xc, 2, beta$concurrent, "*")
  }
  if ("linear" %in% terms) {
    dat$zlin <- gen_z()
    beta$linear <- truth_funs$beta_linear(t_grid)
    eta <- eta + outer(dat$zlin, beta$linear)
  }
  if ("smooth" %in% terms) {
    dat$zsmoo <- gen_z() * diff(z_range)/4  # scale to z_range
    beta$smooth <- truth_funs$f_smooth  # function, not matrix
    sm_contrib <- t(sapply(dat$zsmoo, function(z) truth_funs$f_smooth(z, t_grid)))
    sm_contrib <- sweep(sm_contrib, 2, colMeans(sm_contrib))  # center per t
    eta <- eta + sm_contrib
  }

  # Errors scaled to SNR
  cov_mat <- make_error_cov(t_grid, corr_type, hetero_type)
  eps <- sample_errors(n, cov_mat, error_dist)
  scale_fac <- sqrt(var(as.vector(eta)) / (snr * var(as.vector(eps))))
  eps <- eps * scale_fac

  dat$Y <- I(eta + eps)

  list(
    data = dat, s_grid = s_grid, t_grid = t_grid,
    truth = list(eta = eta, beta = beta, funs = truth_funs, terms = terms),
    rho = if (corr_type == "ar1") 0.4 else NA
  )
}

# Fitting helpers ----

get_family <- function(name) {
 switch(name, gaussian = gaussian(), gaulss = mgcv::gaulss(), scat = mgcv::scat())
}

estimate_rho <- function(fit, dat) {
  r <- as.matrix(dat$Y) - fitted(fit)
  ny <- ncol(r)
  rho <- cor(as.vector(r[, -ny]), as.vector(r[, -1]))
  max(-0.99, min(0.99, ifelse(is.finite(rho), rho, 0)))
}

estimate_hatSigma <- function(fit, dat) {
  r <- as.matrix(dat$Y) - fitted(fit)
  S <- cov(r); S <- 0.5 * (S + t(S))
  S + diag(1e-8 * mean(diag(S)), nrow(S))
}

# Scoring ----

score_fit <- function(fit, truth, s_grid, t_grid, use_sandwich, alpha, terms) {
  coefs <- suppressMessages(coef(fit, sandwich = use_sandwich,
                                  n1 = length(t_grid), n2 = length(s_grid), n3 = 20))
  z <- qnorm(1 - alpha/2)
  results <- list()

  term_map <- c(ff = "X1", concurrent = "Xc", linear = "zlin", smooth = "zsmoo")

  for (term in terms) {
    sm_idx <- grep(term_map[term], names(coefs$smterms), fixed = TRUE)
    if (length(sm_idx) == 0) next
    sm <- coefs$smterms[[sm_idx[1]]]
    est <- sm$coef[, "value"]
    se <- sm$coef[, "se"]

    # Evaluate truth on the EXACT grid from coef()
    if (term == "ff") {
      eval_s <- unique(sm$coef[, 1])
      eval_t <- unique(sm$coef[, 2])
      truth_mat <- truth$funs$beta_ff(eval_s, eval_t)
      truth_mat <- center_beta_ff(truth_mat, eval_s)  # center on eval grid!
      truth_vec <- as.vector(truth_mat)
    } else if (term == "concurrent") {
      truth_vec <- truth$funs$beta_concurrent(sm$coef[, 1])
    } else if (term == "linear") {
      truth_vec <- truth$funs$beta_linear(sm$coef[, 1])
    } else if (term == "smooth") {
      eval_z <- unique(sm$coef[, 1])
      eval_t <- unique(sm$coef[, 2])
      truth_mat <- outer(eval_z, eval_t, function(z, t) truth$funs$f_smooth(z, t))
      truth_mat <- sweep(truth_mat, 2, colMeans(truth_mat))  # center per t
      truth_vec <- as.vector(truth_mat)
    } else next

    covered <- (truth_vec >= est - z*se) & (truth_vec <= est + z*se)
    results[[term]] <- data.frame(
      term_type = term, coverage = mean(covered),
      width = mean(2*z*se), rmse = sqrt(mean((est - truth_vec)^2))
    )
  }

  # E(y) metrics
  pred <- predict(fit, se.fit = TRUE)
  if (is.list(pred)) { fv <- pred$fit; se_fv <- pred$se.fit }
  else { fv <- pred; se_fv <- NA }
  truth_eta <- as.vector(truth$eta)
  fv <- as.vector(fv); se_fv <- as.vector(se_fv)

  if (length(fv) == length(truth_eta) && !any(is.na(se_fv))) {
    cov_ey <- (truth_eta >= fv - z*se_fv) & (truth_eta <= fv + z*se_fv)
    results$Ey <- data.frame(
      term_type = "Ey", coverage = mean(cov_ey),
      width = mean(2*z*se_fv), rmse = sqrt(mean((fv - truth_eta)^2))
    )
  }

  do.call(rbind, results)
}

# Core: fit all methods for one dataset + family ----

fit_and_score <- function(sim, family, methods, alpha, fit_args) {
  dat <- sim$data
  t_grid <- sim$t_grid
  s_grid <- sim$s_grid
  truth <- sim$truth
  terms <- truth$terms
  is_gauss <- (family == "gaussian")

  # Build formula
  fterms <- c(
    if ("ff" %in% terms) "ff(X1, xind = s_grid)" else NULL,
    if ("concurrent" %in% terms) "Xc" else NULL,
    if ("linear" %in% terms) "zlin" else NULL,
    if ("smooth" %in% terms) "s(zsmoo, k = 10)" else NULL
  )
  formula <- as.formula(paste("Y ~", paste(fterms, collapse = " + ")))
  fam_obj <- get_family(family)

  # === FIT ALL MODELS (not in a loop!) ===
  fits <- list()
  times <- list()

  # 1. Base pffr
  t0 <- Sys.time()
  fits$pffr <- tryCatch(
    do.call(pffr, c(list(formula = formula, yind = t_grid, data = dat, family = fam_obj), fit_args)),
    error = function(e) NULL
  )
  times$pffr <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(fits$pffr)) return(NULL)

  # 2. pffr_sandwich = same fit object
  fits$pffr_sandwich <- fits$pffr
  times$pffr_sandwich <- times$pffr

  # 3. Estimate rho/hatSigma from base fit
  rho_est <- estimate_rho(fits$pffr, dat)
  hatSigma_est <- estimate_hatSigma(fits$pffr, dat)

  # 4. pffr_gls_est (gaussian only)
  if (is_gauss && "pffr_gls_est" %in% methods) {
    t0 <- Sys.time()
    fits$pffr_gls_est <- tryCatch(
      do.call(pffr_gls, c(list(formula = formula, yind = t_grid, data = dat,
                               hatSigma = hatSigma_est), fit_args)),
      error = function(e) NULL
    )
    times$pffr_gls_est <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }

  # 5. pffr_ar (gaussian only)
  if (is_gauss && "pffr_ar" %in% methods) {
    t0 <- Sys.time()
    fits$pffr_ar <- tryCatch(
      do.call(pffr, c(list(formula = formula, yind = t_grid, data = dat,
                           family = fam_obj, algorithm = "bam", rho = rho_est), fit_args)),
      error = function(e) NULL
    )
    times$pffr_ar <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }

  # === SCORE ALL MODELS ===
  results <- list()
  for (method in intersect(methods, names(fits))) {
    if (is.null(fits[[method]])) next
    use_sw <- (method == "pffr_sandwich")

    t0 <- Sys.time()
    scored <- tryCatch(
      score_fit(fits[[method]], truth, s_grid, t_grid, use_sw, alpha, terms),
      error = function(e) NULL
    )
    time_score <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (!is.null(scored)) {
      scored$method <- method
      scored$family <- family
      scored$time_fit_sec <- times[[method]]
      scored$time_score_sec <- time_score
      results[[method]] <- scored
    }
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# Run one (dgp, rep) ----

get_families <- function(dgp) {
  fams <- "gaussian"
  if (dgp$hetero_type != "none") fams <- c(fams, "gaulss", "scat")
  else if (dgp$error_dist == "t6") fams <- c(fams, "scat")
  fams
}

run_one <- function(dgp, rep_id, seed, methods, alpha, terms, fit_args) {
  sim <- simulate_dataset(
    n = dgp$n, nxgrid = dgp$nxgrid, nygrid = dgp$nygrid,
    snr = dgp$snr, wiggliness = dgp$wiggliness,
    error_dist = dgp$error_dist, corr_type = dgp$corr_type,
    hetero_type = dgp$hetero_type, seed = seed, terms = terms
  )

  results <- list()
  for (fam in get_families(dgp)) {
    res <- fit_and_score(sim, fam, methods, alpha, fit_args)
    if (!is.null(res)) results[[fam]] <- res
  }
  if (length(results) == 0) return(NULL)

  out <- do.call(rbind, results)
  out$dgp_id <- dgp$dgp_id
  out$rep_id <- rep_id
  out$snr <- dgp$snr
  out$wiggliness <- dgp$wiggliness
  out$error_dist <- dgp$error_dist
  out$corr_type <- dgp$corr_type
  out$hetero_type <- dgp$hetero_type
  rownames(out) <- NULL
  out
}

# Grids ----

make_dgp_grid <- function() {
  expand.grid(
    n = 80L, nxgrid = 35L, nygrid = 45L,
    snr = c(2, 8, 25), wiggliness = c(2, 20),
    error_dist = c("gaussian", "t6"),
    corr_type = c("iid", "ar1", "gauss", "fourier"),
    hetero_type = c("none", "bump"),
    stringsAsFactors = FALSE
  )
}

make_test_grid <- function() {
  data.frame(
    n = 30L, nxgrid = 15L, nygrid = 20L, snr = 10, wiggliness = 2,
    error_dist = c("gaussian", "gaussian", "gaussian", "t6", "t6"),
    corr_type = c("iid", "ar1", "iid", "iid", "ar1"),
    hetero_type = c("none", "none", "bump", "none", "bump"),
    stringsAsFactors = FALSE
  )
}

# Main benchmark ----

run_benchmark <- function(dgp_grid, n_rep = 10L,
                          methods = c("pffr", "pffr_sandwich", "pffr_gls_est", "pffr_ar"),
                          alpha = 0.10, seed = 2024L, n_workers = 1L,
                          output_dir = NULL, terms = c("ff", "linear"),
                          fit_args = list()) {
  require_pkgs()
  dgp_grid$dgp_id <- seq_len(nrow(dgp_grid))

  tasks <- expand.grid(dgp_idx = seq_len(nrow(dgp_grid)), rep_id = seq_len(n_rep))
  tasks$seed <- seed + 1000L * tasks$dgp_idx + tasks$rep_id
  n_tasks <- nrow(tasks)

  cat(sprintf("\n=== BENCHMARK: %d DGPs x %d reps = %d tasks ===\n",
              nrow(dgp_grid), n_rep, n_tasks))
  cat(sprintf("Methods: %s | Terms: %s | Workers: %d\n\n",
              paste(methods, collapse=", "), paste(terms, collapse=", "), n_workers))

  if (!is.null(output_dir) && !dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  do_task <- function(i) {
    dgp <- dgp_grid[tasks$dgp_idx[i], , drop = FALSE]
    run_one(dgp, tasks$rep_id[i], tasks$seed[i], methods, alpha, terms, fit_args)
  }

  start <- Sys.time()
  all_res <- vector("list", n_tasks)

  if (n_workers > 1 && .Platform$OS.type != "windows") {
    # Parallel in chunks for incremental logging
    chunks <- split(seq_len(n_tasks), ceiling(seq_len(n_tasks) / (n_workers * 2)))
    for (ch in chunks) {
      ch_res <- parallel::mclapply(ch, do_task, mc.cores = n_workers)
      for (j in seq_along(ch)) {
        i <- ch[j]
        all_res[[i]] <- ch_res[[j]]
        if (!is.null(ch_res[[j]])) {
          dgp <- dgp_grid[tasks$dgp_idx[i], ]
          cov_str <- paste(sapply(split(ch_res[[j]], ch_res[[j]]$method),
                                  function(m) sprintf("%s:%.2f", m$method[1], mean(m$coverage))),
                           collapse=" | ")
          cat(sprintf("[%d/%d] dgp=%d rep=%d %s/%s/%s | %s\n", i, n_tasks,
                      dgp$dgp_id, tasks$rep_id[i], dgp$corr_type, dgp$hetero_type,
                      dgp$error_dist, cov_str))
        }
      }
      # Incremental save
      if (!is.null(output_dir)) {
        partial <- do.call(rbind, all_res[!sapply(all_res, is.null)])
        if (!is.null(partial)) saveRDS(partial, file.path(output_dir, "results_partial.rds"))
      }
    }
  } else {
    for (i in seq_len(n_tasks)) {
      all_res[[i]] <- tryCatch(do_task(i), error = function(e) NULL)
      if (!is.null(all_res[[i]])) {
        dgp <- dgp_grid[tasks$dgp_idx[i], ]
        cov_str <- paste(sapply(split(all_res[[i]], all_res[[i]]$method),
                                function(m) sprintf("%s:%.2f", m$method[1], mean(m$coverage))),
                         collapse=" | ")
        cat(sprintf("[%d/%d] dgp=%d rep=%d %s/%s/%s | %s\n", i, n_tasks,
                    dgp$dgp_id, tasks$rep_id[i], dgp$corr_type, dgp$hetero_type,
                    dgp$error_dist, cov_str))
      }
      if (!is.null(output_dir) && i %% 5 == 0) {
        partial <- do.call(rbind, all_res[!sapply(all_res, is.null)])
        if (!is.null(partial)) saveRDS(partial, file.path(output_dir, "results_partial.rds"))
      }
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "mins"))
  cat(sprintf("\n=== COMPLETE: %.1f min ===\n", elapsed))

  results <- do.call(rbind, all_res[!sapply(all_res, is.null)])
  rownames(results) <- NULL

  if (!is.null(output_dir) && !is.null(results)) {
    saveRDS(results, file.path(output_dir, "results_final.rds"))
    write.csv(results, file.path(output_dir, "results_final.csv"), row.names = FALSE)
    cat("Saved to:", output_dir, "\n")
  }

  results
}

# Summary ----

summarize_results <- function(res) {
  aggregate(cbind(coverage, width, rmse, time_fit_sec, time_score_sec) ~
              method + family + term_type + snr + wiggliness +
              error_dist + corr_type + hetero_type,
            data = res, FUN = mean)
}

# Auto-run test if sourced ----

if (sys.nframe() == 0) {
  require_pkgs()
  res <- run_benchmark(
    make_test_grid(), n_rep = 1L,
    methods = c("pffr", "pffr_sandwich", "pffr_gls_est", "pffr_ar"),
    terms = c("ff", "linear"),
    output_dir = file.path(getwd(), "ci-benchmark")
  )
  if (!is.null(res)) {
    cat("\n=== SUMMARY ===\n")
    print(summarize_results(res)[, c("method","family","term_type","corr_type",
                                      "hetero_type","coverage","time_fit_sec")])
  }
}
