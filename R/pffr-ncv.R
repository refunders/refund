# NCV neighbourhoods and installed-mgcv behavioural check ---------------------

pffr_ncv_nei <- function(cluster_id) {
  if (
    !is.atomic(cluster_id) ||
      !is.null(dim(cluster_id)) ||
      length(cluster_id) < 2L ||
      anyNA(cluster_id)
  ) {
    stop(
      "NCV cluster ids must be a nonmissing vector with at least two rows.",
      call. = FALSE
    )
  }
  groups <- split(
    seq_along(cluster_id),
    factor(cluster_id, levels = unique(cluster_id))
  )
  if (length(groups) < 2L) {
    stop("NCV requires at least two blocks.", call. = FALSE)
  }
  indices <- as.integer(unlist(groups, use.names = FALSE))
  ends <- as.integer(cumsum(lengths(groups)))
  # mgcv < 1.9-4 uses k/m/i/mi; newer versions use a/ma/d/md.
  list(
    a = indices,
    ma = ends,
    d = indices,
    md = ends,
    k = indices,
    m = ends,
    i = indices,
    mi = ends,
    jackknife = FALSE
  )
}

.pffr_ncv_cache <- new.env(parent = emptyenv())

pffr_ncv_probe <- function() {
  # Deterministic data: neither create nor alter .Random.seed.
  j <- seq_len(160L)
  dat <- data.frame(x = (j * 37 %% 163) / 163, id = rep(1:20, each = 8))
  dat$y <- sin(5 * dat$x) + sin(dat$id * 1.7) + cos(j * 2.3) / 2
  fit <- function(nei)
    mgcv::gam(
      y ~ s(x, bs = "ps", k = 8),
      data = dat,
      method = "NCV",
      nei = nei,
      control = mgcv::gam.control(nthreads = 1L, ncv.threads = 1L)
    )$sp
  list(blocked = fit(pffr_ncv_nei(dat$id)), point = fit(pffr_ncv_nei(j)))
}

pffr_ncv_check_blocks <- function() {
  if (isTRUE(.pffr_ncv_cache$checked)) return(invisible(TRUE))
  if (utils::packageVersion("mgcv") < "1.9.0") {
    stop("pffr NCV requires mgcv >= 1.9.0.", call. = FALSE)
  }
  result <- pffr_ncv_probe()
  if (
    !length(result$blocked) ||
      !length(result$point) ||
      any(!is.finite(c(result$blocked, result$point))) ||
      isTRUE(all.equal(result$blocked, result$point, tolerance = 1e-6))
  ) {
    stop(
      "Installed mgcv ignored the NCV block structure: blocked and pointwise fits do not differ. Update mgcv before using pffr NCV.",
      call. = FALSE
    )
  }
  .pffr_ncv_cache$checked <- TRUE
  invisible(TRUE)
}

pffr_ncv_setup <- function(prep, cluster, blocks, envir) {
  if (any(c("G", "fit") %in% names(prep$dots))) {
    stop("pffr NCV does not support supplying `G` or `fit`.", call. = FALSE)
  }
  pffr_ncv_check_blocks()
  setup_call <- prep$new_call
  setup_call$nei <- NULL
  setup_call$fit <- FALSE
  setup <- eval(setup_call, envir)
  if (inherits(attr(setup$mf, "na.action"), "exclude")) {
    stop(
      "pffr NCV requires na.action = na.omit when rows are missing; na.exclude pads working rows and is incompatible with cluster inference.",
      call. = FALSE
    )
  }
  # gam(fit = FALSE) constructs the real model frame, including na.action,
  # matrix covariates, weights and offsets. Only response omissions are
  # represented in pffr's sandwich/prediction metadata; reject other omissions.
  expected <- seq_len(nrow(prep$pffr_data))
  if (length(prep$missing_indices)) expected <- expected[-prep$missing_indices]
  if (!identical(rownames(setup$mf), rownames(prep$pffr_data)[expected])) {
    stop(
      "pffr NCV cannot align the model frame: rows were removed or reordered beyond missing response values. Remove missing covariates/weights/offsets before fitting; use the default na.action.",
      call. = FALSE
    )
  }
  cid <- build_cluster_id(prep, cluster)
  if (length(cid) != nrow(setup$mf)) {
    stop(
      "pffr NCV cluster ids do not match the fitted model frame.",
      call. = FALSE
    )
  }
  custom <- "nei" %in% names(prep$dots)
  if (custom) {
    message(
      "pffr NCV: using user-supplied `nei` unchanged (indices refer to the retained model-frame rows)."
    )
    nei <- prep$dots$nei
    blocks <- "user"
    n_blocks <- if (is.null(nei)) nrow(setup$mf) else length(nei$ma %||% nei$m)
  } else {
    ids <- if (blocks == "point") seq_along(cid) else cid
    nei <- pffr_ncv_nei(ids)
    n_blocks <- length(nei$ma)
  }
  # Passing G avoids mgcv's nanei() translating already-retained row indices
  # a second time. The original basis and fixed penalties are reused exactly.
  setup$cl$fit <- NULL
  setup$cl$nei <- nei
  list(
    setup = setup,
    info = list(blocks = blocks, n_blocks = n_blocks, nei = nei)
  )
}
