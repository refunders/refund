# Timing benchmark for the fixed-fit cluster-influence degrees-of-freedom
# computation (`pffr_influence_core()` / `pffr_influence_df()`).
#
# Measures, on a realistic Gaussian pffr() fit (ff(X1) + xlin, n_y = 60,
# ff basis 12 x 12, bs.yindex k = 12):
#   (a) influence-core construction  (pffr_influence(), cold cache)
#   (b) coef(sandwich = "cl2", cl2_adjustment = "exact",
#           crit = "satterthwaite", ci = "pointwise")  -- coefficient grids
#   (c) df on an E(Y)-type contrast set: the working-design rows returned by
#       predict(type = "lpmatrix"), through pffr_df_from_context()
#
# Run single-threaded; the timings are dominated by small dense BLAS calls.
#
# Usage:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript inst/benchmarks/df-timing.R \
#     <label> <outfile.csv> [G1,G2] [lp_rows] [reps]
#
# `lp_rows` subsamples the lpmatrix contrast set (NA = all rows); the per-row
# cost is data independent, so the full-grid time is extrapolated linearly and
# reported alongside the measured time.

args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args) >= 1L) args[[1L]] else "current"
outfile <- if (length(args) >= 2L) args[[2L]] else "df-timing.csv"
Gs <- if (length(args) >= 3L) {
  as.integer(strsplit(args[[3L]], ",", fixed = TRUE)[[1L]])
} else {
  c(100L, 200L)
}
lp_rows <- if (length(args) >= 4L) as.integer(args[[4L]]) else NA_integer_
reps <- if (length(args) >= 5L) as.integer(args[[5L]]) else 3L

suppressPackageStartupMessages(library(refund))

.cnt <- new.env(parent = emptyenv())

med_time <- function(expr, reps = 3L) {
  call <- substitute(expr)
  env <- parent.frame()
  el <- vapply(
    seq_len(reps),
    function(i) system.time(eval(call, env))[["elapsed"]],
    numeric(1)
  )
  list(median = stats::median(el), min = min(el), max = max(el), reps = reps)
}

clear_influence_cache <- function(fit) {
  cache <- fit$pffr$Vsandwich_cache
  if (is.environment(cache)) {
    rm(list = ls(cache, all.names = TRUE), envir = cache)
  }
  fit$pffr$Vsandwich <- NULL
  fit
}

count_df_contrasts <- function(expr) {
  call <- substitute(expr)
  env <- parent.frame()
  assign("n", 0L, envir = .cnt)
  assign("calls", 0L, envir = .cnt)
  suppressMessages(trace(
    "pffr_influence_df",
    where = asNamespace("refund"),
    tracer = quote({
      assign("n", get("n", .cnt) + nrow(as.matrix(Xp)), .cnt)
      assign("calls", get("calls", .cnt) + 1L, .cnt)
    }),
    print = FALSE
  ))
  on.exit(suppressMessages(untrace(
    "pffr_influence_df",
    where = asNamespace("refund")
  )))
  eval(call, env)
  list(contrasts = get("n", .cnt), calls = get("calls", .cnt))
}

make_fit <- function(G) {
  dat <- pffr_simulate(
    Y ~ ff(X1) + xlin,
    n = G,
    nxgrid = 40,
    nygrid = 60,
    effects = list(X1 = "cosine", xlin = "dnorm"),
    seed = 4711L
  )
  pffr(
    Y ~ ff(X1, splinepars = list(bs = "ps", k = c(12, 12))) + xlin,
    yind = attr(dat, "yindex"),
    data = dat,
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    sandwich = "cl2",
    cl2_adjustment = "exact"
  )
}

rows <- list()
add_row <- function(...) rows[[length(rows) + 1L]] <<- data.frame(...)

for (G in Gs) {
  message("=== G = ", G, " ===")
  fit <- make_fit(G)
  p <- length(fit$coefficients)
  core <- refund:::pffr_influence(fit, "cl2", cl2_adjustment = "exact")
  ranks <- vapply(core$blocks, function(b) nrow(b$T), integer(1))
  message(
    "p = ",
    p,
    "  G = ",
    core$G,
    "  rank range = ",
    min(ranks),
    "-",
    max(ranks),
    "  sum(r_g) = ",
    sum(ranks)
  )

  # (a) core construction, cold cache
  ta <- med_time(
    {
      fit <- clear_influence_cache(fit)
      refund:::pffr_influence(fit, "cl2", cl2_adjustment = "exact")
    },
    reps
  )
  add_row(
    label = label,
    G = G,
    p = p,
    stage = "core",
    contrasts = NA_integer_,
    rows_timed = NA_integer_,
    seconds = ta$median,
    seconds_min = ta$min,
    seconds_max = ta$max,
    reps = ta$reps,
    seconds_full = ta$median,
    sum_rank = sum(ranks),
    max_rank = max(ranks)
  )
  message("  (a) core: ", signif(ta$median, 4), " s")

  # (b) coefficient grids through coef.pffr(); warm cache first, so the timing
  # is the repeated-use cost (df + interval assembly), not core construction.
  coef_call <- function() {
    coef(
      fit,
      sandwich = "cl2",
      cl2_adjustment = "exact",
      crit = "satterthwaite",
      ci = "pointwise"
    )
  }
  invisible(coef_call())
  nb <- count_df_contrasts(invisible(coef_call()))
  tb <- med_time(invisible(coef_call()), reps)
  add_row(
    label = label,
    G = G,
    p = p,
    stage = "coef_grid",
    contrasts = nb$contrasts,
    rows_timed = nb$contrasts,
    seconds = tb$median,
    seconds_min = tb$min,
    seconds_max = tb$max,
    reps = tb$reps,
    seconds_full = tb$median,
    sum_rank = sum(ranks),
    max_rank = max(ranks)
  )
  message(
    "  (b) coef grids: ",
    signif(tb$median, 4),
    " s for ",
    nb$contrasts,
    " contrasts in ",
    nb$calls,
    " calls"
  )

  # (c) E(Y)-type contrasts: full working-design rows.
  ctx <- refund:::pffr_df_context(fit, "cl2", cl2_adjustment = "exact")
  lp <- suppressWarnings(predict(fit, type = "lpmatrix"))
  n_all <- nrow(lp)
  take <- if (is.na(lp_rows) || lp_rows >= n_all) {
    seq_len(n_all)
  } else {
    unique(round(seq(1, n_all, length.out = lp_rows)))
  }
  lp_sub <- lp[take, , drop = FALSE]
  tc <- med_time(refund:::pffr_df_from_context(ctx, lp_sub), reps)
  add_row(
    label = label,
    G = G,
    p = p,
    stage = "ey_lpmatrix",
    contrasts = n_all,
    rows_timed = nrow(lp_sub),
    seconds = tc$median,
    seconds_min = tc$min,
    seconds_max = tc$max,
    reps = tc$reps,
    seconds_full = tc$median * n_all / nrow(lp_sub),
    sum_rank = sum(ranks),
    max_rank = max(ranks)
  )
  message(
    "  (c) E(Y) df: ",
    signif(tc$median, 4),
    " s for ",
    nrow(lp_sub),
    " of ",
    n_all,
    " rows (full-grid estimate ",
    signif(tc$median * n_all / nrow(lp_sub), 4),
    " s)"
  )

  # (d) Same-process A/B of the two df paths on identical geometry: the cached
  # residualization factor against the general product (the same core with the
  # factor switched off). This machine's load drifts between runs, so the
  # (a)-(c) comparison across two installed versions is noisy; measuring both
  # paths round robin in one process is not.
  if (isTRUE(core$df_precompute)) {
    ab_rows <- lp[
      round(seq(1, n_all, length.out = min(n_all, 2000L))),
      ,
      drop = FALSE
    ]
    gen <- ctx
    gen$core$df_precompute <- FALSE
    el <- matrix(NA_real_, reps, 2L, dimnames = list(NULL, c("gen", "fac")))
    for (i in seq_len(reps)) {
      el[i, "gen"] <- system.time(
        refund:::pffr_df_from_context(gen, ab_rows)
      )[["elapsed"]]
      el[i, "fac"] <- system.time(
        refund:::pffr_df_from_context(ctx, ab_rows)
      )[["elapsed"]]
    }
    mg <- stats::median(el[, "gen"])
    mf <- stats::median(el[, "fac"])
    for (nm in c("gen", "fac")) {
      add_row(
        label = label,
        G = G,
        p = p,
        stage = if (nm == "gen") "ab_general_path" else "ab_factored_path",
        contrasts = nrow(ab_rows),
        rows_timed = nrow(ab_rows),
        seconds = stats::median(el[, nm]),
        seconds_min = min(el[, nm]),
        seconds_max = max(el[, nm]),
        reps = reps,
        seconds_full = stats::median(el[, nm]) * n_all / nrow(ab_rows),
        sum_rank = sum(ranks),
        max_rank = max(ranks)
      )
    }
    message(
      "  (d) interleaved A/B on ",
      nrow(ab_rows),
      " rows: general ",
      signif(mg, 4),
      " s vs factored ",
      signif(mf, 4),
      " s -> ",
      signif(mg / mf, 3),
      "x"
    )
  }
  rm(fit, core, ctx, lp, lp_sub)
  gc()
}

out <- do.call(rbind, rows)
utils::write.csv(out, outfile, row.names = FALSE)
print(out)
