#--------------------------------------
# ff() effective-rank / weak-identifiability guard
#--------------------------------------
#
# ff(check.ident = TRUE) warns when the effective rank of Cov(X(s)) --- the
# number of eigenvalues covering >= 99.5% of the variance --- falls below
# 1.5 * k_s, the marginal basis dimension along s. The historical check only
# fired at the much weaker rank < k_s, so weakly identified designs (part of
# beta(t, s) outside the span of the observed curves) were silent.

# n curves spanned by `rank` smooth components: effective rank ~ `rank`.
low_rank_curves <- function(n, s, rank, sd_decay = 0.9) {
  scores <- matrix(
    rnorm(n * rank, sd = rep(sd_decay^(seq_len(rank) - 1), each = n)),
    n,
    rank
  )
  basis <- sapply(seq_len(rank), function(j) sin(j * pi * s))
  scores %*% t(basis)
}

test_that("ff() warns when the effective rank is below 1.5 * k_s", {
  set.seed(20260804)
  s <- seq(0, 1, length.out = 40)
  X <- low_rank_curves(40, s, rank = 5)

  expect_warning(
    ff(X, xind = s, splinepars = list(bs = "ps", k = c(8, 5))),
    "Effective rank of <X>"
  )
  # message carries both numbers and the threshold
  w <- tryCatch(
    ff(X, xind = s, splinepars = list(bs = "ps", k = c(8, 5))),
    warning = conditionMessage
  )
  expect_match(w, "below 1.5 \\* k = 12")
  expect_match(w, "k = 8 basis functions")
  expect_match(w, "weakly identified")
})

test_that("the hard rank < k_s case is flagged inside the same warning", {
  set.seed(20260805)
  s <- seq(0, 1, length.out = 40)
  # rank ~ 6 (above the separate "very low rank" warning at <= 4) vs k_s = 10
  X <- low_rank_curves(40, s, rank = 6, sd_decay = 1)

  w <- tryCatch(
    ff(X, xind = s, splinepars = list(bs = "ps", k = c(10, 5))),
    warning = conditionMessage
  )
  expect_match(w, "identifiable only through the penalty")
})

test_that("ff() does not warn when the effective rank clears 1.5 * k_s", {
  set.seed(20260806)
  s <- seq(0, 1, length.out = 40)
  # rank ~ 15 against k_s = 5 (needs 7.5)
  X <- low_rank_curves(60, s, rank = 15, sd_decay = 1)

  expect_no_warning(
    ff(X, xind = s, splinepars = list(bs = "ps", k = c(5, 5)))
  )
})

test_that("check.ident = FALSE switches the guard off", {
  set.seed(20260804)
  s <- seq(0, 1, length.out = 40)
  X <- low_rank_curves(40, s, rank = 5)

  expect_no_warning(
    ff(
      X,
      xind = s,
      splinepars = list(bs = "ps", k = c(8, 5)),
      check.ident = FALSE
    )
  )
})

test_that("the guard survives a rank cap set by the number of curves", {
  set.seed(20260807)
  s <- seq(0, 1, length.out = 40)
  # only 6 curves: centred effective rank <= 5 however rich the generator is
  X <- low_rank_curves(6, s, rank = 20, sd_decay = 1)

  w <- tryCatch(
    ff(X, xind = s, splinepars = list(bs = "ps", k = c(8, 5))),
    warning = conditionMessage
  )
  expect_match(w, "cannot exceed min\\(nrow\\(X\\) - 1, ncol\\(X\\)\\) = 5")
})

# rank-5 centred variation around a common mean curve of size `mean_scale`:
# the mean adds one direction to the uncentred SVD (small scale) or swamps the
# variation (large scale), but never changes the rank of Cov(X(s)).
curves_with_mean <- function(
  mean_scale,
  n = 60,
  s = seq(0, 1, length.out = 40)
) {
  basis <- sapply(1:5, function(j) sin(j * pi * s))
  scores <- scale(matrix(rnorm(n * 5), n, 5), scale = FALSE)
  mean_curve <- cos(0.5 * pi * s) + s^2
  scores %*%
    t(basis) +
    mean_scale * matrix(mean_curve, n, length(s), byrow = TRUE)
}

effective_rank <- function(X) {
  ev <- svd(X, nu = 0, nv = 0)$d^2
  min(which(cumsum(ev) / sum(ev) >= 0.995))
}

test_that("the rank guard uses centred X: a mean curve does not add rank", {
  set.seed(1)
  s <- seq(0, 1, length.out = 40)
  X <- curves_with_mean(1, s = s)
  # fixture sanity: uncentred rank 6 would clear 1.5 * 4, centred rank 5 not
  expect_equal(effective_rank(X), 6)
  expect_equal(effective_rank(scale(X, scale = FALSE)), 5)

  w <- tryCatch(
    ff(X, xind = s, splinepars = list(bs = "ps", k = c(4, 5))),
    warning = conditionMessage
  )
  expect_match(w, "Effective rank of <X> is 5, below 1.5 \\* k = 6")
})

test_that("a dominant mean curve does not deflate the centred rank", {
  set.seed(1)
  s <- seq(0, 1, length.out = 40)
  X <- curves_with_mean(10, s = s)
  # uncentred rank 4 would trip the "very low effective rank" warning
  expect_equal(effective_rank(X), 4)
  expect_equal(effective_rank(scale(X, scale = FALSE)), 5)

  warns <- character()
  withCallingHandlers(
    ff(X, xind = s, splinepars = list(bs = "ps", k = c(4, 5))),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("Very low effective rank", warns)))
  expect_match(warns, "Effective rank of <X> is 5", all = FALSE)
})
