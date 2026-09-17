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
  # only 6 curves: effective rank <= 6 no matter how rich the generator is
  X <- low_rank_curves(6, s, rank = 20, sd_decay = 1)

  w <- tryCatch(
    ff(X, xind = s, splinepars = list(bs = "ps", k = c(8, 5))),
    warning = conditionMessage
  )
  expect_match(w, "cannot exceed min\\(nrow\\(X\\), ncol\\(X\\)\\) = 6")
})
