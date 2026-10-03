expect_valid_truth <- function(dat, n, nygrid) {
  truth <- attr(dat, "truth")
  expect_true(is.list(truth), info = "truth should be a list")
  expect_true(
    all(c("eta", "etaTerms", "epsilon") %in% names(truth)),
    info = "truth should contain eta, etaTerms, epsilon"
  )
  expect_equal(
    dim(truth$eta),
    c(n, nygrid),
    info = "truth$eta dimensions should match n x nygrid"
  )
  expect_equal(
    dim(truth$epsilon),
    c(n, nygrid),
    info = "truth$epsilon dimensions should match n x nygrid"
  )
}
