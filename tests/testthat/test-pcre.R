context("Testing pcre")
library(refundDevel)

test_that("pcre works as expected", {
  skip_on_cran()

  residualfunction <- function(t){
   #generate quintic polynomial error functions
       drop(poly(t, 5)%*%rnorm(5, sd=sqrt(2:6)))
   }
   # generate data Y(t) = mu(t) + E(t) + white noise
   set.seed(1122)
   n <- 50
   T <- 30
   t <- seq(0,1, l=T)
   # E(t): smooth residual functions
   E <- t(replicate(n, residualfunction(t)))
   int <- matrix(scale(3*dnorm(t, m=.5, sd=.5) - dbeta(t, 5, 2)), byrow=T, n, T)
   Y <- int + E + matrix(.2*rnorm(n*T), n, T)
   data <- data.frame(Y=I(Y))
   # fit model under independence assumption:
   summary(m0 <- pffr(Y ~ 1, yind=t, data=data))
   # get first 5 eigenfunctions of residual covariance
   # (i.e. first 5 functional PCs of empirical residual process)
   Ehat <- resid(m0)
   fpcE <- fpca.sc(Ehat, npc=5)
   ##expect_equal_to_reference(fpcE, "pcre.fpca.obj.rds")
   expect_is(fpcE, "list")

   efunctions <- fpcE$efunctions
   evalues <- fpcE$evalues
   data$id <- factor(1:nrow(data))
   # refit model with fpc-based residuals
   m1 <- pffr(Y ~ 1 + pcre(id=id, efunctions=efunctions, evalues=evalues, yind=t), yind=t, data=data)
   ## expect_equal_to_reference(m1, "pcre.pffr.obj.rds")
   expect_is(m1, "pffr")

   t1 <- predict(m1, type="terms")
   ## expect_equal_to_reference(t1, "pcre.prediction.obj.rds")
   expect_is(t1, "list")

   expect_is(summary(m1), "summary.pffr")
})

test_that("pcre works for sparse", {
  skip_on_cran()

  residualfunction <- function(t){
    #generate quintic polynomial error functions
    drop(poly(t, 5)%*%rnorm(5, sd=sqrt(2:6)))
  }
  # generate data Y(t) = mu(t) + E(t) + white noise
  set.seed(1122)
  n <- 50
  T <- 30
  t <- seq(0,1, l=T)
  # E(t): smooth residual functions
  E <- t(replicate(n, residualfunction(t)))
  int <- matrix(scale(3*dnorm(t, m=.5, sd=.5) - dbeta(t, 5, 2)), byrow=T, n, T)
  Y <- int + E + matrix(.2*rnorm(n*T), n, T)
  data <- data.frame(Y=I(Y))
  # fit model under independence assumption:
  summary(m0 <- pffr(Y ~ 1, yind=t, data=data))
  # get first 5 eigenfunctions of residual covariance
  # (i.e. first 5 functional PCs of empirical residual process)
  Ehat <- resid(m0)
  fpcE <- fpca.sc(Ehat, npc=5)
  ##expect_equal_to_reference(fpcE, "pcre.fpca.obj.rds")
  expect_is(fpcE, "list")

  efunctions <- fpcE$efunctions
  evalues <- fpcE$evalues

  # make sparse data:
  propmissing <- .5
  nygrid <- T
  missing <- sample(c(rep(T, propmissing*n*nygrid),
    rep(F, n*nygrid-propmissing*n*nygrid)))

  ydata <- data.frame(.obs = rep(1:n, each=nygrid)[!missing],
    .index = rep(t, times=n)[!missing],
    .value = as.vector(t(data$Y))[!missing])

  data <- data.frame(id = factor(1:nrow(data)))

  # refit model with fpc-based residuals
  m1 <- pffr(Y ~ 1 + pcre(id=id, efunctions=efunctions, evalues=evalues, yind=t), yind=t,
    data=data, ydata = ydata)
  ## expect_equal_to_reference(m1, "pcre.pffr.obj.rds")
  expect_is(m1, "pffr")

  t1 <- predict(m1, type="terms")
  ## expect_equal_to_reference(t1, "pcre.prediction.obj.rds")
  expect_is(t1, "list")

  expect_is(summary(m1), "summary.pffr")
})
