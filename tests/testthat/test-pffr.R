context("Testing pffr")
library(refundDevel)

set.seed(9312)
data2 <- pffrSim(scenario="all", n=200)
argvals <- attr(data2, "yindex")
s <- attr(data2, "xindex")
m2 <- pffr(Y ~  ff(X1, xind=s) + #linear function-on-function
##                ff(X2, xind=s) + #linear function-on-function
               xlin  +  #varying coefficient term
               c(te(xte1, xte2)) + #bivariate smooth term in xte1 & xte2, const. over Y-index
               s(xsmoo) + #smooth effect of xsmoo varying over Y-index
               c(xconst), # linear effect of xconst constant over Y-index
       yind=argvals,
          data=data2)

test_that("all major pffr terms are working", {
   ## expect_equal_to_reference(m.pc$coefficients, "pffr.all.coef.rds")
   expect_is(m2, "pffr")
})

test_that("convenience functions are working", {
   expect_is(summary(m2), "summary.pffr")
   ## plot(m2, pers=TRUE)
   expect_is(coef(m2), "list")
   # convenience functions:
   preddata <- pffrSim(scenario="all", n=20)
   expect_is(predict(m2, newdata=preddata), "matrix")

   expect_equal(length(predict(m2, type="terms")), 6L)
   cm2 <- coef(m2)
   expect_is(cm2$pterms, "matrix")
   # str(cm2$smterms, 2)
   # str(cm2$smterms[["s(xsmoo)"]]$coef)
})
 #############################################################################
                                        # sparse data (80% missing on a regular grid):
test_that("pffr with sparse data works", {
   set.seed(88182004)
   data3 <- pffrSim(scenario=c("int", "smoo"), n=100, propmissing=0.8)
   t <- attr(data3, "yindex")
   m3.sparse <- pffr(Y ~ s(xsmoo), data=data3$data, ydata=data3$ydata, yind=t)
   expect_is(summary(m3.sparse), "summary.pffr")
   ## plot(m3.sparse, pers=TRUE, pages=1)
})

set.seed(1122)
n <- 55
S <- 60
n.argvals <- 50
s <- seq(0,1, l=S)
argvals <- seq(0,1, l=n.argvals)

#generate X from a polynomial FPC-basis:
rankX <- 5
Phi <- cbind(1/sqrt(S), poly(s, degree=rankX-1))
lambda <- rankX:1
Xi <- sapply(lambda, function(l)
           scale(rnorm(n, sd=sqrt(l)), scale=FALSE))
X <- Xi %*% t(Phi)

beta.st <- outer(s, argvals, function(s, argvals) cos(2 * pi * s * argvals))

y <- (1/S*X) %*% beta.st + 0.1 * matrix(rnorm(n * n.argvals), nrow=n, ncol=n.argvals)

data <- list(y=y, X=X)

test_that("ffpc terms are working", {
   skip_on_cran()


   # set number of FPCs to true rank of process for this example:
   m.pc <- pffr(y ~ c(1) + 0 + ffpc(X, yind=argvals, decomppars=list(npc=rankX)),
                data=data, yind=argvals)
   ## expect_equal_to_reference(m.pc$coefficients, "pffr.ffpc.coef.rds")
   expect_is(m.pc, "pffr")
   expect_is(summary(m.pc), "summary.pffr")

   expect_is(ffpcplot(m.pc, type="surf", auto.layout=FALSE, theta = 50, phi = 40), "list")
})

test_that("another ff term example", {
   m.ff <- pffr(y ~ c(1) + 0 + ff(X, yind=argvals), data=data, yind=argvals)
   expect_equal_to_reference(m.ff$coefficients, "pffr.ff.coef.rds")

   expect_is(summary(m.ff), "summary.pffr")

   # fits are very similar:
   # expect_equal(fitted(m.pc), fitted(m.ff))  # NOT TRUE
})

test_that("ff limits arg works", {
  set.seed(2121)
  data <- pffrSim(scenario="ff", n=20, SNR=100, limits=function(s, t) s < t)
  t <- attr(data, "yindex")
  s <- attr(data, "xindex")
  m.l <- pffr(Y ~ ff(X1, xind=s, limits="s<t"), yind=t, data=data)
  m.leq <- pffr(Y ~ ff(X1, xind=s, limits="s<=t"), yind=t, data=data)
  m.fun <- pffr(Y ~ ff(X1, xind=s, limits=function(s, t) s < t), yind=t, data=data)
  expect_equal(fitted(m.l), fitted(m.fun))
  expect_that(
    max(abs((fitted(m.l) - fitted(m.leq))/fitted(m.l))) < .05)
})

test_that("weights and offset args work", {
  skip_on_cran()
  set.seed(112)
  n <- 20
  nygrid <- 50
  data <- pffrSim(scenario="lin", n=n, nygrid=nygrid, SNR=100)
  t <- attr(data, "yindex")
  m0 <- pffr(Y ~ xlin, yind=t, data=data)

  vecoffset <- 1:n
  data$Y_vecoffset <- data$Y + vecoffset
  m_vecoffset <- pffr(Y_vecoffset ~ xlin, yind=t, data=data, offset=vecoffset)
  expect_equal(fitted(m0), fitted(m_vecoffset)-vecoffset)

  matoffset <- matrix(rnorm(n*nygrid), n, nygrid)
  data$Y_matoffset <- data$Y + matoffset
  m_matoffset <- pffr(Y_matoffset ~ xlin, yind=t, data=data, offset=matoffset)
  expect_equal(fitted(m0), fitted(m_matoffset) - matoffset)

  expect_error(pffr(Y ~ xlin, yind=t, data=data, offset=vecoffset[-1]))
  expect_error(pffr(Y ~ xlin, yind=t, data=data, offset=matoffset[-1, ]))
})

test_that("sff terms are working", {
  skip_on_cran()

  set.seed(2121)
  data <- pffrSim(scenario="ff", n=20, SNR=100)
  t <- attr(data, "yindex")
  s <- attr(data, "xindex")
  m.s <- pffr(Y ~ sff(X1, xind=s), yind=t, data=data)
  m.lin <- pffr(Y ~ ff(X1, xind=s), yind=t, data=data)
  expect_is(summary(m.s), "summary.pffr")
  expect_that(
    max(abs((fitted(m.lin) - fitted(m.s))/fitted(m.lin))) < .05 )
}


  test_that("ff identifiability stuff works", {
    skip_on_cran()

    set.seed(112)
    n <- 20
    ngrid <- 50
    s <- t <- seq(0, 1, l=ngrid)
    Xefcts <- t(poly(t, 2))
    X <- matrix(rnorm(n*2), n, 2) %*% Xefcts
    L <- matrix(1/ngrid, ncol=ngrid, nrow=n)
    LX <- L*X
    beta.st <- outer(s, t, function(s, t) cos(3*pi*t)*s)
    Y <- LX%*%beta.st
    data <- data.frame(Y=I(Y), X=I(X))

    expect_warning(m <- pffr(Y ~ ff(X, xind=s), yind=t, data=data))
    expect_equal(class(m$smooth[[2]]$margin[[1]]), "pss.smooth")

})
