context("Test fosr")
library(refundDevel)

test_that("fosr Penalized GLS works", {
  skip_on_cran()
  require(fda)
  ## The first two lines, adapted from help(fRegress) in package fda,
  ## set up a functional data object representing daily average
  ## temperatures at 35 sites in Canada
  daybasis25 <- create.fourier.basis(rangeval=c(0, 365), nbasis=25,
                   axes=list('axesIntervals'))
  Temp.fd <- with(CanadianWeather, smooth.basisPar(day.5,
                 dailyAv[,,'Temperature.C'], daybasis25)$fd)

  modmat = cbind(1, model.matrix(~ factor(CanadianWeather$region) - 1))
  constraints = matrix(c(0,1,1,1,1), 1)

  ## Penalized GLS
  glsmod = fosr(fdobj = Temp.fd, X = modmat, con = constraints, method="GLS")
  expect_that(glsmod, is_a("fosr"))
  plot(glsmod, 1)
})
