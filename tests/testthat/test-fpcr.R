context("Testing functional PCR")
library(refund)

test_that("Check that all 3 fpcr calls yield essentially identical estimates", {
   skip_on_cran()

   data(gasoline)

   # Create the requisite functional data objects
   bbasis = create.bspline.basis(c(900, 1700), 40)
   wavelengths = 2*450:850
   nir <- t(gasoline$NIR)
   gas.fd = smooth.basisPar(wavelengths, nir, bbasis)$fd

   # Method 1: Call fpcr with fdobj argument
   gasmod1 = fpcr(gasoline$octane, fdobj = gas.fd, ncomp = 30)
   plot(gasmod1, xlab="Wavelength")

   # Method 2: Call fpcr with explicit signal matrix
   gasmod2 = fpcr(gasoline$octane, xfuncs = gasoline$NIR, ncomp = 30)
   # Method 3: Call fpcr with explicit signal, basis, and penalty matrices
   gasmod3 = fpcr(gasoline$octane, xfuncs = gasoline$NIR,
                  basismat = eval.basis(wavelengths, bbasis),
                  penmat = getbasispenalty(bbasis), ncomp = 30)

   # Check that all 3 calls yield essentially identical estimates
   #expect_equal(gasmod1$fhat, gasmod2$fhat, gasmod3$fhat)
   tmp <- abs(c(gasmod1$fhat-gasmod2$fhat, gasmod1$fhat-gasmod2$fhat))
   expect_less_than(max(tmp), 1e-10)
   # But note that, in general, you'd have to specify argvals in Method 1
   # to get the same coefficient function values as with Methods 2 & 3.
})

test_that("Cross-validation is working", {
   skip_on_cran()
   set.seed(8659)
   data(gasoline)
   cv.gas = fpcr(gasoline$octane, xfuncs = gasoline$NIR,
                  nbasis=seq(20,40,5), ncomp = seq(10,20,5), store.cv = TRUE)
   expect_is(cv.gas, "fpcr")

   ## expect_equal_to_reference(cv.gas$cv.table, "fpcr.cv.rds")
})
