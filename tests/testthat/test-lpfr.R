# context("Testing old lpfr")
# 
# test_that("lpfr works with one predictor", {
#   skip_on_cran()
#   
#   data(DTI)
#   
#   # subset data as needed for this example
#   cca = DTI$cca[which(DTI$case == 1),]
#   rcst = DTI$rcst[which(DTI$case == 1),]
#   DTI = DTI[which(DTI$case == 1),]
#   
#   
#   # note there is missingness in the functional predictors
#   # apply(is.na(cca), 2, mean)
#   # apply(is.na(rcst), 2, mean)
#   
#   
#   # fit two models with single functional predictors and plot the results
#   fit.cca = lpfr(Y=DTI$pasat, subj=DTI$ID, funcs = cca, smooth.cov=FALSE)
#   fit.rcst = lpfr(Y=DTI$pasat, subj=DTI$ID, funcs = rcst, smooth.cov=FALSE)
#   ## expect_equal_to_reference(fit.cca$BetaHat, "lpfr.cca.coef.rds")
#   expect_is(fit.rcst, "list")
#   expect_equal(length(fit.rcst), 10)
# })
# 
# test_that("lpfr works two predictors", {
#   skip_on_cran()
#   
#   data(DTI)
#   
#   # subset data as needed for this example
#   cca = DTI$cca[which(DTI$case == 1),]
#   rcst = DTI$rcst[which(DTI$case == 1),]
#   DTI = DTI[which(DTI$case == 1),]
#   
#   # fit a model with two functional predictors and plot the results
#   fit.cca.rcst = lpfr(Y=DTI$pasat, subj=DTI$ID, funcs = list(cca,rcst),
#                       smooth.cov=FALSE)
#   expect_is(fit.cca.rcst, "list")
#   ## expect_equal_to_reference(fit.cca.rcst, "lpfr.fit.rds")
# })