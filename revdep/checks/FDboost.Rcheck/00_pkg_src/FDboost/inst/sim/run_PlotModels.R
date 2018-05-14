###############################################################################
# run the simulations for boosting of functional data
# code based on code by Fabian Scheipl 
# author: Sarah Brockhaus
###############################################################################

print(R.Version()$version.string)

library(refund)
library(FDboost)
library(splines)
pathResults <- NULL
pathModels <- NULL

library(pryr) # to test for memory consumption

library(plyr)

# only works on Linux -> with try() no error on windows
try(library(doMC))
try(registerDoMC(cores=cores))

source("simfuns.R")


# ###################################### M=100

set.seed(18102012)

settings <- makeSettings(
  M=c(100),
  ni=c(1),
  Gy=c(30),
  Gx=c(30),
  snrEps=c(1,2),
  snrE=c(0),
  snrB=c(2),
  scenario=3,
  balanced=c(TRUE),
  nuisance=c(0),
  rep=1)

length(settings)

usecores <- 1
options(cores=usecores)
M100N1G30 <- try(doSim(settings=settings, cores=usecores))

save(M100N1G30, file=paste(pathResults, "plotModelsM100N1G30.Rdata", sep=""))


######################################

print(sessionInfo())
