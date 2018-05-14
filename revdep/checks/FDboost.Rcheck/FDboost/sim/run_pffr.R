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

# install.packages("pryr")
library(pryr)

source("simfuns.R")

# ##################################### M=50
# 
# set.seed(18102012)
# 
# settingsM50N1G30 <- makeSettings(
#   M=c(50),
#   ni=c(1),
#   Gy=c(30),
#   Gx=c(30),
#   snrEps=c(1,2),
#   snrE=c(0),
#   snrB=c(2),
#   scenario=3,
#   balanced=c(TRUE),
#   nuisance=c(0),
#   rep=1:10)
# 
# length(settingsM50N1G30)
# 
# usecores <- 2
# options(cores=usecores)
# M50N1G30 <- try(doSimPffr(settings=settingsM50N1G30, cores=usecores))
# 
# save(M50N1G30, file=paste(pathResults, "M50N1G30pffr.Rdata", sep=""))
# 
# 
# set.seed(18102012)
# 
# settingsM50N1G100 <- makeSettings(
#   M=c(50),
#   ni=c(1),
#   Gy=c(100),
#   Gx=c(100),
#   snrEps=c(1, 2),
#   snrE=c(0),
#   snrB=c(2),
#   scenario=3,
#   balanced=c(TRUE),
#   nuisance=c(0),
#   rep=1:10)
# 
# length(settingsM50N1G100)
# 
# usecores <- 2
# options(cores=usecores)
# M50N1G100 <- try(doSimPffr(settings=settingsM50N1G100, cores=usecores))
# 
# save(M50N1G100, file=paste(pathResults, "M50N1G100pffr.Rdata", sep=""))


##################################### M=100

set.seed(18102012)

settingsM100N1G30 <- makeSettings(
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
  rep=1:10)

length(settingsM100N1G30)

usecores <- 2
options(cores=usecores)
M100N1G30 <- try(doSimPffr(settings=settingsM100N1G30, cores=usecores))

save(M100N1G30, file=paste(pathResults, "M100N1G30pffr.Rdata", sep=""))
#save(M50N1G30, file=paste(pathResults, "test.Rdata", sep=""))


set.seed(18102012)

settingsM100N1G100 <- makeSettings(
  M=c(100),
  ni=c(1),
  Gy=c(100),
  Gx=c(100),
  snrEps=c(1,2),
  snrE=c(0),
  snrB=c(2),
  scenario=3,
  balanced=c(TRUE),
  nuisance=c(0),
  rep=1:10)

length(settingsM100N1G100)

usecores <- 2
options(cores=usecores)
M100N1G100 <- try(doSimPffr(settings=settingsM100N1G100, cores=usecores))

save(M100N1G100, file=paste(pathResults, "M100N1G100pffr.Rdata", sep=""))


# ##################################### M=200
# 
# set.seed(18102012)
# 
# settingsM200N1G30 <- makeSettings(
#   M=c(200),
#   ni=c(1),
#   Gy=c(30),
#   Gx=c(30),
#   snrEps=c(1,2),
#   snrE=c(0),
#   snrB=c(2),
#   scenario=3,
#   balanced=c(TRUE),
#   nuisance=c(0),
#   rep=1:10)
# 
# length(settingsM200N1G30)
# 
# usecores <- 2
# options(cores=usecores)
# M200N1G30 <- try(doSimPffr(settings=settingsM200N1G30, cores=usecores))
# 
# save(M200N1G30, file=paste(pathResults, "M200N1G30pffr.Rdata", sep=""))
# #save(M50N1G30, file=paste(pathResults, "test.Rdata", sep=""))
# 
# 
# set.seed(18102012)
# 
# settingsM200N1G100 <- makeSettings(
#   M=c(200),
#   ni=c(1),
#   Gy=c(100),
#   Gx=c(100),
#   snrEps=c(1,2),
#   snrE=c(0),
#   snrB=c(2),
#   scenario=3,
#   balanced=c(TRUE),
#   nuisance=c(0),
#   rep=1:10)
# 
# length(settingsM200N1G100)
# 
# usecores <- 2
# options(cores=usecores)
# M200N1G100 <- try(doSimPffr(settings=settingsM200N1G100, cores=usecores))
# 
# save(M200N1G100, file=paste(pathResults, "M200N1G100pffr.Rdata", sep=""))
# 


##################################### M=500

set.seed(18102012)

settingsM500N1G30 <- makeSettings(
  M=c(500),
  ni=c(1),
  Gy=c(30),
  Gx=c(30),
  snrEps=c(1,2),
  snrE=c(0),
  snrB=c(2),
  scenario=3,
  balanced=c(TRUE),
  nuisance=c(0),
  rep=1:10)

length(settingsM500N1G30)

usecores <- 2
options(cores=usecores)
M500N1G30 <- try(doSimPffr(settings=settingsM500N1G30, cores=usecores))

save(M500N1G30, file=paste(pathResults, "M500N1G30pffr.Rdata", sep=""))


set.seed(18102012)

settingsM500N1G100 <- makeSettings(
  M=c(500),
  ni=c(1),
  Gy=c(100),
  Gx=c(100),
  snrEps=c(1, 2),
  snrE=c(0),
  snrB=c(2),
  scenario=3,
  balanced=c(TRUE),
  nuisance=c(0),
  rep=1:10)

length(settingsM500N1G100)

usecores <- 2
options(cores=usecores)
M500N1G100 <- try(doSimPffr(settings=settingsM500N1G100, cores=usecores))

save(M500N1G100, file=paste(pathResults, "M500N1G100pffr.Rdata", sep=""))


######################################

print(sessionInfo())