###############################################################################
# run the simulations for FAMM of models with functional historical effects
# code based on code by Fabian Scheipl 
# author: Sarah Brockhaus
###############################################################################

rm(list=ls())

print(R.Version()$version.string)

## library(refundDevel)
library(refund)
library(FDboost)
library(splines)
library(MASS)
library(Matrix)

# setwd("../")
# path <- getwd()
# 
# source(paste0(getwd(),"/Code/boosting4.R"))    # functions for simulations
# 
# pathResults <- paste0(path, "/results/")
# pathModels <- paste0(path, "/models/")

source("boosting4.R")
pathResults <- NULL
pathModels <- NULL


getwd()

###################################### M=30

set.seed(18102012)

arginS <- "smooth"
penaltyS <- "ps"

set1 <- makeSettings(
  dgpsettings=list(M=c(30), 
                   ni=c(1),
                   p=c(0.8),
                   Gy=c(27),
                   Gx=c(100),
                   snrEps=c(2),
                   snrE=c(0),
                   snrB=c(2),
                   scenario=c(2),
                   k=c(0, 1, 2),
                   type=c("lines"), 
                   balanced=c(TRUE),
                   nuisance=c(18),  # 10
                   a=c("pen1coef4", "pen2coef4"),
                   regularS=c(TRUE), 
                   regularT=c(FALSE),
                   rep=1:20),  ########### rep=1:20)
  algorithms=list(addNoise=c(FALSE),
                  centerX=c(TRUE),
                  penaltyS=c("ps","pss"),
                  diffPen=c(1,2),
                  inS=c("smooth"))
)

length(set1)

usecores <- 5
options(cores=usecores)

# FAMM on set1
pffr1 <- try(doSimPffr(settings=set1, cores=usecores))
#pffr1
save(pffr1, file=paste(pathResults, "pffr1nuisance.Rdata", sep=""))



set.seed(18102012)

arginS <- "smooth"
penaltyS <- "ps"

set2 <- makeSettings(
  dgpsettings=list(M=c(30), 
                   ni=c(1),
                   p=c(0.8),
                   Gy=c(27),
                   Gx=c(100),
                   snrEps=c(2),
                   snrE=c(0),
                   snrB=c(2),
                   scenario=c(2),
                   k=c(5, 10),
                   type=c("bsplines"),
                   balanced=c(TRUE),
                   nuisance=c(8),  # 10
                   a=c("pen1coef4", "pen2coef4"),
                   regularS=c(TRUE), 
                   regularT=c(FALSE),
                   rep=1:20),
  algorithms=list(addNoise=c(FALSE),
                  centerX=c(TRUE),
                  penaltyS=c("ps","pss"),
                  diffPen=c(1,2),
                  inS=c("smooth"))
)

length(set2)

usecores <- 10
options(cores=usecores)

# FAMM on set2
pffr2 <- try(doSimPffr(settings=set2, cores=usecores))
#pffr2
save(pffr2, file=paste(pathResults, "pffr2nuisance.Rdata", sep=""))



set.seed(18102012)

arginS <- "smooth"
penaltyS <- "ps"

set3 <-  makeSettings(
  dgpsettings=list(M=c(30), 
                   ni=c(1),
                   p=c(0.8),
                   Gy=c(27),
                   Gx=c(100),
                   snrEps=c(2),
                   snrE=c(0),
                   snrB=c(2),
                   scenario=c(2),
                   k=c(5, 10),
                   type=c("local", "end"), ##, "start"
                   balanced=c(TRUE),
                   nuisance=c(8),  # 10
                   a=c("pen1coef4", "pen2coef4"),
                   regularS=c(TRUE), 
                   regularT=c(FALSE),
                   rep=1:20),
  algorithms=list(addNoise=c(FALSE),
                  centerX=c(TRUE),
                  penaltyS=c("ps","pss"),
                  diffPen=c(1,2),
                  inS=c("smooth"))
)
  
length(set3)

usecores <- 10
options(cores=usecores)

# FAMM on set3
pffr3 <- try(doSimPffr(settings=set3, cores=usecores))
#pffr3
save(pffr3, file=paste(pathResults, "pffr3nuisance.Rdata", sep=""))


set.seed(18102012)

arginS <- "smooth"
penaltyS <- "ps"


# set4 <-  makeSettings(
#   dgpsettings=list(M=c(30), 
#                    ni=c(1),
#                    p=c(0.8),
#                    Gy=c(27),
#                    Gx=c(100),
#                    snrEps=c(2),
#                    snrE=c(0),
#                    snrB=c(2),
#                    scenario=c(2),
#                    k=c(4),
#                    type=c("fourier"), 
#                    balanced=c(TRUE),
#                    nuisance=c(0),  # 10
#                    a=c("pen1coef4", "pen.1coef4"),
#                    regularS=c(TRUE), 
#                    regularT=c(FALSE),
#                    rep=1:20),
#   algorithms=list(addNoise=c(FALSE),
#                   centerX=c(TRUE, FALSE),
#                   penaltyS=c("ps","pss"),
#                   diffPen=c(1,2),
#                   inS=c("smooth"))
# )
# 
# length(set4)
# 
# usecores <- 10
# options(cores=usecores)
# 
# # FAMM on set4
# pffr4 <- try(doSimPffr(settings=set4, cores=usecores))
# #pffr4
# save(pffr4, file=paste(pathResults, "pffr4.Rdata", sep=""))
# 


######################################

print(sessionInfo())


