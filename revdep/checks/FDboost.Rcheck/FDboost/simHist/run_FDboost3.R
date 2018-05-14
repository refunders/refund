###############################################################################
# run the simulations for boosting of models with functional historical effects
# code based on code by Fabian Scheipl 
# author: Sarah Brockhaus
###############################################################################

rm(list=ls())

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

doStabsel <- TRUE


###################################### M=30, Gy=27, but p=0.8!

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
                   nuisance=c(8),  # add 8 nuisance variables 
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

length(set1)

usecores <- 10
options(cores=usecores)

# boosting on set1
res1n <- try(doSimFDboost(settings=set1))
#res1n
save(res1n, file=paste(pathResults, "res1nuisance.Rdata", sep=""))



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

# boosting on set2
res2n <- try(doSimFDboost(settings=set2))
#res2n
save(res2n, file=paste(pathResults, "res2nuisance.Rdata", sep=""))


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

# boosting on set3
res3n <- try(doSimFDboost(settings=set3))
#res3n
save(res3n, file=paste(pathResults, "res3nuisance.Rdata", sep=""))



### do not fit the fourier-settings
# set.seed(18102012)
# 
# arginS <- "smooth"
# penaltyS <- "ps"
# 
# 
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
#                    nuisance=c(8),  # 10
#                    a=c("pen1coef4", "pen2coef4"),
#                    regularS=c(TRUE), 
#                    regularT=c(FALSE),
#                    rep=1:20),
#   algorithms=list(addNoise=c(FALSE),
#                   centerX=c(TRUE),
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
# # boosting on set4
# res4n <- try(doSimFDboost(settings=set4))
# #res4n
# save(res4n, file=paste(pathResults, "res4nuisance.Rdata", sep=""))
# 


# ######################################

print(sessionInfo())

