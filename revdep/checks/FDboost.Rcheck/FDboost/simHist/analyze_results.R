###############################################################################
# Author: Sarah Brockhaus
# Code inspired by Fabian Scheipl: lfpr3_analysis.R, 
# online appendix of paper: 
# Scheipl, F., Staicu, A.-M. and Greven, S. (2015). 
# Functional Additive Mixed Models. 
# Journal of Computational and Graphical Statistics
###############################################################################

library(ggplot2)
library(plyr)
library(RColorBrewer)

print("boosting_analysis.R")

# setwd("../results")
print(getwd())

# results FDboost
load("res1.Rdata", verbose=TRUE)
load("res2.Rdata", verbose=TRUE)
load("res3.Rdata", verbose=TRUE)

res123 <- rbind(res1, res2, res3)
res123[ , c("ONEx..t", "ONEx..tlong", paste0("X", c(1, 10, 2:9)))] <- NA
res123$stabsel <- "no"

# results FDboost with nuisance and stability selection
load("res1nuisance.Rdata", verbose=TRUE)
load("res2nuisance.Rdata", verbose=TRUE)
load("res3nuisance.Rdata", verbose=TRUE)

res123n <- rbind(res1n, res2n, res3n)
res123n$stabsel <- "yes"


# results FDboost with nuisance but without stability selection
load("res1nu.Rdata", verbose=TRUE)
load("res2nu.Rdata", verbose=TRUE)
load("res3nu.Rdata", verbose=TRUE)

res123nu <- rbind(res1nu, res2nu, res3nu)
res123nu$stabsel <- "no"


# results PFFR
load("pffr1.Rdata", verbose=TRUE)
pffr1 <- pffr1[,-1]
load("pffr2.Rdata", verbose=TRUE)
pffr2 <- pffr2[,-1]
load("pffr3.Rdata", verbose=TRUE)
pffr3 <- pffr3[,-1]

pffr123 <-  rbind(pffr1, pffr2, pffr3)
pffr123[ , c("ONEx..t", "ONEx..tlong", paste0("X", c(1, 10, 2:9)))] <- NA
pffr123$stabsel <- NA


#########################################################
## check that columns are the same
all(names(res123)==names(pffr123))
### merge results without nuisance variables into one data matrix
res <- rbind(res123, pffr123)

res$mod <- res$model
res$mod[res$inS=="linear"] <- 3
res$mod <- factor(res$mod, levels=1:2, labels = c("FAMM", "FDboost") )

res$typek <- paste(res$type, res$k, sep="-")

levels_typek <- c("local-5","local-10", 
                  "bsplines-5","bsplines-10",
                  "end-5","end-10",
                  #"start-5","start-10",
                  #"fourier-5", "fourier-10",
                  #"fourierLin-5", "fourierLin-10",
                  "lines-0","lines-1","lines-2")

res$typek <- factor(res$typek, levels=levels_typek)
res$dbeta <- 1
res$dbeta[res$a=="pen2coef4"] <- 2
res$dbeta <- factor(res$dbeta)

res$nuis <- paste( res$nuisance, "nuisance")

res$penAll <- paste(res$penaltyS, res$diffPen, sep="-")
res$aPen <- 1
res$aPen[res$a == "pen1coef4s"] <- 2
res$aPen <- factor(res$aPen)

#### do a grouping of the relative errors
res$relmsefx1b_grouped <- cut(res$relmsefx1b, breaks = c(0, 0.1, 1, 10, Inf))
table(res$relmsefx1b_grouped, useNA="ifany")

res$relmsefx2b_grouped <- cut(res$relmsefx2b, breaks = c(0, 0.1, 1, 10, Inf))
table(res$relmsefx2b_grouped, useNA="ifany")



#########################################################
### merge results with and without nuisance variables
## check that columns are the same
all(names(res123)==names(res123n))
all(names(res123)==names(res123nu))

### do not use stability selection
resNu <- rbind(res123, res123nu, pffr123)   ## , res123n

resNu$mod <- resNu$model
resNu$mod[resNu$inS=="linear"] <- 3
resNu$mod <- factor(resNu$mod, levels=1:2, labels = c("FAMM", "FDboost") )

resNu$typek <- paste(resNu$type, resNu$k, sep="-")

levels_typek <- c("local-5","local-10", 
                  "bsplines-5","bsplines-10",
                  "end-5","end-10",
                  #"start-5","start-10",
                  #"fourier-5", "fourier-10",
                  #"fourierLin-5", "fourierLin-10",
                  "lines-0","lines-1","lines-2")

resNu$typek <- factor(resNu$typek, levels=levels_typek)

resNu$dbeta <- 1
resNu$dbeta[resNu$a=="pen2coef4"] <- 2
resNu$dbeta <- factor(resNu$dbeta)

resNu$nuis <- paste( resNu$nuisance, "nuisance")


############################################################################
### look at estimation errors in the different areas of the coefficient surface
### diffSurface <- true_bst1 - est_bst1
diffbeta <- res[ ,grepl("diffSurface", names(res))]
res <- res[ ,!grepl("diffSurface", names(res))]

save(res, file="boosting.Rdata")


############################################################################

myplot <- function(perc, ...){
  perc <- matrix(perc, ncol=40)
  persp(perc, theta=30, phi=30, ticktype="detailed", xlab="\n s", ylab="\n t", zlab="", 
        nticks=4, ...)
} 


### estimation is especially bad at the edges 
plotSubset <- function(subset, sub=""){
  
  mse25 <- apply(diffbeta[subset,]^2, 2, quantile, probs=0.25)
  mse50 <- apply(diffbeta[subset,]^2, 2, quantile, probs=0.5)
  mse75 <- apply(diffbeta[subset,]^2, 2, quantile, probs=0.75)
  mse95 <- apply(diffbeta[subset,]^2, 2, quantile, probs=0.95)
  
  par(mfrow=c(2,2), mar=c(1, 0, 1, 0), cex=1.5, cex.main=1)
  #myplot(mse, main="mean MSE")
  #myplot(mse5, main="p5 MSE")
  myplot(mse25, main="p25 MSE")
  myplot(mse50, main="p50 MSE")
  myplot(mse75, main="p75 MSE")
  myplot(mse95, main="p95 MSE")
  
  title(sub=sub, cex.sub=1, line=-2, outer=TRUE)
  
}


#plotSubset(subset=res$diffPen==1)

pdf("diffSurfaceM30.pdf", width=10, height=10)
plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="lines" & res$k==0, sub="penalty=1, lines-0")
plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="lines" & res$k==1, sub="penalty=1, lines-1")
plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="lines" & res$k==2, sub="penalty=1, lines-2")

plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="bsplines" & res$k==5, sub="penalty=1, bsplines-5")
plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="bsplines" & res$k==10, sub="penalty=1, bsplines-10")

plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="local" & res$k==5, sub="penalty=1, local-5")
plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="local" & res$k==10, sub="penalty=1, local-10")

#plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="start" & res$k==5, sub="penalty=1, start-5")
#plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="start" & res$k==10, sub="penalty=1, start-10")

plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="end" & res$k==5, sub="penalty=1, end-5")
plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="end" & res$k==10, sub="penalty=1, end-10")

#plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="fourier" & res$k==5, sub="penalty=1, fourier-5")
#plotSubset(subset=res$diffPen==1 & res$M==30 & res$type=="fourier" & res$k==10, sub="penalty=1, fourier-10")
dev.off()


############################################################################
#### look at condition number of D_s
res$logCondDs <- as.numeric(levels(res$logCondDs))[res$logCondDs]
summary(res$logCondDs)
summary(res$logCondDs[is.finite(res$logCondDs)])

## for which settings is the condition number finite?
res$logCondDsFinite <- res$logCondDs
res$logCondDsFinite[is.finite(res$logCondDs)] <- 0
res$logCondDsFinite[!is.finite(res$logCondDs)] <- 1
## table for % of finite condition numbers per setting
round(with(res, prop.table(table(typek, logCondDsFinite), margin=1))*100, 2)

res$logCondDs[!is.finite(res$logCondDs)] <- 25
res$logCondDs_grouped <- cut(res$logCondDs, breaks = c(0, log10(10^6), Inf), right=FALSE) 
summary(res$logCondDs_grouped)
plot(res$logCondDs, col=res$type)
abline(h=log10(10^6))

res$condNr <- res$logCondDs_grouped


#### look at condition number of submatrices
tused <- round((1:27-1)^2/(27-1)^2*15 + 1, 2)[1:26]
namescondNrt <- paste0("logCondDs_hist.", tused )

for(i in 1:length(namescondNrt)){
  res[,namescondNrt[i]]  <- as.numeric(levels(res[,namescondNrt[i]]))[res[,namescondNrt[i]]]
}

## look at course of the condition number over t
## only makes sense for setting that have not Inf as condition number anyway
# logCondDsFinite
# typek            0      1
# local-5     100.00   0.00
# local-10    100.00   0.00
# bsplines-5    8.75  91.25
# bsplines-10 100.00   0.00
# end-5         0.00 100.00
# end-10        0.00 100.00
# lines-0       0.00 100.00
# lines-1       0.00 100.00
# lines-2       3.75  96.25


pdf("logCondTime.pdf")
par(mfrow=c(1,1), mar=c(2, 3, 1, 1), cex=2)
for(i in c(1:4,09 )){
  temp <- res[res$typek==levels_typek[i], namescondNrt]
  matplot(tused, t(temp), type="l", xlab="t", main=levels_typek[i],  
          ylab=expression(kappa[j](t)), col=1)
  abline(h=6, col=2)
  rug(tused, 0.01)
}
dev.off()

## local-5, local-10 and bsplines-10 show similar behaviour
## local-5, local-10: k_j(t) < 6, for t > 2
## bsplines-10: k_j(t) < 6, for t > 5

## bsplines-5 and lines-2 have a lot of infinite values
## all condition number are > 6


############################################################################
### look at overlap measures
### use the maximum of the cumulative overlap with 10 parts, successive overlap
namesOverlap <- paste0("cumOverlapKe", 1:27)
for(i in 1:length(namesOverlap)){
  res[,namesOverlap[i]]  <- as.numeric(levels(res[,namesOverlap[i]]))[res[,namesOverlap[i]]]
}

## overlap does not make sense for shrinkage penalty
# summary(res[res$penaltyS=="ps", namesOverlap]) 

tused27 <- round((1:27-1)^2/(27-1)^2*15 + 1, 2)


pdf("overlapTime.pdf")
par(mfrow=c(1,1), mar=c(2, 3, 1, 1), cex=2)
# par(mfrow=c(3,3))
for(i in 1:length(levels_typek)){
  temp <- res[res$typek==levels_typek[i], namesOverlap]
  penTemp <- res$diffPen[res$typek==levels_typek[i]]
  ## set infinite values to 25
  temp[(as.matrix(temp))==5] <- 3
  matplot(tused27, t(temp)+runif(320*27, min=-0.05, max=0.05), type="l", xlab="t", main=levels_typek[i],  
          ylab="overlap", col=penTemp, ylim=c(-0.1, 1.5))
  abline(h=1, col=4)
  rug(tused, 0.01)
  legend("topleft", fill=1:2, title="penalty order", legend=1:2)
}
dev.off()

res$overlapMaxSuc10 <- suppressWarnings(apply(res[, namesOverlap], 1, max, na.rm=TRUE))
res$overlapMaxSuc10[(!is.finite(res$overlapMaxSuc10))] <- NA
summary(res$overlapMaxSuc10)
res$overMaxSuc10 <- cut(res$overlapMaxSuc10, breaks = c(0, 1, Inf), right=FALSE)
summary(res$overMaxSuc10)  
## maximal kernel overlap is computed in check_ident() as well 
all(res$overlapKe == res$overMaxSuc10)

## NA for local-5, local-10, bslines-10, lines-0, lines-1, lines-2 in combination with pss
## do not look at kernel overlap for pss-penalty, just for ps-penalty!
#table(res$overMaxSuc10, res$penAll, res$typek, useNA="ifany")
table(res$overMaxSuc10, res$penaltyS, useNA="ifany")

### use the kernel overlap with the design matrix in s-direction, D_s
res[,"overlapKeComplete"]  <- as.numeric(levels(res[,"overlapKeComplete"]))[res[,"overlapKeComplete"]]
summary(res$overlapKeComplete)
res$overlapKeComplete <- cut(res$overlapKeComplete, breaks = c(0, 1, Inf), right=FALSE)
summary(res$overlapKeComplete)
table(res$typek, res$overlapKeComplete, res$penaltyS)

### delete extra overlap-measures
#res <- res[ , !grepl("overlapKe", names(res)) | names(res) == "overlapKeXPbot2"]

############################################################################

## look at range of errors for response 
summary(res$relmsey[res$mod=="FAMM"])
summary(res$relmsey[res$mod=="FDboost"])

########## use relmsefx1b
# the correlation between goodness of fit for y and beta is very small
ddply(res, .(mod), summarize, cor.y.beta = cor(relmsey, relmsefx1b, method="spearman")) 
ddply(res, .(mod), summarize, cor.y.beta = cor(relmsey, relmsefx2b, method="spearman")) 

ddply(res, .(mod), summarize, cor.y.beta = cor(relmsey, relmsefx1b, method="spearman")) 
ddply(res, .(mod), summarize, cor.y.beta = cor(relmsey, relmsefx2b, method="spearman"))

# compute rank correlation between reliMSE beta and reliMSE Y for each combination
tab <- ddply(res, .(type, k, M, penaltyS, diffPen, model), summarize, 
             cor.y.beta = cor(relmsey, relmsefx1b, method="spearman")) 
tab
#sort(tab$cor.y.beta)


############################################################################
#### mosaicplot for flagged and grouped rel. error

#### define the variable flagged using overMaxSuc10
res$flagged <- NA
res$flagged[res$condNr=="[0,6)"] <- 0
res$flagged[res$condNr=="[6,Inf)"] <- 1
res$flagged[res$overMaxSuc10=="[1,Inf)" & res$condNr=="[6,Inf)"] <- 2
res$flagged <- factor(res$flagged, levels=0:2, labels=c("no", "cond" ,"yes"))
## only look at ps penalty!!!!
res$flagged[res$penaltyS=="pss"] <- NA

table(res$flagged, useNA="always")
with(res[res$penaltyS=="ps", ], table(res$condNr, res$overMaxSuc10, res$flagged, useNA="always"))


## & M==30 & mod=="FDboost" & a=="pen1coef4" & penaltyS=="ps"
### Figure 2
pdf("flagged.pdf", width=11, height=7)
par(mfrow=c(1,1), mar=c(2, 3, 1, 1), cex=2)
tabFlag1 <- with(res[!is.na(res$flagged) & res$M==30 & res$mod=="FDboost" & res$a=="pen1coef4" & res$penaltyS=="ps",], 
                 table(relmsefx1b_grouped, flagged, useNA="ifany"))
mosaicplot(tabFlag1, main="", col=c("white","grey80", "grey60"), xlab=bquote(reliMSE(beta[1])))
dev.off()

prop.table(tabFlag1, 2)


############################################################################
## boxplots for goodness of identifiability measures

### add fake data to res, so that all boxes have the same width and location 
### for overMaxSuc10
temp <- with(res, expand.grid(M=30, mod=unique(mod), a=unique(a),
                              overMaxSuc10=na.omit(unique(overMaxSuc10)), 
                              #overKeXPbot2=na.omit(unique(overKeXPbot2)),
                              diffPen=unique(diffPen), condNr=unique(condNr), 
                              typek=unique(typek), penaltyS="ps"))

## only keep combinations in temp that do not exist in res
tempPaste <- apply(temp, 1, paste, collapse="_", sep="")
resPaste <- apply(res[,c("M","mod","a","overMaxSuc10","diffPen","condNr","typek", "penaltyS")], 
                  1, paste, collapse="_", sep="")
sum(!tempPaste %in% resPaste)
length(tempPaste)
temp <- temp[!tempPaste %in% resPaste, ]

temp$penAll <- paste(temp$penaltyS, temp$diffPen, sep="-")

help <- rbind.fill(res, temp)
help$relmsefx1b[is.na(help$relmsefx1b)] <- 10^8

rangefx1 <- range(res$relmsefx1b, na.rm=TRUE)

with(res[res$penaltyS=="ps", ], table(overMaxSuc10, useNA="ifany"))

## check for missings in relmsefx1b
summary(res$relmsefx1b)
## as there are no missing values, it is ok to set all missing values to 0.1

penalty_names <- list(
  '1'="pen 1",
  '2'="pen 2"
)


help$overMaxSuc10Nice <- "1"
help$overMaxSuc10Nice[help$overMaxSuc10=="[1,Inf)"] <- "2"

help$overMaxSuc10Nice <- factor(help$overMaxSuc10Nice, levels=1:2, labels=c("<1", ">1"))

help$d <- factor(help$diffPen, levels=1:2, labels=1:2) # diffPen is penalty order of estimation!

### delte the impossible combination: true for global kernel overlap!!
help <- help[!(help$condNr == "[0,6)" & help$overMaxSuc10=="[1,Inf)"), ]

## local-5 and local-10 have maximal kernel overlap>1, but no problem in condition number!!

## look at condition number and kernel overlap
## do separate plots for condition number< 10^6 and condition number > 10^6
pdf("condNrOverlap0.pdf", width=5, height=6) 
par(mfrow=c(1,1), cex=1.5)
####### use overMaxSuc10
## do plot for condition number < 10^6
## typek %in% c("bsplines-10","local-5","local-10"
print(p1 <- ggplot(subset(help, condNr == "[0,6)" & !is.na(overMaxSuc10Nice) & M==30 & 
                      mod=="FDboost" & a=="pen1coef4" & penaltyS=="ps"), 
             aes(y=relmsefx1b, fill=typek, colour=typek, x=overMaxSuc10Nice)) + # & k>3
        geom_hline(aes(yintercept=0.1)) +
        geom_boxplot(outlier.size=.6) +
        facet_grid(~penAll) + 
        scale_y_log10() +  
        coord_cartesian(ylim = rangefx1) +                      
        scale_fill_manual(name = "", 
                          values=c(brewer.pal(8, "Paired")[c(1:4, 7:8)], brewer.pal(8, "PRGn")[1:3] )) +              
        scale_colour_manual(name = "", values=c(rep(c("black", "grey60"), 3), rep("grey30", 3)) ) +
        labs(title=bquote(kappa[1]<10^6) ) +
        ylab(bquote(reliMSE(beta[1]))) +  
        xlab("kernel overlap") +
        scale_x_discrete(labels=c("<1")) + 
        theme(text=element_text(size = 25, colour=1), plot.title=element_text(size = 25), 
              axis.text=element_text(size = 20, colour=1)) + 
        theme(legend.position = "none") # no legend 
)
dev.off()

pdf("condNrOverlap1.pdf", width=11, height=6) # 
print(p2 <- ggplot(subset(help, condNr == "[6,Inf)" & !is.na(overMaxSuc10Nice) & M==30 & 
                      mod == "FDboost" & a =="pen1coef4" & penaltyS=="ps"), 
             aes(y=relmsefx1b, fill=typek, colour=typek, x=overMaxSuc10Nice)) + # & k>3
        geom_hline(aes(yintercept=0.1)) + 
        geom_boxplot(outlier.size=.6) +
        facet_grid(~penAll) + 
        scale_y_log10() +  
        coord_cartesian(ylim = rangefx1) +                     
        scale_fill_manual(name = "", 
                          values=c(brewer.pal(8, "Paired")[c(1:4, 7:8)], brewer.pal(8, "PRGn")[1:3] )) +              
        scale_colour_manual(name = "", values=c(rep(c("black", "grey60"), 3), rep("grey30", 3)) ) +
        labs(title=expression(kappa[1]~"">=10^6)) +
        ylab("") + 
        xlab("kernel overlap") +
        scale_x_discrete(labels=c("<1", expression("">=1))) +
        theme(text=element_text(size = 25, colour=1), plot.title=element_text(size = 25), 
              axis.text=element_text(size = 20, colour=1), axis.text.y = element_blank()) +  
        theme(legend.title=element_blank() ) # no legend title
)
dev.off()

### Create Figure 1 of paper as pdf and as eps
if(FALSE){
  require(gridExtra)
  
  setEPS()
  postscript("Fig1.eps", width=16, height=6)
  grid.arrange(p1, p2, ncol=2, widths=c(5, 11))
  dev.off()
  
  pdf("Fig1.pdf", width=16, height=6)
  grid.arrange(p1, p2, ncol=2, widths=c(5, 11))
  dev.off()
}



############################################################################
### do at plot for the different penalties
with(res, table(typek, condNr, overMaxSuc10))

### Figure 3
# look at different penalties
pdf("penalties.pdf", width=10, height=5)
par(mfrow=c(1,1))
## order of penalty
print(ggplot(subset(res, !is.na(relmsefx1b) & M==30 & mod=="FDboost"), 
             aes(y=relmsefx1b, fill=typek, color=typek, x=dbeta)) + 
        geom_boxplot(outlier.size=.6) +
        facet_wrap(  ~ penAll, ncol=2) +
        scale_y_log10() + 
        scale_fill_manual(name = "", 
                          values=c(brewer.pal(8, "Paired")[c(1:4, 7:8)], brewer.pal(8, "PRGn")[1:3] )) +              
        scale_colour_manual(name = "", values=c(rep(c("black", "grey60"), 3), rep("grey30", 3)) ) +
        ylab(bquote(reliMSE(beta[1]))) + 
        xlab(bquote(d[beta])) + 
        theme(text=element_text(size = 25, colour=1), plot.title=element_text(size = 25), 
              axis.text=element_text(size = 20, colour=1)) +
        geom_hline(aes(yintercept=0.1)))

## order of penalty
print(ggplot(subset(res, !is.na(relmsefx1b) & M==30 & mod=="FDboost"), 
             aes(y=relmsefx1b, fill=penAll, color=penAll, x=dbeta)) + 
        geom_boxplot(outlier.size=.6) +
        facet_wrap(  ~ typek, ncol=3) +
        scale_y_log10() + 
        scale_fill_manual(name="penalty", values=c("blue", "lightblue", "darkgreen", "lightgreen")) +                    
        scale_colour_manual(name="penalty", values=c("black", "grey40", "black", "grey40")) + 
        ylab(bquote(reliMSE(beta[1]))) + 
        xlab(bquote(d[beta])) + 
        theme(text=element_text(size = 25, colour=1), plot.title=element_text(size = 25), 
              axis.text=element_text(size = 20, colour=1)) +
        geom_hline(aes(yintercept=0.1)))
dev.off()


############################################################################
### do a plot comparing FDboost and pffr without nuisance variables

### focusing on "best" settings
pdf("reliMSEfx1fbest.pdf", width=10, height=5)
print(ggplot(subset(res, !is.na(relmsefx1b) & !is.na(relmsefx1b & M==30) & penaltyS=="ps" & diffPen==1), 
             aes(y=relmsefx1b, fill=typek, color=typek, x=dbeta)) + # & k>3
        geom_boxplot(outlier.size=.6) +
        facet_wrap(  ~ mod, ncol=2) +  ## nuisance
        scale_y_log10() + 
        scale_fill_manual(name = "", 
                          values=c(brewer.pal(8, "Paired")[c(1:4, 7:8)], brewer.pal(8, "PRGn")[1:3] )) +              
        scale_colour_manual(name = "", values=c(rep(c("black", "grey60"), 3), rep("grey30", 3)) ) +
        ylab(bquote(reliMSE(beta[1]))) + 
        xlab(bquote(d[beta])) + 
        theme(text=element_text(size = 25, colour=1), plot.title=element_text(size = 25), 
              axis.text=element_text(size = 20, colour=1)) +
        geom_hline(aes(yintercept=0.1)))

dev.off()


############################################################################
### do a plot comparing FDboost and pffr WITH nuisance variables

### focusing on "best" settings
pdf("reliMSEfx1fbestNuisance.pdf", width=10, height=6)

print(ggplot(subset(resNu, !is.na(relmsefx1b) & !is.na(relmsefx1b & M==30) & penaltyS=="ps" & diffPen==1), 
             aes(y=relmsefx1b, fill=typek, color=typek, x=dbeta)) + # & k>3
        geom_boxplot(outlier.size=.6) + 
        facet_wrap( nuis ~ mod, ncol=2, nrow=2, drop=FALSE) +  ## nuisance
        scale_y_log10() + 
        scale_fill_manual(name = "", 
                          values=c(brewer.pal(8, "Paired")[c(1:4, 7:8)], brewer.pal(8, "PRGn")[1:3] )) +              
        scale_colour_manual(name = "", values=c(rep(c("black", "grey60"), 3), rep("grey30", 3)) ) +
        ylab(bquote(reliMSE(beta[1]))) + 
        xlab(bquote(d[beta])) + 
        theme(text=element_text(size = 25, colour=1), plot.title=element_text(size = 25), 
              axis.text=element_text(size = 20, colour=1)) +
        geom_hline(aes(yintercept=0.1)))

dev.off()


################ look at results with stability seleciton
table(res123n$X1>0, useNA="ifany")/nrow(res123n) # always selected
table(res123n$X2>0, useNA="ifany")/nrow(res123n) # always selected

table(res123n$X3>0, useNA="ifany")/nrow(res123n)
table(res123n$X4>0, useNA="ifany")/nrow(res123n)
table(res123n$X5>0, useNA="ifany")/nrow(res123n)
table(res123n$X6>0, useNA="ifany")/nrow(res123n)
table(res123n$X7>0, useNA="ifany")/nrow(res123n)
table(res123n$X8>0, useNA="ifany")/nrow(res123n)
table(res123n$X9>0, useNA="ifany")/nrow(res123n)
table(res123n$X10>0, useNA="ifany")/nrow(res123n)




################ look at computation time

# ## use the data containing boosting with and without stability selection 
temp123 <- res123[ , c("time.elapsed", "M", "type", "k", "a", "nuisance", "model", "stabsel")]
temp_pffr123 <- pffr123[ , c("time.elapsed", "M", "type", "k", "a", "nuisance", "model", "stabsel")]
temp_pffr123$stabsel <- "no"
temp123n <- res123n[ , c("time.elapsed", "M", "type", "k", "a", "nuisance", "model", "stabsel")]
temp123nu <- res123nu[ , c("time.elapsed", "M", "type", "k", "a", "nuisance", "model", "stabsel")]

resTime <- rbind(temp123, temp_pffr123, temp123n, temp123nu)

resTime$mod <- resTime$model
resTime$mod <- factor(resTime$mod, levels=1:2, labels = c("FAMM", "FDboost") )

resTime$typek <- paste(resTime$type, resTime$k, sep="-")

levels_typek <- c("local-5","local-10", 
                  "bsplines-5","bsplines-10",
                  "end-5","end-10",
                  #"start-5","start-10",
                  #"fourier-5", "fourier-10",
                  #"fourierLin-5", "fourierLin-10",
                  "lines-0","lines-1","lines-2")

resTime$typek <- factor(resTime$typek, levels=levels_typek)

resTime$dbeta <- 1
resTime$dbeta[resTime$a=="pen2coef4"] <- 2
resTime$dbeta <- factor(resTime$dbeta)

resTime$N <- factor(resTime$M)
resTime$nuis <- paste( resTime$nuisance, "nuisance")
resTime$stabsel_plot <- paste0("stability selection: ", resTime$stabsel)

pdf("computation_time.pdf", width=12, height=5)
print(pt <- ggplot(resTime, 
                   aes(y=time.elapsed, fill=mod, colour=mod, x=mod)) + # & k>3
        geom_boxplot(outlier.size=.6) +
        facet_grid(~ nuis + stabsel_plot) + 
        scale_y_continuous(breaks=c(1, 2, 5, 10, 20, 60, 120, 300, 600, 1200, 2700, 
                                    5400, 10800, 21600, 43200), 
                           trans="log10",
                           labels=c("1s", "2s", "5s", "10s", "20s", "1 min", "2 min", "5 min", 
                                    "10 min", "20 min", "45 min", "90 min", "3h", "6h", "12h")) + 
        scale_fill_manual(name = "",  values=c("white", "grey80")) +              
        scale_colour_manual(name = "", values=c("grey40", "black")) +
        ylab("computation time") + 
        xlab("estimation algorithm") +
        theme(text=element_text(size = 25, colour=1), plot.title=element_text(size = 25), 
              axis.text=element_text(size = 20, colour=1)) +  
        theme(legend.title=element_blank() ) # no legend title
)
dev.off() 


