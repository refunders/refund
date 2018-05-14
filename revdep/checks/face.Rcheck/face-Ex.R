pkgname <- "face"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "face-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('face')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cor.face")
### * cor.face

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cor.face
### Title: Extraction of correlation and mean from a 'face.sparse' object
### Aliases: cor.face
### Keywords: ~face.sparse

### ** Examples

# See the examples for "face.sparse".



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cor.face", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("face.sparse")
### * face.sparse

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: face.sparse
### Title: Fast covariance estimation for sparse functional data
### Aliases: face.sparse
### Keywords: ~face.sparse

### ** Examples


## Not run: 
##D ##########################
##D #### CD4 data example
##D ##########################
##D 
##D require(refund)
##D data(cd4)
##D n <- nrow(cd4)
##D T <- ncol(cd4)
##D 
##D id <- rep(1:n,each=T)
##D t <- rep(-18:42,times=n)
##D y <- as.vector(t(cd4))
##D sel <- which(is.na(y))
##D 
##D 
##D ## organize data and apply FACEs
##D data <- data.frame(y=log(y[-sel]),
##D argvals = t[-sel],
##D subj = id[-sel])
##D data <- data[data$y>4.5,]
##D fit_face <- face.sparse(data,argvals.new=(-20:40))
##D 
##D data.h <- data
##D tnew <- fit_face$argvals.new
##D 
##D ## scatter plots
##D Xlab <- "Months since seroconversion"
##D Ylab <- "log (CD4 count)"
##D par(mfrow=c(1,1),mar = c(4.5,4.5,3,2))
##D id <- data.h$subj
##D uid <- unique(id)
##D plot(data.h$argvals,data.h$y,
##D type = "n", ylim = c(4.5,8),
##D xlab = Xlab, ylab = Ylab,
##D cex.lab = 1.25,cex.axis=1.25,cex.main = 1.25)
##D 
##D for(i in 1:10){
##D seq <- which(id==uid[i])
##D lines(data.h$argvals[seq],data.h$y[seq],lty=1,col="gray",lwd=1,type="l")
##D #points(data.h$argvals[seq],data.h$y[seq],col=1,lty=1,pch=1)
##D }
##D 
##D Sample <- seq(10,50,by=10)
##D for(i in Sample){
##D seq <- which(id==uid[i])
##D lines(data.h$argvals[seq],data.h$y[seq],lty=1,col="black",lwd=1,type="l")
##D }
##D lines(tnew,fit_face$mu.new,lwd=2,lty=2,col="red")
##D 
##D ## plots of variance/correlation functions
##D 
##D Cov <- fit_face$Chat.new
##D Cov_diag <- diag(Cov)
##D Cor <- fit_face$Cor.new
##D 
##D par(mfrow=c(1,2),mar=c(4.5,4.1,3,4.5))
##D 
##D 
##D plot(tnew,Cov_diag,type="l",
##D xlab = Xlab, ylab="",main= "CD4: variance function",
##D #ylim = c(0.8,1.5),
##D cex.axis=1.25,cex.lab=1.25,cex.main=1.25,lwd=2)
##D 
##D require(fields)
##D image.plot(tnew,tnew,Cor,
##D xlab=Xlab, ylab = Xlab,
##D main = "CD4: correlation function",
##D cex.axis=1.25,cex.lab=1.25,cex.main=1.25,
##D axis.args = list(at = c(0,0.2,0.4,0.6,0.8,1.0)),
##D legend.shrink=0.75,legend.line=-1.5)
##D 
##D 
##D ## prediction of several subjects
##D 
##D par(mfrow=c(2,2),mar=c(4.5,4.5,3,2))
##D Sample <- c(30,40,50,60)
##D for(i in 1:4){
##D sel <- which(id==uid[Sample[i]])
##D dati <- data.h[sel,]
##D 
##D seq <- -20:40
##D k <- length(seq)
##D dati_pred <- data.frame(y = rep(NA,nrow(dati) + k ),
##D argvals = c(rep(NA,nrow(dati)),seq),
##D subj=rep(dati$subj[1],nrow(dati) + k )
##D )
##D 
##D dati_pred[1:nrow(dati),] <- dati
##D yhat2 <- predict(fit_face,dati_pred)
##D 
##D data3 <- dati
##D Ylim <- range(c(data3$y,yhat2$y.pred))
##D 
##D plot(data3$argvals,data3$y,xlab=Xlab,ylab=Ylab, main = paste("Male ",i,sep=""),
##D ylim = c(4,8.5),
##D cex.lab=1.25,cex.axis = 1.25,cex.main = 1.25,pch=1,xlim=c(-20,40))
##D 
##D Ord <- nrow(dati) + 1:k
##D lines(dati_pred$argvals[Ord],yhat2$y.pred[Ord],col="red",lwd=2)
##D lines(dati_pred$argvals[Ord],
##D yhat2$y.pred[Ord] - 1.96*yhat2$se.pred[Ord], col="red",lwd=1,lty=2)
##D lines(dati_pred$argvals[Ord],
##D yhat2$y.pred[Ord] + 1.96*yhat2$se.pred[Ord], col="red",lwd=1,lty=2)
##D 
##D lines(tnew,fit_face$mu.new,lty=3,col="black",lwd=2)
##D legend("bottomleft",c("mean","prediction"),lty=c(3,1),col=1:2,lwd=2,bty="n")
##D }
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("face.sparse", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.face.sparse")
### * predict.face.sparse

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.face.sparse
### Title: Subject-specific curve prediction from a face.sparse fit
### Aliases: predict.face.sparse
### Keywords: ~face.sparse ~prediction

### ** Examples


#See the examples for "face.sparse".



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.face.sparse", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.pspline.face")
### * predict.pspline.face

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.pspline.face
### Title: Mean prediction from a P-spline smoothing fit
### Aliases: predict.pspline.face
### Keywords: ~Pspline ~prediction

### ** Examples

#See the examples for "pspline".



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.pspline.face", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pspline")
### * pspline

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pspline
### Title: Univariate P-spline smoothing
### Aliases: pspline
### Keywords: ~Pspline

### ** Examples

## cd4 data
require(refund)
data(cd4)
n <- nrow(cd4)
T <- ncol(cd4)

id <- rep(1:n,each=T)
t <- rep(-18:42,times=n)
y <- as.vector(t(cd4))
sel <- which(is.na(y))

## organize data
data <- data.frame(y=log(y[-sel]),
argvals = t[-sel],
subj = id[-sel])
data <- data[data$y>4.5,]

## smooth
fit <- pspline(data)

## plot
plot(data$argvals,fit$mu.new,type="p")

## prediction
pred <- predict(fit,quantile(data$argvals,c(0.2,0.6)))
pred



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pspline", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("select.knots")
### * select.knots

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: select.knots
### Title: Knots selection for P-spline smoothing
### Aliases: select.knots
### Keywords: ~Pspline

### ** Examples

t <- rnorm(100)
knots <- select.knots(t)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("select.knots", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
