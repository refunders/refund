\name{face.sparse}
\alias{face.sparse}

\title{
Fast covariance estimation for sparse functional data
}
\description{
The function is to estimate the mean and covariance function
from a cluster of functions/longitudinal observations.}
\usage{
face.sparse(data, newdata = NULL,
             center = TRUE, argvals.new = NULL,
             knots = 7, 
             p = 3, m = 2, lambda = NULL, lambda_mean = NULL,
             search.length = 14, 
             lower = -3, upper = 10, lower2 = -3, upper2 = 5,
             calculate.scores = FALSE,pve=0.99,two_step=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data frame with three arguments:
(1) \code{argvals}: observation times;
(2) \code{subj}: subject indices;
(3) \code{y}: values of observations.
Missing values not allowed.
}
  \item{newdata}{
of the same strucutre as \code{data}; defaults to NULL, then no prediction. 
}
  \item{center}{
logical. If TRUE, then Pspline smoothing of the population mean will be conducted and
subtracted from the data before covariance smoothing; if FALSE, then
the population mean will be just 0s.
}
  \item{argvals.new}{
a vector of observation time points to evaluate mean function, covariance function, error variance
and etc. If NULL,
then 100 equidistant points in the range of \code{data$argvals}.
}
  \item{knots}{
 the number of knots for B-spline basis functions to be used; defaults to 7.
    The resulting number of basis functions is the number of interior knots plus the degree of B-splines.
}


\item{p}{
the degrees of B-splines; defaults to 3.
}
\item{m}{
the order of differencing penalty; defaults to 2.
}
\item{lambda}{
the value of the smoothing parameter for covariance smoothing; defaults to NULL.
}
  \item{lambda_mean}{
the value of the smoothing parameter for mean smoothing; defaults to NULL.
}

\item{search.length}{
the number of equidistant (log scale) smoothing parameters to search; defaults to 14.
}
\item{lower, upper}{
bounds for log smoothing parameter for first step of estimation; defaults are -3 and 10, respectively.
}
\item{lower2, upper2}{
bounds for log smoothing parameter for second step of estimation; defaults are \code{lower} and 5, respectively.
}
\item{calculate.scores}{
if TRUE, scores will be calculated.
}
\item{pve}{
Defaults 0.99. To select the number of eigenvalues by percentage of variance.
}
\item{two_step}{
if TRUE, a two-step estimation procedure will be applied.
}
}
\details{
This is a generalized version of bivariate P-splines (Eilers and Marx, 2003) for
covariance smoothing of sparse functional or longitudinal data.
It uses tensor product B-spline basis functions and
employes a differencing penalty on the assosciated parameter matrix.
The only smoothing parameter in the method is selected
by leave-one-subject-out cross validation and is implemented with a fast algorithm.

There are two steps for estimation. During the first step, the objective function to minimize
is the penalized least squares on empirical estimates of covariance function. During the second step,
the covariance between the empirical estimates (depending on the estimates of 
covariance function) are accounted and thus a generalized penalized
least squares are minimized.

If \code{center} is TRUE, then a population mean will be calculated and is smoothed by
univariate P-spline smoothing:\code{pspline} (Eilers and Marx, 1996). This univariate
smoothing uses leave-one-subject-out cross validation to select the smoothing parameter.

The knots are "equally-spaced", the differencing penalty in Eilers and Marx (2003)
is used.
                                                
If the functional data are observed at the same grid for each function/curve
and can be organized into
a data matrix, then \code{fpca.face} in the package \code{refund}
should instead be used. \code{fpca.face} allows a small percentage (less than
30 percent) of 
missing data in the data matrix.
}
\value{
\item{newdata}{Input}
\item{y.pred,mu.pred,Chat.diag.pred, var.error.pred}{Predicted/estimated objects at
 \code{newdata$argvals}}
\item{Theta}{Estimated parameter matrix}
\item{argvals.new}{Vector of time points to evaluate population parameters}
\item{mu.new, Chat.new,Cor.new,Cor.raw.new,Chat.raw.diag.new, var.error.new}{Estimated objects at
\code{argvals.new}}
\item{eigenfunctions, eigenvalues}{Estimated eigenfunctions (scaled eigenvector) and eigenvalues at \code{argvals.new}}
\item{mu.hat,var.error.hat}{Estimated objects at  \code{data$argvals}}
\item{calculate.scores,scores}{if \code{calculate.scores} is TRUE (default to FALSE) and \code{newdata}
is not NULL, then predicted scores will be calculated.}
\item{\dots}{\dots}
}
\references{
Luo Xiao, Cai Li, William Checkley and Ciprian Crainiceanu, Fast covariance estimation
for sparse functional data, Stat. Comput., \doi{doi:10.1007/s11222-017-9744-8}.

Paul Eilers and Brian  Marx, Multivariate calibration with temperature
interaction using two-dimensional penalized signal regression,
Chemometrics and Intelligent Laboratory Systems 66 (2003), 159-174.

Paul Eilers and Brian Marx, Flexible smoothing with B-splines and penalties,
Statist. Sci., 11, 89-121, 1996.
 
Simon N. Wood, P-splines with derivative based penalties and tensor product 
  smoothing of unevenly distributed data, Stat. Comput., \doi{10.1007/s11222-016-9666-x}.
}
\author{
Luo Xiao <lxiao5@ncsu.edu> and Cai Li <cli9@ncsu.edu>
}

\seealso{
\code{fpca.face} and \code{fpca.sc} in \code{refund}
}
\examples{

\dontrun{
##########################
#### CD4 data example
##########################

require(refund)
data(cd4)
n <- nrow(cd4)
T <- ncol(cd4)

id <- rep(1:n,each=T)
t <- rep(-18:42,times=n)
y <- as.vector(t(cd4))
sel <- which(is.na(y))


## organize data and apply FACEs
data <- data.frame(y=log(y[-sel]),
argvals = t[-sel],
subj = id[-sel])
data <- data[data$y>4.5,]
fit_face <- face.sparse(data,argvals.new=(-20:40))

data.h <- data
tnew <- fit_face$argvals.new

## scatter plots
Xlab <- "Months since seroconversion"
Ylab <- "log (CD4 count)"
par(mfrow=c(1,1),mar = c(4.5,4.5,3,2))
id <- data.h$subj
uid <- unique(id)
plot(data.h$argvals,data.h$y,
type = "n", ylim = c(4.5,8),
xlab = Xlab, ylab = Ylab,
cex.lab = 1.25,cex.axis=1.25,cex.main = 1.25)

for(i in 1:10){
seq <- which(id==uid[i])
lines(data.h$argvals[seq],data.h$y[seq],lty=1,col="gray",lwd=1,type="l")
#points(data.h$argvals[seq],data.h$y[seq],col=1,lty=1,pch=1)
}

Sample <- seq(10,50,by=10)
for(i in Sample){
seq <- which(id==uid[i])
lines(data.h$argvals[seq],data.h$y[seq],lty=1,col="black",lwd=1,type="l")
}
lines(tnew,fit_face$mu.new,lwd=2,lty=2,col="red")

## plots of variance/correlation functions

Cov <- fit_face$Chat.new
Cov_diag <- diag(Cov)
Cor <- fit_face$Cor.new

par(mfrow=c(1,2),mar=c(4.5,4.1,3,4.5))


plot(tnew,Cov_diag,type="l",
xlab = Xlab, ylab="",main= "CD4: variance function",
#ylim = c(0.8,1.5),
cex.axis=1.25,cex.lab=1.25,cex.main=1.25,lwd=2)

require(fields)
image.plot(tnew,tnew,Cor,
xlab=Xlab, ylab = Xlab,
main = "CD4: correlation function",
cex.axis=1.25,cex.lab=1.25,cex.main=1.25,
axis.args = list(at = c(0,0.2,0.4,0.6,0.8,1.0)),
legend.shrink=0.75,legend.line=-1.5)


## prediction of several subjects

par(mfrow=c(2,2),mar=c(4.5,4.5,3,2))
Sample <- c(30,40,50,60)
for(i in 1:4){
sel <- which(id==uid[Sample[i]])
dati <- data.h[sel,]

seq <- -20:40
k <- length(seq)
dati_pred <- data.frame(y = rep(NA,nrow(dati) + k ),
argvals = c(rep(NA,nrow(dati)),seq),
subj=rep(dati$subj[1],nrow(dati) + k )
)

dati_pred[1:nrow(dati),] <- dati
yhat2 <- predict(fit_face,dati_pred)

data3 <- dati
Ylim <- range(c(data3$y,yhat2$y.pred))

plot(data3$argvals,data3$y,xlab=Xlab,ylab=Ylab, main = paste("Male ",i,sep=""),
ylim = c(4,8.5),
cex.lab=1.25,cex.axis = 1.25,cex.main = 1.25,pch=1,xlim=c(-20,40))

Ord <- nrow(dati) + 1:k
lines(dati_pred$argvals[Ord],yhat2$y.pred[Ord],col="red",lwd=2)
lines(dati_pred$argvals[Ord],
yhat2$y.pred[Ord] - 1.96*yhat2$se.pred[Ord], col="red",lwd=1,lty=2)
lines(dati_pred$argvals[Ord],
yhat2$y.pred[Ord] + 1.96*yhat2$se.pred[Ord], col="red",lwd=1,lty=2)

lines(tnew,fit_face$mu.new,lty=3,col="black",lwd=2)
legend("bottomleft",c("mean","prediction"),lty=c(3,1),col=1:2,lwd=2,bty="n")
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~face.sparse}
