\name{fpca.sc}
\alias{fpca.sc}
\title{
Functional principal components analysis by smoothed covariance
}
\description{
Decomposes functional observations using functional principal components analysis. A mixed model framework is used to estimate scores and obtain variance estimates.
}
\usage{
fpca.sc(Y, Y.pred=NULL, nbasis = 10, pve = .99, npc = NULL, var = FALSE, 
        simul = FALSE, sim.alpha = .95, useSymm = FALSE, makePD = FALSE)
}
\arguments{
  \item{Y }{matrix of observed functions for which estimates and covariance matrices are desired.}
  \item{Y.pred }{if desired, a matrix of observed functions that will be estimated using the FPC decomposition of Y.}
  \item{nbasis }{number of splines used in the estimation of the mean function and the bivariate smoothing of the covariance matrix.}
  \item{pve }{proportion of variance explained used to choose the number of principal components to be included in the expansion.}
  \item{npc }{prespecified value for the number of principal components to be included in the expansion (if given, this overrides 'pve').}
  \item{var }{\code{TRUE} or \code{FALSE} indicating whether model-based estimates for the variance of FPCA expansions should be computed.}
  \item{simul }{\code{TRUE} or \code{FALSE}, indicating whether critical values for simultaneous confidence intervals should be estimated.}
  \item{sim.alpha }{alpha level of the simultaneous intervals.}
  \item{useSymm }{\code{TRUE} or \code{FALSE}, indicating whether to do the smoothing based on the upper triagonal of the empiricial covariance or not.
      Can save computation time for large data and allows for covariance surfaces that are very peaked on the diagonal.}
  \item{makePD }{\code{TRUE} or \code{FALSE}, indicating whether to enforce positive definiteness of the estimated surface.}
}
\details{
  This function computes a FPC decomposition for a collection of observed curves. The curves may be sparsely observed or measured with error. After the decomposition is estimated, a mixed model framework is used to estimate curve-specific scores and variances.

FPCA via kernel smoothing of the covariance function, with the diagonal treated separately, was proposed in Staniswalis and Lee (1998) and much extended by Yao et al. (2005), who introduced the "PACE" method. \code{fpca.sc} uses penalized splines to smooth the covariance function, as developed by Di et al. (2009) and Goldsmith et al. (2012).
}
\value{
  \item{Yhat }{\eqn{n \times d} matrix whose rows are the estimates (truncated FPC expansions) of the curves in Y.}
  \item{scores }{\eqn{n \times npc} matrix of estimated principal component scores.}
  \item{mu }{estimated mean function.}
  \item{efunctions }{\eqn{d \times npc} matrix of estimated eigenfunctions of the functional covariance operator, i.e., the FPC basis functions.}
  \item{evalues }{estimated eigenvalues of the covariance operator, i.e., variances of FPC scores.}
  \item{npc }{number of FPCs: either the supplied \code{npc}, or the minimum number of basis functions needed to explain proportion \code{pve} of the variance in the observed curves.}
  \item{sigma2 }{estimated measurement error variance.}
  \item{diag.var }{the diagonal elements of the covariance matrices for each estimated curve.}
  \item{VarMats }{a list containing the estimated covariance matrices for each curve in Y.}
  \item{crit.val }{estimated critical values for constructing simultaneous confidence intervals.}
}
\references{
Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009). Multilevel functional principal component analysis. \emph{Annals of Applied Statistics}, 3, 458--488.

Goldsmith, J., Greven, S., and Crainiceanu, C. (2012). Corrected confidence bands for functional data using principal components. Submitted.

Staniswalis, J. G., and Lee, J. J. (1998). Nonparametric regression analysis of longitudinal data. \emph{Journal of the American Statistical Association}, 93, 1403--1418.

Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis for sparse longitudinal data. \emph{Journal of the American Statistical Association}, 100, 577--590.
}
\author{
Jeff Goldsmith \email{jeff.goldsmith@columbia.edu}, Sonja Greven \email{sonja.greven@stat.uni-muenchen.de}, and Lan Huo \email{Lan.Huo@nyumc.org}
}
\examples{
\dontrun{
data(cd4)

Fit.MM = fpca.sc(cd4, var = TRUE, simul = TRUE)

# for one subject, examine curve estimate, pointwise and simultaneous itervals
EX = 1
EX.MM =  cbind(Fit.MM$Yhat[EX,], 
      Fit.MM$Yhat[EX,] + 1.96 * sqrt(Fit.MM$diag.var[EX,]), 
      Fit.MM$Yhat[EX,] - 1.96 * sqrt(Fit.MM$diag.var[EX,]),
      Fit.MM$Yhat[EX,] + Fit.MM$crit.val[EX] * sqrt(Fit.MM$diag.var[EX,]),
      Fit.MM$Yhat[EX,] - Fit.MM$crit.val[EX] * sqrt(Fit.MM$diag.var[EX,]))


par(mfrow=c(1,3))

# plot data for one subject, with curve and interval estimates
d = as.numeric(colnames(cd4))
plot(d[which(!is.na(cd4[EX,]))], cd4[EX,which(!is.na(cd4[EX,]))], type = 'o', pch = 19,
  cex=.75, ylim = range(0, 3400), xlim = range(d), xlab = "Months since seroconversion", 
    lwd = 1.2, ylab = "Total CD4 Cell Count", main = "Est. & CI - Sampled Data")

matpoints(d, EX.MM, col = 4, type = 'l', lwd = c(2, 1, 1, 1, 1), lty = c(1,1,1,2,2))

# plot estimated mean function
plot(d, Fit.MM$mu, type = 'l', xlab = "Months since seroconversion",
  ylim = range(0, 3400), ylab = "Total CD4 Cell Count", main = "Est. Mean Function")

# plot the first estimated basis function
plot(d, Fit.MM$efunctions[,1], type = 'l', xlab = "Months since seroconversion",
  ylab = "Total CD4 Cell Count", main = "First Est. Basis Function")
}}