\name{fpca.ssvd}
\alias{fpca.ssvd}
\title{Smoothed functional principal component analysis via iterative penalized rank one SVDs}
\usage{
  fpca.ssvd(X, npc = NA, center = TRUE, maxiter = 15,
    tol = 1e-04, diffpen = 3, gridsearch = TRUE,
    alphagrid = 1.5^(-20:40), lower.alpha = 1e-05,
    upper.alpha = 1e+07, verbose = FALSE)
}
\arguments{
  \item{X}{data matrix (rows: observations; columns: grid
  of eval. points)}

  \item{npc}{how many smooth SVs to try to extract, if \code{NA}
  (the default) the hard thresholding rule of Donoho and
  Gavish (2013) is used (see Details, References).}

  \item{center}{center \code{X} so that its column-means
  are 0? Defaults to \code{TRUE}}

  \item{maxiter}{how many iterations of the power algorithm
  to perform at most (defaults to 15).}

  \item{tol}{convergence tolerance for power algorithm
  (defaults to 1e-4).}

  \item{diffpen}{difference penalty order controlling the
  desired smoothness of the right singular vectors,
  defaults to 3 (i.e., deviations from local quadratic
  polynomials).}

  \item{gridsearch}{use \code{\link[stats]{optimize}} or a
  grid search to find GCV-optimal smoothing parameters?
  Defaults to \code{TRUE}.}

  \item{alphagrid}{grid of smoothing parameter values for
  grid search}

  \item{lower.alpha}{lower limit for for smoothing
  parameter if \code{!gridsearch}}

  \item{upper.alpha}{upper limit for smoothing parameter if
  \code{!gridsearch}}

  \item{verbose}{generate graphical summary of progress and
  diagnostic messages?  defaults to \code{FALSE}}
}
\value{
  a list with entries \code{smooth}, the rank \code{npc}
  SVD of the smooth part of \code{X}, and \code{noise}, the
  rank \code{min(dim(X)) - npc} SVD of the non-smooth part
  of \code{X}, and \code{mean}, the column means of
  \code{X} (or a vector of zeroes if \code{!center}).
}
\description{
  Implements the algorithm of Huang, Shen and Buja (2008) for finding smooth right singular vectors of
  a matrix \code{X} containing (contaminated) evaluations
  of functional random variables on a regular, equidistant
  grid. If the number of smooth SVs to extract is not
  specified, the function hazards a guess for the
  appropriate number based on the asymptotically optimal
  truncation threshold under the assumption of a low rank
  matrix contaminated with i.i.d. Gaussian noise with
  unknown variance derived in Donoho and Gavish (2013).
}
\examples{
## as in Sec. 6.2 of Huang, Shen and Buja (2008):
 set.seed(2678695)
 n <- 101
 m <- 101
 s1 <- 20
 s2 <- 10
 s <- 4
 t <- seq(-1, 1, l=m)
 v1 <- t + sin(pi*t)
 v2 <- cos(3*pi*t)
 V <- cbind(v1/sqrt(sum(v1^2)), v2/sqrt(sum(v2^2)))
 U <- matrix(rnorm(n*2), n, 2)
 D <- diag(c(s1^2, s2^2))
 eps <- matrix(rnorm(m*n, sd=s), n, m)
 X <- U\%*\%D\%*\%t(V) + eps

 smoothSV <- fpca.ssvd(X, verbose=TRUE)

 layout(t(matrix(1:6, nr=2)))
 clrs <- sapply(rainbow(n), function(c)
           do.call(rgb, as.list(c(col2rgb(c)/255, .1))))
 matplot(V, type="l", lty=1, col=1:2, xlab="",
         main="V: true", bty="n")
 matplot(smoothSV$smooth$v, type="l", lty=1, col=1:5, xlab="",
         main="V: estimate", bty="n")
 plot(U[,1]*s1^2, smoothSV$smooth$u[,1]*smoothSV$smooth$d[1], col=1,
         xlab="true U", ylab="estimated U",
         main="U[,1]: true vs. estimated", bty="n"); abline(c(0,-1), col="grey")
 plot(U[,2]*s2^2, smoothSV$smooth$u[,2]*smoothSV$smooth$d[2], col=1,
         xlab="true U", ylab="estimated U",
         main="U[,2]: true vs. estimated", bty="n"); abline(c(0,-1), col="grey")
 matplot(1:m, t(U\%*\%D\%*\%t(V)), type="l", lty=1, col=clrs, xlab="", ylab="",
         main="true smooth X", bty="n")
 matplot(1:m, with(smoothSV,
                 mean + t(smooth$u\%*\%diag(smooth$d)\%*\%t(smooth$v))),
         xlab="", ylab="",
         type="l", lty=1,col=clrs, main="estimated smooth X", bty="n")
}
\author{
  Fabian Scheipl
}
\references{
  Huang, J. Z., Shen, H., & Buja, A. (2008).  Functional
  principal components analysis via penalized rank one
  approximation.  \emph{Electronic Journal of Statistics}, 2,
  678--695.

  Donoho, D.L., and Gavish, M. (2013). The optimal hard
  threshold for singular values is 4/sqrt(3).  eprint
  arXiv:1305.5870. Available from
  \url{http://arxiv.org/abs/1305.5870}.
}

