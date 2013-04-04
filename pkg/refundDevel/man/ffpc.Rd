\name{ffpc}
\alias{ffpc}
\title{Construct a PC-based function-on-function regression term}
\usage{
  ffpc(X, yind, xind = seq(0, 1, length = ncol(X)),
    splinepars = list(bs = "ps", m = c(2, 1), k = 8),
    center = TRUE,
    decomppars = list(pve = 0.99, useSymm = TRUE),
    npc.max = 15)
}
\arguments{
  \item{X}{an n by \code{ncol(xind)} matrix of function
  evaluations \eqn{X_i(s_{i1}),\dots, X_i(s_{iS})};
  \eqn{i=1,\dots,n}.}

  \item{yind}{matrix (or vector) of indices of evaluations
  of \eqn{Y_i(t)}}

  \item{xind}{matrix (or vector) of indices of evaluations
  of \eqn{X_i(t)}, defaults to \code{seq(0, 1,
  length=ncol(X))}.}

  \item{splinepars}{optional arguments supplied to the
  \code{basistype}-term. Defaults to a cubic B-spline with
  first difference penalties and 8 basis functions for each
  \eqn{\tilde \beta_k(t)}.}

  \item{center}{center \code{X} so that the mean of
  \eqn{X(s)} is zero for each index value \eqn{s}?}

  \item{decomppars}{paramters for the FPCA performed with
  \code{\link{fpca.sc}}.}

  \item{npc.max}{maximal number \eqn{K} of FPCs to use,
  regardless of \code{decomppars}; defaults to 15}
}
\value{
  a list containing the necessary information to construct
  a term to be included in a \code{mgcv::gam}-formula.
}
\description{
  Defines a term \eqn{\int X_i(s)\beta(t,s)ds} for
  inclusion in an \code{mgcv::gam}-formula (or \code{bam}
  or \code{gamm} or \code{gamm4:::gamm}) as constructed by
  \code{\link{pffr}}. In contrast to \code{\link{ff}},
  \code{ffpc} does an FPCA decomposition \eqn{X(s) \approx
  \sum^K_{k=1} \xi_{ik} \Phi_k(s)} using
  \code{\link{fpca.sc}}. Because \deqn{\int
  X_i(s)\beta(t,s)ds = \sum^K_{k=1} \xi_{ik} \int \Phi_k(s)
  \beta(s,t) ds = \sum^K_{k=1} \xi_{ik} \tilde \beta_k(t),}
  i.e., the function-on-function term can be represented as
  a sum of \eqn{K} univariate functions \eqn{\tilde
  \beta_k(t)} in \eqn{t} each multiplied by the FPC
  loadings \eqn{\xi_{ik}}. The truncation paramter \eqn{K}
  is chosen as described in \code{\link{fpca.sc}}. To
  reduce model complexity, the \eqn{\tilde \beta_k(t)} all
  have a single joint smoothing parameter (in \code{mgcv},
  they get the same \code{id}).\cr
}
\details{
  Using this instead of \code{ff()} can be beneficial if
  the covariance operator of the \eqn{X_i(s)} has low
  effective rank (i.e., if \eqn{K} is small). If the
  covariance operator of the \eqn{X_i(s)} is of (very) high
  rank, i.e., if eqn{K} is large, \code{ffpc()} will not be
  very efficient.

  Please see \code{\link[refund]{pffr}} for details on
  model specification and implementation. \cr \code{ffpc()}
  IS AN EXPERIMENTAL FEATURE AND NOT WELL TESTED YET -- USE
  AT YOUR OWN RISK.
}
\examples{
\dontrun{
set.seed(1122)
n <- 55
S <- 60
T <- 50
s <- seq(0,1, l=S)
t <- seq(0,1, l=T)
df <- 10
B <- bs(s, df=df)
X <- t(B \%*\% matrix(rt(n*df, df=3), nrow=df))
beta.st <- outer(s,t, test2)

y <- (1/S*X)\%*\%beta.st

data <- list(y=y, X=X)
# set number of FPCs to true rank of process for this example:
m.pc <- pffr(y ~ c(1) + 0 + ffpc(X, yind=t, decomppars=list(npc=df)),
        data=data, yind=t)
summary(m.pc)
m.ff <- pffr(y ~ c(1) + 0 + ff(X, yind=t), data=data, yind=t)
summary(m.ff)

# plot implied coefficient surfaces:
betatilde <- predict.gam(m.pc, newdata=data.frame(t.vec=t, X.PC1=1, X.PC2=1,
        X.PC3=1, X.PC4=1, X.PC5=1, X.PC6=1, X.PC7=1, X.PC8=1, X.PC9=1,
        X.PC10=1, X.PC10=1), type="terms")
layout(t(1:3))
persp(t, s, t(beta.st), theta=30, phi=40, main="Truth")
plot(m.ff, select=1, pers=TRUE, zlim=range(beta.st), theta=30, phi=40,
        main="ff()")
# betatilde_k (t) = \\int Phi_k(s) beta(s,t) ds \\approx 1/S t(Phi_k) \%*\% beta(s,t)
# --> beta(s,t) = S * Phi \%*\% [betatilde_1(t) ... betatilde_K(t)]
persp(t, s, t(S * (m.pc$pffr$ffpc[[1]]$PCMat \%*\% t(betatilde))), zlim=range(beta.st),
   theta=30, phi=40, main="ffpc()")
}
}
\author{
  Fabian Scheipl
}

