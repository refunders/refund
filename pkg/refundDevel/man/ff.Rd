\name{ff}
\alias{ff}
\title{Construct a function-on-function regression term}
\usage{
ff(X, yind = NULL, xind = seq(0, 1, l = ncol(X)), basistype = c("te",
  "t2", "s"), integration = c("simpson", "trapezoidal", "riemann"),
  L = NULL, limits = NULL, splinepars = list(bs = "ps", m = list(c(2, 1),
  c(2, 1))), check.ident = TRUE)
}
\arguments{
  \item{X}{an n by \code{ncol(xind)} matrix of function
  evaluations \eqn{X_i(s_{i1}),\dots, X_i(s_{iS})};
  \eqn{i=1,\dots,n}.}

  \item{yind}{\emph{DEPRECATED} used to supply matrix (or
  vector) of indices of evaluations of \eqn{Y_i(t)}, no
  longer used.}

  \item{xind}{matrix (or vector) of indices of evaluations
  of \eqn{X_i(s)}; i.e. matrix with rows
  \eqn{(s_{i1},\dots,s_{iS})}}

  \item{basistype}{defaults to "\code{\link[mgcv]{te}}",
  i.e. a tensor product spline to represent
  \eqn{\beta(t,s)}.  Alternatively, use \code{"s"} for
  bivariate basis functions (see \code{mgcv}'s
  \code{\link[mgcv]{s}}) or \code{"t2"} for an alternative
  parameterization of tensor product splines (see
  \code{mgcv}'s \code{\link[mgcv]{t2}}).}

  \item{integration}{method used for numerical integration.
  Defaults to \code{"simpson"}'s rule for calculating
  entries in \code{L}.  Alternatively and for
  non-equidistant grids, \code{"trapezoidal"} or
  \code{"riemann"}. \code{"riemann"} integration is always
  used if \code{limits} is specified}

  \item{L}{optional: an n by \code{ncol(xind)} matrix
  giving the weights for the numerical integration over
  \eqn{s}.}

  \item{limits}{defaults to NULL for integration across the
  entire range of \eqn{X(s)}, otherwise specifies the
  integration limits \eqn{s_{hi, i}, s_{lo, i}}: either one
  of \code{"s<t"} or \code{"s<=t"} for \eqn{(s_{hi, i},
  s_{lo, i}) = (0, t)} or a function that takes \code{s} as
  the first and \code{t} as the second argument and returns
  TRUE for combinations of values \code{(s,t)} if \code{s}
  falls into the integration range for the given \code{t}.
  This is an experimental feature and not well tested yet;
  use at your own risk.}

  \item{splinepars}{optional arguments supplied to the
  \code{basistype}-term. Defaults to a cubic tensor product
  B-spline with marginal first difference penalties, i.e.
  \code{list(bs="ps", m=list(c(2, 1), c(2,1)))}. See
  \code{\link[mgcv]{te}} or \code{\link[mgcv]{s}} in
  \pkg{mgcv} for details}

  \item{check.ident}{check identifiability of the model
  spec. See Details and References. Defaults to TRUE.}
}
\value{
a list containing \itemize{ \item \code{call} a "call" to
\code{\link[mgcv]{te}} (or \code{\link[mgcv]{s}},
\code{\link[mgcv]{t2}}) using the appropriately constructed
covariate and weight matrices \item \code{data} a list
containing the necessary covariate and weight matrices }
}
\description{
Defines a term \eqn{\int^{s_{hi, i}}_{s_{lo, i}}
X_i(s)\beta(t,s)ds} for inclusion in an
\code{mgcv::gam}-formula (or \code{bam} or \code{gamm} or
\code{gamm4:::gamm4}) as constructed by \code{\link{pffr}}.
\cr Defaults to a cubic tensor product B-spline with
marginal first order differences penalties for
\eqn{\beta(t,s)} and numerical integration over the entire
range \eqn{[s_{lo, i}, s_{hi, i}] = [\min(s_i), \max(s_i)]}
by using Simpson weights. Can't deal with any missing
\eqn{X(s)}, unequal lengths of \eqn{X_i(s)} not (yet?)
possible. Unequal integration ranges for different
\eqn{X_i(s)} should work. \eqn{X_i(s)} is assumed to be
numeric (duh...).
}
\details{
If \code{check.ident==TRUE} and \code{basistype!="s"} (the
default), the routine checks for overlap between the
kernels of Cov\eqn{(X(s))} and the marginal penalty over
\code{s}, as well as for sufficient rank of
Cov\eqn{(X(s))}. If there is overlap of the kernels, the
penalty and basis for the covariate direction are changed
to B-splines with a modified `shrinkage' penalty, see
\code{\link{smooth.construct.pss.smooth.spec}}. A warning
is given if the effective rank of Cov\eqn{(X(s))} (defined
the number of eigenvalues accounting for at least 0.995 of
the total variance in \eqn{X_i(s)}) is lower than the
dimension of the basis for the covariate direction. See
reference for details. Using an \code{\link{ffpc}}-term may
be preferable if \eqn{X_i(s)} is of very low rank.
}
\author{
Fabian Scheipl, Sonja Greven
}
\references{
For background on \code{check.ident}:\cr Scheipl, F., &
Greven, S. (2012). Identifiability in penalized
function-on-function regression models.  LMU Munich,
Department of Statistics: Technical Report 125.
\url{http://epub.ub.uni-muenchen.de/13060/}
}
\seealso{
\code{mgcv}'s \code{\link[mgcv]{linear.functional.terms}}
}

