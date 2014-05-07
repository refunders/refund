\name{coef.pffr}
\alias{coef.pffr}
\title{Get estimated coefficients from a pffr fit}
\usage{
\method{coef}{pffr}(object, raw = FALSE, se = TRUE, freq = FALSE,
  sandwich = FALSE, seWithMean = TRUE, n1 = 100, n2 = 40, n3 = 20,
  Ktt = NULL, ...)
}
\arguments{
  \item{object}{a fitted \code{pffr}-object}

  \item{raw}{logical, defaults to FALSE. If TRUE, the
  function simply returns \code{object$coefficients}}

  \item{se}{logical, defaults to TRUE. Return estimated
  standard error of the estimates?}

  \item{freq}{logical, defaults to FALSE. If FALSE, use
  posterior variance \code{object$Vp} for variability
  estimates, else use \code{object$Ve}. See
  \code{\link[mgcv]{gamObject}}}

  \item{sandwich}{logical, defaults to FALSE. Use a
  Sandwich-estimator for approximate variances? See
  Details.  THIS IS AN EXPERIMENTAL FEATURE, USE A YOUR OWN
  RISK.}

  \item{seWithMean}{logical, defaults to TRUE. Include
  uncertainty about the intercept/overall mean in standard
  errors returned for smooth components?}

  \item{n1}{see below}

  \item{n2}{see below}

  \item{n3}{\code{n1, n2, n3} give the number of gridpoints
  for 1-/2-/3-dimensional smooth terms used in the marginal
  equidistant grids over the range of the covariates at
  which the estimated effects are evaluated.}

  \item{Ktt}{(optional) an estimate of the covariance
  operator of the residual process \eqn{\epsilon_i(t) \sim
  N(0, K(t,t'))}, evaluated on \code{yind} of
  \code{object}. If not supplied, this is estimated from
  the crossproduct matrices of the observed residual
  vectors. Only relevant for sandwich CIs.}

  \item{...}{other arguments, not used.}
}
\value{
If \code{raw==FALSE}, a list containing \itemize{ \item
\code{pterms} a matrix containing the parametric /
non-functional coefficients (and, optionally, their se's)
\item \code{smterms} a named list with one entry for each
smooth term in the model. Each entry contains \itemize{
\item \code{coef} a matrix giving the grid values over the
covariates, the estimated effect (and, optionally, the
se's).  The first covariate varies the fastest.  \item
\code{x, y, z} the unique gridpoints used to evaluate the
smooth/coefficient function/coefficient surface \item
\code{xlim, ylim, zlim} the extent of the x/y/z-axes \item
\code{xlab, ylab, zlab} the names of the covariates for the
x/y/z-axes \item \code{dim} the dimensionality of the
effect \item \code{main} the label of the smooth term (a
short label, same as the one used in \code{summary.pffr})
}}
}
\description{
Returns estimated coefficient functions/surfaces
\eqn{\beta(t), \beta(s,t)} and estimated smooth effects
\eqn{f(z), f(x,z)} or \eqn{f(x, z, t)} and their point-wise
estimated standard errors. Not implemented for smooths in
more than 3 dimensions.
}
\details{
The \code{seWithMean}-option corresponds to the
\code{"iterms"}-option in \code{\link[mgcv]{predict.gam}}.
The \code{sandwich}-options works as follows: Assuming that
the residual vectors \eqn{\epsilon_i(t), i=1,\dots,n} are
i.i.d. realizations of a mean zero Gaussian process with
covariance \eqn{K(t,t')}, we can construct an estimator for
\eqn{K(t,t')} from the \eqn{n} replicates of the observed
residual vectors. The covariance matrix of the stacked
observations vec\eqn{(Y_i(t))} is then given by a
block-diagonal matrix with \eqn{n} copies of the estimated
\eqn{K(t,t')} on the diagonal. This block-diagonal matrix
is used to construct the "meat" of a sandwich covariance
estimator, similar to Chen et al. (2012), see reference
below.
}
\author{
Fabian Scheipl
}
\references{
Chen H., Wang Y., Paik C.M., Choi A. (2012). A marginal
approach to reduced-rank penalized spline smoothing for
multilevel data. \emph{Journal of the American Statistical
Association}, under revision.
\url{http://www.columbia.edu/~yw2016/Marginal Spline6.pdf}
}
\seealso{
\code{\link[mgcv]{plot.gam}},
\code{\link[mgcv]{predict.gam}} which this routine is based
on.
}

