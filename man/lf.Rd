% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lf.R
\name{lf}
\alias{lf}
\title{Construct an FLM regression term}
\usage{
lf(
  X,
  argvals = NULL,
  xind = NULL,
  integration = c("simpson", "trapezoidal", "riemann"),
  L = NULL,
  presmooth = NULL,
  presmooth.opts = NULL,
  ...
)
}
\arguments{
\item{X}{functional predictors, typically expressed as an \code{N} by \code{J} matrix,
where \code{N} is the number of columns and \code{J} is the number of
evaluation points. May include missing/sparse functions, which are
indicated by \code{NA} values. Alternatively, can be an object of class
\code{"fd"}; see \code{[fda]{fd}}.}

\item{argvals}{indices of evaluation of \code{X}, i.e. \eqn{(t_{i1},.,t_{iJ})} for
subject \eqn{i}. May be entered as either a length-\code{J} vector, or as
an \code{N} by \code{J} matrix. Indices may be unequally spaced. Entering
as a matrix allows for different observations times for each subject. If
\code{NULL}, defaults to an equally-spaced grid between 0 or 1 (or within
\code{X$basis$rangeval} if \code{X} is a \code{fd} object.)}

\item{xind}{same as argvals. It will not be supported in the next version of refund.}

\item{integration}{method used for numerical integration. Defaults to \code{"simpson"}'s rule
for calculating entries in \code{L}. Alternatively and for non-equidistant grids,
\code{"trapezoidal"} or \code{"riemann"}.}

\item{L}{an optional \code{N} by \code{ncol(argvals)} matrix giving the weights for the numerical
integration over \code{t}. If present, overrides \code{integration}.}

\item{presmooth}{string indicating the method to be used for preprocessing functional predictor prior
to fitting. Options are \code{fpca.sc}, \code{fpca.face}, \code{fpca.ssvd}, \code{fpca.bspline}, and
\code{fpca.interpolate}. Defaults to \code{NULL} indicating no preprocessing. See
\code{{create.prep.func}}.}

\item{presmooth.opts}{list including options passed to preprocessing method
\code{{create.prep.func}}.}

\item{...}{optional arguments for basis and penalization to be passed to
\code{mgcv::s}. These could include, for example,
\code{"bs"}, \code{"k"}, \code{"m"}, etc. See \code{[mgcv]{s}} for details.}
}
\value{
a list with the following entries
  \item{\code{call}}{a \code{call} to \code{te} (or \code{s}, \code{t2}) using the appropriately
    constructed covariate and weight matrices}
  \item{\code{argvals}}{the \code{argvals} argument supplied to \code{lf}}
  \item{\code{L}}{the  matrix of weights used for the integration}
  \item{\code{xindname}}{the name used for the functional predictor variable in the \code{formula}
    used by \code{mgcv}}
  \item{\code{tindname}}{the name used for \code{argvals} variable in the \code{formula} used by \code{mgcv}}
  \item{\code{LXname}}{the name used for the \code{L} variable in the \code{formula} used by \code{mgcv}}
  \item{\code{presmooth}}{the \code{presmooth} argument supplied to \code{lf}}
  \item{\code{prep.func}}{a function that preprocesses data based on the preprocessing method specified in \code{presmooth}. See
    \code{{create.prep.func}}}
}
\description{
Defines a term \eqn{\int_{T}\beta(t)X_i(t)dt} for inclusion in an \code{mgcv::gam}-formula (or
\code{{bam}} or \code{{gamm}} or \code{gamm4:::gamm}) as constructed by
\code{{pfr}}, where \eqn{\beta(t)} is an unknown coefficient
function and \eqn{X_i(t)} is a functional predictor on the closed interval
\eqn{T}. See
\code{{smooth.terms}} for a list of basis and penalty options; the
default is thin-plate regression splines, as this is the default option
for \code{[mgcv]{s}}.
}
\examples{
data(DTI)
DTI1 <- DTI[DTI$visit==1 & complete.cases(DTI),]

# We can apply various preprocessing options to the DTI data
fit1 <- pfr(pasat ~ lf(cca, k=30), data=DTI1)
fit2 <- pfr(pasat ~ lf(cca, k=30, presmooth="fpca.sc",
                       presmooth.opts=list(nbasis=8, pve=.975)), data=DTI1)
fit3 <- pfr(pasat ~ lf(cca, k=30, presmooth="fpca.face",
                       presmooth.opts=list(m=3, npc=9)), data=DTI1)
fit4 <- pfr(pasat ~ lf(cca, k=30, presmooth="fpca.ssvd"), data=DTI1)
fit5 <- pfr(pasat ~ lf(cca, k=30, presmooth="bspline",
                       presmooth.opts=list(nbasis=8)), data=DTI1)
fit6 <- pfr(pasat ~ lf(cca, k=30, presmooth="interpolate"), data=DTI1)

# All models should result in similar fits
fits <- as.data.frame(lapply(1:6, function(i)
  get(paste0("fit",i))$fitted.values))
names(fits) <- c("none", "fpca.sc", "fpca.face", "fpca.ssvd", "bspline", "interpolate")
pairs(fits)

}
\references{
Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011).
Penalized functional regression. \emph{Journal of Computational and Graphical
Statistics}, 20(4), 830-851.

Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D. (2012). Longitudinal
penalized functional regression for cognitive outcomes on neuronal tract
measurements. \emph{Journal of the Royal Statistical Society: Series C},
61(3), 453-469.
}
\seealso{
\code{{pfr}}, \code{{af}}, mgcv's \code{{smooth.terms}}
 and \code{{linear.functional.terms}}; \code{{pfr}} for additional examples
}
\author{
Mathew W. McLean \email{mathew.w.mclean@gmail.com}, Fabian Scheipl,
  and Jonathan Gellar
}
