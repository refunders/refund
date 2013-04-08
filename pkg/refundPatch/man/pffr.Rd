\name{pffr}
\alias{pffr}
\title{Penalized function-on-function regression}
\usage{
  pffr(formula, yind, data = NULL, algorithm = NA,
    method = "REML", tensortype = c("te", "t2"),
    bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 20, m = c(2, 1)), ...)
}
\arguments{
  \item{formula}{a formula with special terms as for
  \code{\link[mgcv]{gam}}, with additional special terms
  \code{\link{ff}()} and \code{c()}}

  \item{yind}{a vector with length equal to the number of
  columns of the matrix of functional responses giving the
  vector of evaluation points \eqn{(t_1, \dots ,t_{G})}.
  If \code{formula} contains an \code{\link{ff}}-term which
  specifies \code{yind} this is used. If neither is given,
  \code{yind} is \code{1:ncol(<response>)}.}

  \item{algorithm}{the name of the function used to
  estimate the model. Defaults to \code{\link[mgcv]{gam}}
  if the matrix of functional responses has less than
  \code{2e5} data points and to \code{\link[mgcv]{bam}} if
  not. "gamm" (see \code{\link[mgcv]{gamm}}) and "gamm4"
  (see \code{\link[gamm4]{gamm4}}) are valid options as
  well.}

  \item{data}{an (optional) \code{data.frame} or a named
  list containing the data. The functional response and
  functional covariates have to be supplied as n by <no. of
  evaluations> matrices, i.e. each row is one functional
  observation. The model is then fitted with the data in
  long format, i.e., the rows of the matrix of the
  functional response evaluations \eqn{Y_i(t)} are stacked
  into one long vector and the covariates are
  expanded/repeated correspondingly. This means the models
  get quite big fairly fast, since the effective number of
  rows in the design matrix is number of observations times
  number of evaluations of \eqn{Y(t)} per observation.\cr}

  \item{method}{Defaults to \code{"REML"}-estimation,
  including of unknown scale. See \code{\link[mgcv]{gam}}
  for details.}

  \item{bs.yindex}{a named (!) list giving the parameters
  for spline bases on the index of the functional response.
  Defaults to \code{list(bs="ps", k=5, m=c(2, 1))}, i.e. 5
  cubic B-splines bases with first order difference
  penalty.}

  \item{bs.int}{a named (!) list giving the parameters for
  the spline basis for the global functional intercept.
  Defaults to \code{list(bs="ps", k=20, m=c(2, 1))}, i.e.
  20 cubic B-splines bases with first order difference
  penalty.}

  \item{tensortype}{which type of tenor product splines to
  use. One of "\code{\link[mgcv]{te}}" or
  "\code{\link[mgcv]{t2}}", defaults to \code{te}}

  \item{...}{additional arguments that are valid for
  \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.
  \code{weights, subset, offset} are not yet implemented!}
}
\value{
  a fitted \code{pffr}-object, which is a
  \code{\link[mgcv]{gam}}-object with some additional
  information in an \code{pffr}-entry. If \code{algorithm}
  is \code{"gamm"} or \code{"gamm4"}, only the \code{$gam}
  part of the returned list is modified in this way.
}
\description{
  Implements additive regression for functional and scalar
  covariates and functional responses. This function is a
  wrapper for \code{mgcv}'s \code{\link[mgcv]{gam}} and its
  siblings to fit models of the general form \cr
  \eqn{E(Y_i(t)) = g(\mu(t) + \int X_i(s)\beta(s,t)ds +
  f(z_{1i}, t) + f(z_{2i}) + z_{3i} \beta_3(t) + \dots }\cr
  with a functional (but not necessarily continuous)
  response \eqn{Y(t)}, response function \eqn{g},
  (optional) smooth intercept \eqn{\mu(t)}, (multiple)
  functional covariates \eqn{X(t)} and scalar covariates
  \eqn{z_1}, \eqn{z_2}, etc.
}
\section{Details}{
  The routine can estimate \enumerate{ \item (nonlinear,
  and possibly multivariate) effects of (one or multiple)
  scalar covariates \eqn{z} that vary smoothly over the
  index \eqn{t} of \eqn{Y(t)} (e.g. \eqn{f(z_{1i}, t)},
  specified in the \code{formula} simply as \code{~s(z1)}),
  \item (nonlinear) effects of scalar covariates that are
  constant over \eqn{t} (e.g. \eqn{f(z_{2i})}, specified as
  \code{~c(s(z2))}, or \eqn{\beta_2 z_{2i}}, specified as
  \code{~c(z2)}), \item linear functional effects of scalar
  (numeric or factor) covariates that vary smoothly over
  \eqn{t} (e.g. \eqn{z_{3i} \beta_3(t)}, specified as
  \code{~z3}), \item function-on-function regression terms
  (e.g. \eqn{\int X_i(s)\beta(s,t)ds}, specified as
  \code{~ff(X, yindex=t, xindex=s)}, see \code{\link{ff}}).
  } Use the \code{c()}-notation to denote model terms that
  are constant over the index of the functional
  response.\cr

  Internally, univariate smooth terms without a
  \code{c()}-wrapper are expanded into bivariate smooth
  terms in the original covariate and the index of the
  functional response. Bivariate smooth terms (\code{s(),
  te()} or \code{t2()}) without a \code{c()}-wrapper are
  expanded into trivariate smooth terms in the original
  covariates and the index of the functional response.
  Linear terms for scalar covariates or categorical
  covariates are expanded into varying coefficient terms,
  varying smoothly over the index of the functional
  response. For factor variables, a separate smooth
  function with its own smoothing parameter is estimated
  for each level of the factor.\cr Functional random
  intercepts \eqn{B_{0g(i)}(t)} for a grouping variable
  \code{g} can be specified via \code{~s(g, bs="re")}),
  functional random slopes \eqn{u_i B_{1g(i)}(t)} in a
  numeric variable \code{u} via \code{~s(g, u,
  bs="re")}).\cr The marginal spline basis used for the
  index of the the functional response is specified via the
  \emph{global} argument \code{bs.yindex}. If necessary,
  this can be overriden for any specific term by supplying
  a \code{bs.yindex}-argument to that term in the formula,
  e.g. \code{~s(x, bs.yindex=list(bs="tp", k=7))} would
  yield a tensor product spline for which the marginal
  basis for the index of the response are 7 cubic
  thin-plate spline functions overriding the global default
  for the basis and penalty on the index of the response
  given by the \emph{global} \code{bs.yindex}-argument .\cr
  Use \code{~-1 + c(1) + ...} to specify a model with only
  a constant and no functional intercept. \cr

  The functional response and functional covariates have to
  be supplied as n by <no. of evaluations> matrices, i.e.
  each row is one functional observation. The model is then
  fitted with the data in long format, i.e., the rows of
  the matrix of the functional response evaluations
  \eqn{Y_i(t)} are stacked into one long vector and the
  covariates are expanded/repeated correspondingly. This
  means the models get quite big fairly fast, since the
  effective number of rows in the design matrix is number
  of observations times number of evaluations of \eqn{Y(t)}
  per observation.\cr

  Note that pffr overrides the default identifiability
  constraints (\eqn{\sum_{i,t} \hat f(z_i, x_i, t) = 0})
  implemented in \code{mgcv} for tensor product terms whose
  marginal terms include the index of the functional
  response \eqn{t}.  Instead, \eqn{\sum_i \hat f(z_i, x_i,
  t) = 0} for all \eqn{t} is enforced, so that effects
  varying over \eqn{t} can be interpreted as deviations
  from the global functional intercept. I recommend using
  centered scalar covariates for terms like \eqn{z
  \beta(t)} (\code{~z}) and centered functional covariates
  with \eqn{\sum_i X_i(t) = 0} for all \eqn{t} in
  \code{ff}-terms so that the global functional intercept
  can be interpreted as the global mean function.
}
\examples{

\dontrun{
###############################################################################
# univariate model:
# Y(t) = f(t)  + \\int X1(s)\\beta(s,t)ds + eps
data1 <- pffrSim(scenario="1")
t <- attr(data1, "yindex")
s <- attr(data1, "xindex")
m1 <- pffr(Y ~ ff(X1, yind=t, xind=s), data=data1)
summary(m1)
plot(m1, pers=TRUE)

###############################################################################
# multivariate model:
# Y(t) = f0(t)  + \\int X1(s)\\beta1(s,t)ds + \\int X2(s)\\beta2(s,t)ds +
#		xlin \\beta3(t) + f1(xte1, xte2) + f2(xsmoo, t) + beta4 xconst + eps
data2 <- pffrSim(scenario="2", n=200)
t <- attr(data2, "yindex")
s <- attr(data2, "xindex")
m2 <- pffr(Y ~  ff(X1, xind=s, yind=t) + #linear function-on-function
                ff(X2, xind=s, yind=t) + #linear function-on-function
                xlin  +  #varying coefficient term
                c(te(xte1, xte2)) + #bivariate smooth term in xte1 & xte2, const. over Y-index
                s(xsmoo) + #smooth effect of xsmoo varying over Y-index
                c(xconst), # linear effect of xconst constant over Y-index
        data=data2)
summary(m2)
plot(m2, pers=TRUE)
}
}
\author{
  Fabian Scheipl, Sonja Greven
}
\references{
  Ivanescu, A., Staicu, A.-M., Scheipl, F. and Greven, S.
  (2012).  Penalized function-on-function regression.
  (submitted)
  \url{http://biostats.bepress.com/jhubiostat/paper240/}

  Scheipl, F., Staicu, A.-M. and Greven, S. (2012).
  Additive Mixed Models for Correlated Functional Data.
  (submitted) \url{http://arxiv.org/abs/1207.5947}
}
\seealso{
  \code{\link[mgcv]{smooth.terms}} for details of
  \code{mgcv} syntax and available spline bases
}

