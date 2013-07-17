\name{predict.pfr}
\alias{predict.pfr}
\title{
Prediction for penalized functional regression
}

\description{
Given a \code{\link{pfr}} object and new data, produces fitted values on both a
population and subject-specific scale.
}

\usage{
  \method{predict}{pfr}(object, new.data=NULL, levels=NULL, ...)
}

\arguments{
  \item{object}{an object returned by \code{\link{pfr}}.}
  \item{new.data}{the covariate and functional predictor values for which predicted values are desired.}
  \item{levels}{currently not supported; both population- and subject-level fitted values are returned.}
  \item{...}{additional arguments.}
}

\details{
Predicting for new data in a penalized functional regression setting
takes care.  Utilizing \code{\link{fpca.sc}} in both \code{\link{pfr}} and
\code{predict.pfr} allows for the correct estimation of new scores
using the original fit's basis functions.  In the spirit of
\code{\link[nlme]{predict.lme}}, we can estimate for population and (if the subjects
are in both the original and new data) subject-specific levels.
}

\value{A list with components
   \item{fitted.vals.level.0 }{predicted outcomes at population level}
  \item{fitted.vals.level.1 }{predicted outcomes at subject-specific
    level (if applicable)}
}

\references{
Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011).
Penalized functional regression. \emph{Journal of Computational and Graphical Statistics}, 20(4), 830-851.

Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D. (2012).
Longitudinal penalized functional regression for cognitive outcomes on
neuronal tract measurements. \emph{Journal of the Royal Statistical Society: Series C}, 61(3), 453--469.

Swihart, B. J., Goldsmith, J., and Crainiceanu, C. M. (2012).
Testing for functional effects. Johns Hopkins University Dept. of Biostatistics Working Paper 247, available at \url{http://biostats.bepress.com/jhubiostat/paper247}
}

\seealso{\code{\link{pfr}}, \code{\link{rlrt.pfr}}, \code{\link[nlme]{predict.lme}}}
\author{
  Jeff Goldsmith <jeff.goldsmith@columbia.edu>
  and Bruce Swihart <bswihart@jhsph.edu>
}

\examples{
\dontrun{
##################################################################
#########               DTI Data Example                 #########
##################################################################


## load and reassign the data;
data(DTI2)
Y  <- DTI2$pasat ## PASAT outcome
id <- DTI2$id    ## subject id
W1 <- DTI2$cca   ## Corpus Callosum
W2 <- DTI2$rcst  ## Right corticospinal
V  <- DTI2$visit ## visit

## prep scalar covariate
visit.1.rest <- matrix(as.numeric(V > 1), ncol=1)
covar.in <- visit.1.rest 

## kz -- put it too high and the computation time bogs down.
pfr.obj <- pfr(Y = Y, covariates=covar.in, funcs = list(W1), subj = id, kz = 10, kb = 50)
## first column is not using random effects; second column is
head(cbind(pfr.obj$fitted.vals.level.0, pfr.obj$fitted.vals.level.1))

## do some predictions with predict.pfr()
## run predict.pfr() on all rows of data -- hope for same results as pfr fit above.
same.inputs <- predict.pfr(pfr.obj, new.data = list(subj=id, covariates=covar.in, funcs=list(W1)))
head(cbind(same.inputs[[1]], same.inputs[[2]]))
## run predict.pfr() on first 300 rows of data
subset.inputs <- predict.pfr(pfr.obj,
                             new.data = list(subj=head(id,300),
                                             covariates=head(covar.in,300),
                                             funcs=list(W1[1:300,])))
head(cbind(subset.inputs[[1]], subset.inputs[[2]]))
## compare the first 300 vs. the original full data
plot(subset.inputs[[1]], same.inputs[[1]][1:300]); abline(a=0,b=1,col="blue")
summary(subset.inputs[[1]] - same.inputs[[1]][1:300])

## try one where we have different subjects for predictions
## ids 7:12 are new
## and we just double-down the covariates and funcs from the first 6 to the second.
## check level.1 and level.0 predictions
subset.inputs <- predict.pfr(pfr.obj,
                             new.data = list(subj=c(head(id,6), 7:12),
                                             covariates=rbind(head(covar.in,6),head(covar.in,6)),
                                             funcs=list(rbind(W1[1:6,],W1[1:6,]))))
## first 6 have both levels; the last 6 do not have subject-specific
head(cbind(subset.inputs[[1]], subset.inputs[[2]]),12)

## see that whether predicted on subset or full data, the predicted
## values are the same for individuals regardless of who is in
## prediction set
## compare the first 6 vs. the original full data: 6 on the 45 degree line and 6 off.
plot(subset.inputs[[1]], same.inputs[[1]][1:12]); abline(a=0,b=1,col="blue")
summary(subset.inputs[[1]] - same.inputs[[1]][1:12])
}
}