\name{plot.fpcr}
\Rdversion{1.1}
\alias{plot.fpcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default plotting for functional principal component regression output}
\description{Inputs an object created by \code{\link{fpcr}}, and plots the estimated coefficient function.}
\usage{
\method{plot}{fpcr}(x, se=TRUE, col=1, lty=c(1,2,2), xlab="", 
                    ylab="Coefficient function", \dots) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{"\link{fpcr}"}.}
  \item{se}{
if \code{TRUE} (the default), upper and lower lines are added at 2 standard errors (in the Bayesian sense; see Wood, 2006) above and below the coefficient function estimate.  If a positive number is supplied, the standard error is instead multiplied by this number.}
  \item{col}{
color for the line(s).  This should be either a number, or a vector of length 3 for the coefficient function estimate, lower bound, and upper bound, respectively.}
  \item{lty}{
line type(s) for the coefficient function estimate, lower bound, and upper bound.}
\item{xlab, ylab}{x- and y-axis labels.}
\item{\dots}{
other arguments passed to the underlying plotting function.}
}

\value{
None; only a plot is produced.}
\references{
Wood, S. N. (2006). \emph{Generalized Additive Models: An Introduction with R}. Boca Raton, FL: Chapman & Hall.}

\author{
Philip Reiss \email{phil.reiss@nyumc.org}}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{fpcr}}, which includes an example.
}


