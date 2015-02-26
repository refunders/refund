\name{plot.fosr}
\Rdversion{1.1}
\alias{plot.fosr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default plotting of function-on-scalar regression objects}
\description{
Plots the coefficient function estimates produced by \code{fosr()}.}
\usage{
\method{plot}{fosr}(x, split = NULL, titles = NULL, xlabel = "", ylabel = "Coefficient function", 
     set.mfrow = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{"\link{fosr}"}.}
  \item{split}{
value, or vector of values, at which to divide the set of coefficient functions into groups, each plotted on a different scale.  E.g., if set to 1, the first function is plotted on one scale, and all others on a different (common) scale.  If \code{NULL}, all functions are plotted on the same scale.}
  \item{titles}{
character vector of titles for the plots produced, e.g., names of the corresponding scalar predictors.}
  \item{xlabel}{
label for the x-axes of the plots.}
  \item{ylabel}{
label for the y-axes of the plots.}
  \item{set.mfrow}{
logical value: if \code{TRUE}, the function will try to set an appropriate value of the \code{mfrow} parameter for the plots.  Otherwise you may wish to set \code{mfrow} outside the function call.}
  \item{\dots}{
graphical parameters (see \code{\link{par}}) for the plot.}
}



\author{
Philip Reiss \email{phil.reiss@nyumc.org}}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fosr}}, which includes examples.
}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
