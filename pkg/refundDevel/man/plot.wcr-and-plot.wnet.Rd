\name{plot.wcr/wnet}
\alias{plot.wnet}
\alias{plot.wcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default plotting for wavelet-domain scalar-on-function regression
}
\description{
Plots the coefficient function/image estimates produced by \code{\link{wcr}} and \code{\link{wnet}}.
}
\usage{
\method{plot}{wcr}(x, xlabel = "", ylabel = "Coefficient function", which.dim = 1, slices = NULL, 
set.mfrow = TRUE, image.axes = FALSE, ...)
\method{plot}{wnet}(x, xlabel = "", ylabel = "Coefficient function", which.dim = 1, slices = NULL, 
set.mfrow = TRUE, image.axes = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{"\link{wcr}"} or \code{"\link{wnet}"}.
}
  \item{xlabel, ylabel}{
for 1D functional predictors, x- and y-axis labels.
}
  \item{which.dim, slices}{
for 3D image predictors, the dimension (1, 2 or 3) and slices to use for plotting; see Details. 
}
  \item{set.mfrow}{
logical value: for 3D predictors, if \code{TRUE} (the default), the function will try to set an appropriate value of the \code{mfrow} plotting parameter to display the number of slices specified.  Otherwise you may wish to set \code{mfrow} outside the function call.
}
\item{image.axes}{for 2D and 3D predictors, the \code{axes} argument passed to \code{\link{image}}.}
  \item{\dots}{
additional parameters passed to \code{\link{plot}} or \code{\link{image}}.
}
}
\details{
As an example of how \code{which.dim} and \code{slices} are used, suppose we set \code{which.dim=2} and \code{slices=7:9}. Then three 2D slices of the coefficient image estimate \code{x$fhat} are displayed: \code{x$fhat[ , 7, ]}, \code{x$fhat[ , 8, ]}, and \code{x$fhat[ , 9, ]}. 
}

\author{
Lan Huo \email{lan.huo@nyumc.org}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wcr}} and \code{\link{wnet}}; the latter includes examples.
}
