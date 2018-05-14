\name{predict.face.sparse}
\alias{predict.face.sparse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subject-specific curve prediction from a face.sparse fit
}
\description{
Predict subject-specific curves based on a fit from "face.sparse".
}
\usage{
\method{predict}{face.sparse}(object, newdata,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted object from the R function  "face.sparse".
}
\item{newdata}{
a data frame with three arguments: 
(1) \code{argvals}: observation times;
(2) \code{subj}: subject indices;
(3) \code{y}: values of observations.
NA values are allowed in "y" but not in the other two.
}
\item{\dots}{further arguments passed to or from other methods.}
}
\details{
This function makes prediction based on observed data for each subject. So for each subject,
it requires at least one observed data. For the time points prediction is desired but
no observation is available, just make the corresponding \code{data$y} as NA.
}
\value{
\item{object}{A "face.sparse" fit}
\item{newdata}{Input}
\item{y.pred,mu.pred,Chat.pred, Chat.diag.pred, var.error.pred}{Predicted/estimated objects at
the observation time points in \code{newdata}}
\item{scores}{if \code{calculate.scores} in \code{object} is TRUE (typically FALSE),
then predicted scores will be calculated.}
\item{\dots}{\dots}
}

\references{
Luo Xiao, Cai Li, William Checkley and Ciprian Crainiceanu, Fast covariance estimation
for sparse functional data, Stat. Comput., \doi{10.1007/s11222-017-9744-8}.
}
\author{
Luo Xiao <lxiao5@ncsu.edu>
}


\examples{

#See the examples for "face.sparse".
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~face.sparse}
\keyword{~prediction}% __ONLY ONE__ keyword per line
