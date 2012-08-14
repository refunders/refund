\name{longitudinal PEER simulation}
\alias{PEER.Sim1}
\alias{Q}
\docType{data}
\title{
Simulated longitudinal data with functional predictor and scalar response, and structural information associated with predictor function }

\description{The data consist of simulated observations given by \code{PEER.Sim1} and structural information contained in matrix \code{Q}.
 For specific details, please refer to Kundu et. al. (2012).
}

\usage{
data(PEER.Sim1)
data(Q)
}

\format{\code{PEER.Sim1} is a data frame of simulated observations from 100 subjects, each observed at 4 distinct time points. At each time point a bumpy predictor profile is generated randomly, and the scalar response variable is generated using a time-varying regression function and subject intercept. 
The data frame is made up of 
\describe{
\item{\code{id} }{subject ID number;}
\item{\code{t} }{subject-specific time of measurement;}
\item{\code{W.1,\dots,W.100} }{functional predictor profile; and}
\item{\code{Y} }{scalar response.}}

\code{Q} is a 9 x 100 matrix where each row provides structural information about the functional predictor profile for data \code{PEER.Sim1}.}

\seealso{\code{\link{lpeer}}}

\references{
Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012). Longitudinal functional models with structured penalties. (please contact J. Harezlak at harezlak@iupui.edu)
}
\keyword{datasets}