\name{PEER.Sim, Q}
\alias{PEER.Sim}
\alias{Q}
\docType{data}
\title{
Simulated longitudinal data with functional predictor and scalar response, and structural information associated with predictor function }

\description{
\code{PEER.Sim} contains simulated observations from 100 subjects, each observed at 4 distinct timepoints. At each timepoint bumpy predictor profile is generated randomly and the scalar response variable is generated considering a time-varying regression function and subject intercept. Accompanying the functional predictor and scalar response are the subject ID numbers and time of measurements. 

\code{Q} represents the 7 x 100 matrix where each row provides structural information about the functional predictor profile for data \code{PEER.Sim}. For specific details about the simulation and Q matrix, please refer to Kundu et. al. (2012).
}

\usage{
data(PEER.Sim)
data(Q)
}

\format{\code{PEER.Sim} is a data frame with 400 observations (100 subjects,
    4 time points per subject) and 4 variables:
\describe{
  \item{id}{Subject ID number}
  \item{t}{Time of measurement}
  \item{Y}{Scalar response}
  \item{W}{A \code{400 x 100} matrix containing the functional predictor profile}
}
\code{Q} is a \code{7 x 100} matrix containing structural information about
  \code{PEER.Sim$W}.
}

\references{
Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012). Longitudinal functional models with structured penalties. (please contact J. Harezlak at harezlak@iupui.edu)
}