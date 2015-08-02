##' Simulated longitudinal data with functional predictor and scalar response,
##' and structural information associated with predictor function
##'
##' \code{PEER.Sim} contains simulated observations from 100 subjects, each
##' observed at 4 distinct timepoints. At each timepoint bumpy predictor
##' profile is generated randomly and the scalar response variable is generated
##' considering a time-varying regression function and subject intercept.
##' Accompanying the functional predictor and scalar response are the subject
##' ID numbers and time of measurements.
##'
##' \code{Q} represents the 7 x 100 matrix where each row provides structural
##' information about the functional predictor profile for data
##' \code{PEER.Sim}. For specific details about the simulation and Q matrix,
##' please refer to Kundu et. al. (2012).
##'
##'
##' @name PEER.Sim
##' @aliases PEER.Sim Q
##' @docType data
##' @format The data frame \code{PEER.Sim} is made up of subject ID
##' number(\code{id}), subject-specific time of measurement (\code{t}),
##' functional predictor profile (\code{W.1-W.100}) and scalar response
##' (\code{Y})
##' @references Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012).
##' Longitudinal functional models with structured penalties. (please contact
##' J. Harezlak at \email{harezlak@@iupui.edu})
NULL
