\name{DTI2}
\alias{DTI2}
\docType{data}
\title{
Diffusion Tensor Imaging: tract profiles and outcomes
}
\description{
A DTI dataset used in Swihart et al 2012.
Fractional anisotropy tract profiles for the corpus callosum (cca) and
the right corticospinal tract (rcst). Accompanying the tract profiles
are the subject ID numbers, visit number, and Paced Auditory Serial Addition Test
(pasat) score. We thank Dr. Daniel Reich for making this dataset
available.
}
\usage{data(DTI2)}
\format{
  A data frame made up of
  \describe{
    \item{\code{cca} }{a 340 x 93 matrix of fractional anisotropy tract profiles from the corpus callosum;}
    \item{\code{rcst} }{a 340 x 55 matrix of fractional anisotropy tract profiles from the right corticospinal tract;}
    \item{\code{id} }{numeric vector of subject ID numbers;}
    \item{\code{visit} }{numeric vector of the subject-specific visit numbers;}
    \item{\code{pasat }}{numeric vector containing the PASAT score at each visit.}
  }
}
\references{
Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011). Penalized functional regression. \emph{Journal of Computational and Graphical Statistics}, 20(4), 830-851.

Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich,
D. (2012). Longitudinal penalized functional regression for cognitive
outcomes on neuronal tract measurements. \emph{Journal of the Royal
  Statistical Society: Series C}, 61(3), 453--469.

Swihart, Bruce J., Goldsmith, Jeff; and Crainiceanu, Ciprian M. (July 2012).
Testing For Functional Effects 
\emph{Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper 247.}
http://biostats.bepress.com/jhubiostat/paper247

}
\keyword{datasets}

