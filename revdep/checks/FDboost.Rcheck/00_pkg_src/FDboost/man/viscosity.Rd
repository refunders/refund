\name{viscosity}
\alias{viscosity}
\docType{data}
\title{ Viscosity of resin over time}
\description{

  In an experimental setting the viscosity of resin was measured over time
  to asses the curing process depending on 5 binary factors (low-high).

}
\usage{data("viscosity")}
\format{
  A data list with 64 observations on the following 7 variables.
  \describe{
    \item{\code{visAll}}{viscosity measures over all available time points}
    \item{\code{timeAll}}{time points of viscosity measures}
    \item{\code{T_C}}{ temperature of tools}
    \item{\code{T_A}}{temperature of resin}
    \item{\code{T_B}}{temperature of curing agent}
    \item{\code{rspeed}}{rotational speed}
    \item{\code{mflow}}{mass flow}
  }
}
\details{
The aim is to determine factors that affect the curing process in the mold.
The desired viscosity-curve has low values in the beginning followed
by a sharp increase.
Due to technical reasons the measuring method of the rheometer has to be
changed in a certain range of viscosity. The first observations are measured
by rotation of a blade giving observations every two seconds,
the later observations are measured through oscillation of a blade giving
observations every ten seconds. In the later observations the resin is quite
hard so the measurements should be interpreted as a qualitative measure of hardening.
}
\source{
   Wolfgang Raffelt, Technical University of Munich, Institute for Carbon Composites
}
\examples{

 data("viscosity", package = "FDboost")
 ## set time-interval that should be modeled
 interval <- "101"

 ## model time until "interval" and take log() of viscosity
 end <- which(viscosity$timeAll==as.numeric(interval))
 viscosity$vis <- log(viscosity$visAll[,1:end])
 viscosity$time <- viscosity$timeAll[1:end]
 # with(viscosity, funplot(time, vis, pch=16, cex=0.2))

 ## fit median regression model with 100 boosting iterations,
 ## step-length 0.4 and smooth time-specific offset
 ## the factors are in effect coding -1, 1 for the levels
 mod <- FDboost(vis ~ 1 + bols(T_C, contrasts.arg = "contr.sum", intercept=FALSE)
                + bols(T_A, contrasts.arg = "contr.sum", intercept=FALSE),
                timeformula=~bbs(time, lambda=100),
                numInt="equal", family=QuantReg(),
                offset=NULL, offset_control = o_control(k_min = 9),
                data=viscosity, control=boost_control(mstop = 100, nu = 0.4))
 summary(mod)

}
\keyword{datasets}
