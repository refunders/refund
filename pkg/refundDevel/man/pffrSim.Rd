\name{pffrSim}
\alias{pffrSim}
\title{Simulate example data for pffr()}
\usage{
  pffrSim(scenario = 1, n = 100, nxgrid = 40, nygrid = 60,
    SNR = 10)
}
\arguments{
  \item{scenario}{see description}

  \item{n}{number of observations}

  \item{nxgrid}{number of evaluation points of functional
  covariates}

  \item{nygrid}{number of evaluation points of the
  functional response}

  \item{SNR}{the signal-to-noise ratio for the generated
  data: empirical variance of the additive predictor
  divided by variance of the errors.}
}
\value{
  a named list with the simulated data, and the true
  components of the predictor etc as attributes.
}
\description{
  Simulates example data for \code{\link{pffr()}}. Scenario
  1 generates data from a simple univariate model
  \eqn{Y_i(t) = \mu(t) +\int X_i(s)\beta(s,t)ds +
  \epsilon_i(t)}. Scenario 2 generates data from a complex
  multivariate model \eqn{Y_i(t) = \mu(t) + \int
  X_{1i}(s)\beta_1(s,t)ds + \int X_{2i}(s)\beta_2(s,t)ds +
  x_3 \beta_3(t) + f(x_4, x_5) + f(x6, t) + \beta_4 x_7 +
  \epsilon_i(t)}.
}
\details{
  See source code for details.
}

