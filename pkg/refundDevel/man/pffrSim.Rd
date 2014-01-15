\name{pffrSim}
\alias{pffrSim}
\title{Simulate example data for pffr}
\usage{
  pffrSim(scenario = "all", n = 100, nxgrid = 40,
    nygrid = 60, SNR = 10, propmissing = 0)
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

  \item{propmissing}{proportion of missing data in the
  response, default = 0. See Details.}
}
\value{
  a named list with the simulated data, and the true
  components of the predictor etc as attributes.
}
\description{
  Simulates example data for \code{\link{pffr}}. Scenario
  "all" generates data from a complex multivariate model
  \eqn{Y_i(t) = \mu(t) + \int X_{1i}(s)\beta_1(s,t)ds +
  \int X_{2i}(s)\beta_2(s,t)ds + \text{xlin} \beta_3(t) +
  f(\text{xte1}, \text{xte2}) + f(\text{xsmoo}, t) +
  \beta_4 \text{xconst} + \epsilon_i(t)}. Scenarios "ff",
  "te", "smoo", "lin", "const" generate data from simpler
  models containing only the respective term(s) in the
  model equation given above. Specifiying a vector-valued
  scenario will generate data from a combination of the
  respective terms. See source code for details.\cr
  Sparse/irregular response trajectories can be generated
  by setting \code{propmissing} to something greater than 0
  (and smaller than 1). The return object then also
  includes a \code{ydata}-item with the sparsified data.
}

