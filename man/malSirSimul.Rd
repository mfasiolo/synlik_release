\name{malSirSimul}
\alias{malSirSimul}
\title{Simulator for the S2EI2 malaria model of Bhadhra et at. 2011.}
\usage{
  malSirSimul(procparam, obsparam, splineparam, initStates,
    population, rainfall, nSimul, nObs, nSteps = 30,
    obsInterval = 1/12, T0 = 0,
    splineSettings = list(nBasis = as.integer(6), degree = as.integer(2), period = 1))
}
\arguments{
  \item{procparam:}{matrix of size (nSimul by 10), the i-th
  row is the set of process parameters used for the i-th
  simulation. The columns represent, in order, the _LOG_ of
  parameters mu_EI1, mu_I1I2, mu_I2S2, mu_S2S1, mu_I1S1, q,
  covarCoeff, sigma, c and tau.}

  \item{obsparam:}{matrix of size (nSimul by 2), the i-th
  row is the set of parameters of the observational process
  used for the i-th simulation. The columns represent, in
  order, the _LOG_ of parameters rho and psiSquare.}

  \item{splineparam:}{matrix of size (nSimul by
  splineSettings[["nBasis"]]), the i-th row is the set of
  spline coefficients used for the i-th simulation.}

  \item{initStates:}{matrix of size (nSimul by 7), the i-th
  row is the set of initial values used for the i-th
  simulation. The columns represent, in order, the states
  S1, S2, E, I1, I2, lambda1 and lambda2.  NB: S1, S2, E,
  I1 and I2 have to be integers that sum up to
  population_[0]}

  \item{population:}{vector of length (nObs) representing
  the total population at each month. In this model I'm
  assuming that the population stays constant between
  months (which is roughly true).}

  \item{rainfall:}{vector of length (nObs * nSteps)
  representing the rainfall at each step.}

  \item{nSteps:}{finite difference step in which one month
  is divided (mininum nStep = 1, i.e. you directly step
  from one month to the next).}

  \item{nSimul:}{number of simulated paths}

  \item{obsInterval:}{the length of the interval between 2
  observation (which will be subdivided in nSteps).
  Typically it is going to be 1/12 so each observation
  represents a month.}

  \item{nObs:}{number of observations (months) of each
  simulated trajectory.}

  \item{T0:}{initial time of the simulations. It is used to
  set up the spline basis.}

  \item{splineSettings:}{named list where [["nBasis"]],
  [["degree"]] and [["period"]] are the (self-explanatory)
  settings used to set up the spline basis for the seasonal
  component of the model.}
}
\value{
  A named list where S1, S2, E, I1, I2, lambda1 and lambda2
  are matrices of size (nSimul by nObs) each row
  representing a simulated trajectory of a particular
  state.
}
\description{
  Simulator for the S2EI2 malaria model of Bhadhra et at.
  2011.
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

