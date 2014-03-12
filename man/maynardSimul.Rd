\name{maynardSimul}
\alias{maynardSimul}
\title{Simulates from the Maynard-Smith map}
\usage{
  maynardSimul(param, nsim, extraArgs, randInit = TRUE,
    initVal = 1, ...)
}
\arguments{
  \item{nObs}{Length of each simulated time serie.}

  \item{nsim}{Number of simulations from the model.}

  \item{A}{vector of length 4 or a matrix of size nsim by
  4.  If param is a vector each of the nsim simulations
  will use the same parameters, if it's a matrix the i-th
  simulation will use the i-th row of param as parameters.}

  \item{extraArgs}{Only needed for compatibility with the
  "synlik" package.}

  \item{nBurn}{Number of initial steps to be discarded
  before saving the following nObs steps.}

  \item{randInit}{If TRUE the initial values of each of the
  nsim run will be simulated using runif(nsim, 0, 1), if
  FALSE each of the nsim simulations will start from
  initVal.}

  \item{initVal}{If randInit == FALSE each of the nsim
  simulations will start from this value.}

  \item{...}{Only needed for compatibility with the
  "synlik" package.}
}
\value{
  A matrix of size nsim by nObs where each row is a
  trajectory simulated from the model.
}
\description{
  Simulates from the Maynard-Smith map
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

