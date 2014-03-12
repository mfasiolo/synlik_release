\name{blowSimul}
\alias{blowSimul}
\title{Simulates from the blowfly model}
\usage{
  blowSimul(param, nsim, extraArgs, ...)
}
\arguments{
  \item{nObs}{Length of each simulated time serie.}

  \item{nsim}{Number of simulations from the model.}

  \item{A}{vector of length 2 or a matrix of size nsim by
  2.  If param is a vector each of the nsim simulations
  will use the same parameters, if it's a matrix the i-th
  simulation will use the i-th row of param as parameters.}

  \item{extraArgs}{Only needed for compatibility with the
  "synlik" package.}

  \item{nBurn}{Number of initial steps to be discarded
  before saving the following nObs steps.}
}
\description{
  Simulates from the blowfly model
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com> and Simon Wood.
}

