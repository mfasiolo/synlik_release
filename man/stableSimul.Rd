\name{stableSimul}
\alias{stableSimul}
\title{Simulates from the alpha stable family of distribution}
\usage{
  stableSimul(param, nsim, extraArgs, ...)
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

  \item{...}{Only needed for compatibility with the
  "synlik" package.}
}
\value{
  A matrix of size nsim by nObs where each row is a sample
  from the distribution.
}
\description{
  Simulates from the alpha stable family of distribution
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

