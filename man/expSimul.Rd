\name{expSimul}
\alias{expSimul}
\title{Simulates from the exponential distribution}
\usage{
  expSimul(param, nsim, extraArgs, nObs, ...)
}
\arguments{
  \item{nObs}{Length of each simulated time serie.}

  \item{nsim}{Number of simulations from the model.}

  \item{param}{A vector of length 1 or a matrix of size
  nsim by 1.  If param is a vector each of the nsim
  simulations will use the same parameters, if it's a
  matrix the i-th simulation will use the i-th row of param
  as parameters.}

  \item{extraArgs}{Only needed for compatibility with the
  "synlik" package.}
}
\value{
  A matrix of size nsim by nObs where each row is a sample
  from the distribution.
}
\description{
  Simulates from the exponential distribution
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

