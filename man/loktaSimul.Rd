\name{loktaSimul}
\alias{loktaSimul}
\title{Simulates from the Lokta-volterra model}
\usage{
  loktaSimul(param, nsim, extraArgs, nBurn, totTime = 15,
    nSteps = 5, randInit = TRUE, x0 = 1, y0 = 0.5,
    onlyPrey = TRUE, ...)
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

  \item{totTime}{Total period of time to which the
  simulations correspond (ex: 15 years).}

  \item{nSteps}{Number of finite differences steps between
  2 observations.}

  \item{randInit}{random initialization for prey and
  predator?}

  \item{x0}{Initial state for the prey.}

  \item{y0}{Initial state for the predator.}

  \item{onlyPrey}{If TRUE only the prey simulations are
  returned.}

  \item{...}{Only needed for compatibility with the
  "synlik" package.}
}
\value{
  If onlyPrey == FALSE it returns a list where [["prey"]]
  and [["predator"]] are two matrices of size nsim by nObs,
  if onlyPrey == TRUE it return a matrix containing only
  the simulations for the prey.
}
\description{
  Simulates from the Lokta-volterra model
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

