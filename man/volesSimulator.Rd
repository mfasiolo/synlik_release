\name{volesSimulator}
\alias{volesSimulator}
\title{Simulates from the voles models of Turchin and Ellner (2000).}
\usage{
  volesSimulator(param, nsim, nObs, extraArgs = NULL,
    model = "full", nBurn = 100, nSteps = 50, T0 = 0, ...)
}
\arguments{
  \item{nObs}{Length of each simulated time serie.}

  \item{nsim}{Number of simulations from the model.}

  \item{param}{If model == "full" it must be a vector of
  length 9 or a matrix of size nsim by 9, if model ==
  "standard" it must be a vector of length 7 or a matrix of
  size nsim by 7.  If param is a vector each of the nsim
  simulations will use the same parameters, if it's a
  matrix the i-th simulation will use the i-th row of param
  as parameters.}

  \item{extraArgs}{Only needed for compatibility with the
  "synlik" package.}

  \item{model}{If model == "full" the function will
  simulate from the full model described in eq. 2 of
  Turchin and Ellner (2000), model == "standard" the
  function will simulate from the full model described in
  eq. 6 of Turchin and Ellner (2000).}

  \item{nBurn}{Number of initial steps to be discarded
  before saving the following nObs steps.}

  \item{nSteps}{Number of finite differences steps between
  2 observations (an interval which corresponds to a
  semester).}

  \item{...}{Only needed for compatibility with the
  "synlik" package.}
}
\value{
  A list where [[""voles"]] and [["weasels"]] are two
  matrices of size nsim by nObs where each row is a
  trajectory simulated from the model.
}
\description{
  Simulates from the voles models of Turchin and Ellner
  (2000).
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}
\references{
  Turchin and Ellner (2000), LIVING ON THE EDGE OF CHAOS:
  POPULATION DYNAMICS OF FENNOSCANDIAN VOLES.
}

