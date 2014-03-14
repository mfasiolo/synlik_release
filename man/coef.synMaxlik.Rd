\name{coef.synMaxlik}
\alias{coef.synMaxlik}
\title{Get the final parameters' estimates from a synMaxLik object}
\usage{
  \method{coef.synMaxlik}{synMaxlik #} (object, Lag = 10,
    tol = 1e-12, ...)
}
\arguments{
  \item{object}{An object of class "synMaxlik".}

  \item{Lag}{Number of values over which the Hessian is
  averaged to give it's final estimate. Should be a number
  between 1 and object@nIter.}

  \item{tol}{Tolerance used to determine whether the
  estimated Hessian is positive definite (PD). If
  smallest_eigen > largest_eigen * tol than the Hessian is
  considered to be PD.}
}
\value{
  A matrix where each row contains a 95% confidence
  interval for each parameter. If the estimated Hessian is
  not PD then the boundaries of the CI will be set to NA.
}
\description{
  This function extract the final estimate for each
  parameter from object@resultPar and calculates confidence
  intervals (if the estimated Hessian is positive
  definite).
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com> #
}

