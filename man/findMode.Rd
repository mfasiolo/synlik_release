\name{findMode}
\alias{findMode}
\title{Finding the mode of the empirical saddelpoint density}
\usage{
  findMode(X, init = NULL, decay = 0.5, method = "BFGS",
    sadTol = 1e-06, marginal = FALSE, ...)
}
\arguments{
  \item{X}{n by d matrix containing the data.}

  \item{init}{d-dimensional vector containing the starting
  point for the optimization. By default it is equal to
  colMeans(X).}

  \item{decay}{Rate at which the SPA falls back on a normal
  density. Should be a positive number, by default set to
  0.5.}

  \item{method}{Optimization method used by stats::optim(),
  see ?optim for details. By default it is "BFGS".}

  \item{sadTol}{Tolerance used to assess the convergence of
  the rootfinding routine used to fit the saddlepoint
  density. Default value is 1e-6.}

  \item{...}{Extra arguments to be passed to the
  optimization routine stats::optim.}
}
\value{
  A list where \code{mode} is the location of mode of the
  empirical saddlepoint density, while the other entries
  are the same as for stats::optim.
}
\description{
  Given a sample from a d-dimensinal distribution the
  routine find the mode of the corresponding empirical
  saddlepoint density.
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>.
}

