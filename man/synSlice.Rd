\name{synSlice}
\alias{synSlice}
\title{Plot slices of the synthetic log-likelihood.}
\usage{
  synSlice(object, ranges, nsim, param = object@param,
    pairs = FALSE, draw = TRUE, trans = NULL,
    multicore = FALSE, ncores = detectCores() - 1,
    cluster = NULL, ...)
}
\arguments{
  \item{object}{("synlik") object.}

  \item{ranges}{ranges of values along which we want the
  slices. If \code{length(parName) == 1} than range has a
  vector, while if \code{length(parName) == 2} it have to
  be a named list of 2 vectors (ex: \code{list("alpha" =
  1:10, "beta" = 10:1)}).}

  \item{nsim}{Number of simulations used to evaluate the
  synthetic likelihood at each location.}

  \item{param}{Named vector containing the value of the ALL
  parameters (including the sliced one). Parameters that
  are not in \code{parName} will be fixed to the values in
  \code{param}.}

  \item{draw}{If TRUE the slice will be plotted.}

  \item{trans}{Named vector of transformations to be
  applied to the parameters in \code{parName} before
  plotting {ex: \code{trans = c(s = "exp", d = "exp")}}/}

  \item{multicore}{If TRUE the object@simulator and
  object@summaries functions will be executed in parallel.
  That is the nsim simulations will be divided in multiple
  cores.}

  \item{ncores}{Number of cores to use if multicore ==
  TRUE.}

  \item{cluster}{An object of class c("SOCKcluster",
  "cluster"). This allowes the user to pass her own
  cluster, which will be used if multicore == TRUE. The
  user has to remember to stop the cluster.}

  \item{...}{additional arguments to be passed to
  synlikEval, see \code{\link{synlikEval}}.}
}
\value{
  Either a vector or matrix of log-synthetic likelihood
  estimates, depending on whether \code{length(parNames)
  ==} 1 or 2.
}
\description{
  Plot slices of the synthetic log-likelihood.
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

