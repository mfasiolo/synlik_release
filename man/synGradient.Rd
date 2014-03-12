\name{synGradient}
\alias{synGradient}
\title{Estimate gradient and hessian of the synthetic likelihood.}
\usage{
  synGradient(object, param, nsim, covariance,
    addRegr = TRUE, constr = list(), multicore = FALSE,
    ncores = detectCores() - 1, cluster = NULL,
    tolVar = 10 * .Machine$double.eps, ...)
}
\arguments{
  \item{object}{an object of class "synlik".}

  \item{param}{(numeric) vector containing the current
  values of the model's parameters.}

  \item{nsim}{(integer) number of simulated statistics.}

  \item{covariance}{(matrix) covariance matrix used to
  simulate the parameters.}

  \item{addRegr}{(logical).  If FALSE the statistics
  calculated by object@summaries will be used (SL
  approach); if TRUE the simulated parameters will be
  regressed on the statistics and the fitted values of the
  paramaters given the _observed_ statistics will be used
  as statistics (referred to as SL+ approach).}

  \item{constr}{(named list) used to impose constraints on
  the parameters.  Composed of 3 elements: [["indexes"]] =
  (numeric integers) indexes of the elements to check;
  [["upper"]] = (numeric) upper bounds for the elements in
  "indexes"; [["lower"]] = (numeric) lower bounds for the
  elements in "indexes".}

  \item{multicore}{(logical) if TRUE the object@simulator
  and object@summaries functions will be executed in
  parallel. That is the nsim simulations will be divided in
  multiple cores.}

  \item{ncores}{(integer) number of cores to use if
  multicore == TRUE.}

  \item{cluster}{an object of class c("SOCKcluster",
  "cluster"). This allowes the user to pass her own
  cluster, which will be used if multicore == TRUE. The
  user has to remember to stop the cluster.}

  \item{...}{additional arguments to be passed to
  object@simulator and object@summaries.  In general I
  would avoid using it and including in those two function
  everything they need.}
}
\value{
  a list containing: ["gradient"] = (numeric) estimated
  gradient of the synthetic log-likelihood at currPar;
  ["hessian"] = (matrix) estimated hessian of the synthetic
  log-likelihood at currPar; ["llk"] = (scalar) estimated
  value of the synthetic log-likelihood at currPar.
}
\description{
  Estimate gradient and hessian of the synthetic
  likelihood.
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

