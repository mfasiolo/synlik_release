\name{synMaxlik}
\alias{synMaxlik}
\title{Synthetic likelihood maximization routine.}
\usage{
  synMaxlik(object, nIter, nsim, initCov,
    initPar = object@param, addRegr = TRUE,
    constr = list(), control = list(), multicore = FALSE,
    ncores = detectCores() - 1, cluster = NULL,
    verbose = FALSE, ...)
}
\arguments{
  \item{object}{("synlik") object.}

  \item{nIter}{(integer) numer of iterations.}

  \item{nsim}{(integer) numer of simulations from the model
  at each step.}

  \item{initCov}{(matrix) initial covariance matrix used to
  simulate the paramters at each step.}

  \item{initPar}{(numeric) vector of initial values of the
  parameters.}

  \item{addRegr}{(logical) if FALSE the statistics
  calculated by object@summaries will be used (SL approach)
  if TRUE the simulated parameters will be regressed on the
  statistics and the fitted values of the paramaters given
  the _observed_ statistics will be used as statistics (SL+
  approach)}

  \item{constr}{(named list) of 3 elements: [["indexes"]] =
  (numeric integers) indexes of the elements to check;
  [["upper"]] = (numeric) upper bounds for the elements in
  "indexes"; [["lower"]] = (numeric) lower bounds for the
  elements in "indexes".}

  \item{control}{Named list of control setting for the
  optimization routine.}

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

  \item{verbose}{(logical) if TRUE lots of things will be
  printed.}

  \item{...}{additional arguments to be passed to
  object@simulator and object@summaries.  In general I
  would avoid using it and including in those two function
  everything they need.}
}
\value{
  object of class "synOptim": see info about that class.
}
\description{
  Synthetic likelihood maximization routine.
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

