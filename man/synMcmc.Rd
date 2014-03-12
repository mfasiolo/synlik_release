\name{synMcmc}
\alias{synMcmc}
\title{MCMC parameter estimation for "synlik" objects.}
\usage{
  synMcmc(object, initPar, nIter, nsim, propCov,
    burnIn = 0, priorFun = function(param) {     1 },
    targetRate = NULL, recompute = FALSE,
    multicore = !is.null(cluster), cluster = NULL,
    control = list(), ...)
}
\arguments{
  \item{object}{An object of class "synlik".}

  \item{initPar}{Vector of initial parameters where the
  MCMC chain will start.}

  \item{nIter}{Number of MCMC iterations.}

  \item{burnIn}{Number of initial MCMC iterations which
  will be discarded.}

  \item{priorFun}{Function that takes a vector of
  parameters as input and gives either 0 or 1 as output
  (uniform prior).}

  \item{propCov}{Matrix representing the covariance matrix
  to be used to perturb the parameters at each step of the
  MCMC chain.}

  \item{nsim}{Number of simulated statistics at each step.}

  \item{sameSeed}{If TRUE the synthetic likelihood will be
  evaluated at the current and proposed positions in the
  parameter space using the same seed (thus doubling the
  computational effort). If FALSE the likelihood of the
  current position won't be re-estimated.}

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
  synlikEval, see \code{\link{synlikEval}}.}
}
\value{
  An object of class "synMcmc".
}
\description{
  MCMC parameter estimation for "synlik" objects.
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

