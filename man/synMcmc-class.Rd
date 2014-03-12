\name{synMcmc-class}
\alias{synMcmc-class}
\title{\code{synMcmc-class}}
\usage{
  .synMcmc(...)
}
\description{
  \describe{ \item{initPar}{Vector of initial parameters
  where the MCMC chain will start (\code{numeric}).}
  \item{nIter}{Number of MCMC iterations (\code{numeric}).}
  \item{burnIn}{Number of initial MCMC iterations which are
  discarded (\code{numeric}).} \item{priorFun}{Function
  that takes a vector of parameters as input and gives
  either 0 or 1 as output (uniform prior)
  (\code{function}).} \item{propCov}{Matrix representing
  the covariance matrix to be used to perturb the
  parameters at each step of the MCMC chain
  (\code{matrix}).} \item{nsim}{Number of simulations from
  the simulator at each step of the MCMC algorithm
  (\code{numeric}).} \item{sameSeed}{If TRUE the synthetic
  likelihood will be evaluated at the current and proposed
  positions in the parameter space using the same seed
  (thus doubling the computational effort). If FALSE the
  likelihood of the current position won't be re-estimated
  (\code{logical}).} \item{multicore}{If TRUE the
  object@simulator and object@summaries functions will be
  executed in parallel. That is the nsim simulations will
  be divided in multiple cores (\code{logical}).}
  \item{ncores}{Number of cores to use if multicore == TRUE
  (\code{numeric}).} \item{acceptRate}{Acceptance rate of
  the MCMC chain, between 0 and 1 (\code{numeric}).}
  \item{mcmcChain}{Matrix of size nIter by length(initPar)
  where the i-th row contains the position of the MCMC
  algorithm in the parameter space at the i-th
  (\code{matrix}).} \item{LoglikChain}{Vector of nIter
  elements where the i-th element is contains the estimate
  of the synthetic likelihood at the i-th iteration
  (\code{numeric}).} }
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

