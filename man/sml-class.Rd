\name{sml-class}
\alias{sml-class}
\title{\code{sml-class}}
\usage{
  .sml(...)
}
\description{
  \describe{ \item{initPar}{Vector of initial parameters
  from which the optimization starts (\code{numeric}).}
  \item{nIter}{Number of iterations of the optimization
  routine (\code{numeric}).} \item{nsim}{Number of
  simulations from the simulator at each step of the
  optimization routine (\code{numeric}).}
  \item{initCov}{Covariance matrix used to simulate the
  parameter at each step of the optimization routine
  (\code{matrix}).} \item{constr}{Named list of contraints
  on the parameters, it has 3 elements: [["indexes"]] =
  (numeric integers) indexes of the elements to check;
  [["upper"]] = (numeric) upper bounds for the elements in
  "indexes"; [["lower"]] = (numeric) lower bounds for the
  elements in "indexes" (\code{list}).} \item{multicore}{If
  TRUE the \code{object@simulator} and
  \code{object@summaries} functions will be executed in
  parallel. That is the nsim simulations will be divided in
  multiple cores (\code{logical}).} \item{ncores}{Number of
  cores to use if multicore == TRUE (\code{numeric}).}
  \item{estim}{Matrix typically nIter by length(initPar)
  where the i-th row contains the estimate of the
  parameters at the i-th iteration (\code{matrix}).}
  \item{simPar}{Matrix typically nIter * nP by
  length(initPar) where the i-th row contains the estimate
  of the parameters at the i-th iteration (\code{matrix}).}
  \item{resultHess}{List of nIter elements where the i-th
  element is contains the estimate of the Hessian of the
  synthetic likelihood at the i-th iteration
  (\code{list}).} \item{resultCovar}{List of nIter elements
  where the i-th element is contains the estimate of the
  covariance matrix of the parameters at the i-th iteration
  (\code{list}).} \item{resultLoglik}{Vector of nIter
  elements where the i-th element is contains the estimate
  of the synthetic likelihood at the i-th iteration
  (\code{numeric}).} }
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

