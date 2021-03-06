\name{smcmc-class}
\alias{smcmc-class}
\title{\code{smcmc-class}}
\description{
  Object representing the results of MCMC estimation on an
  object of class \code{synlik}, from which it inherits.  }

  \section{Slots}{ \describe{ \item{initPar}{Vector of
  initial parameters where the MCMC chain will start
  (\code{numeric}).} \item{niter}{Number of MCMC iterations
  (\code{integer}).} \item{nsim}{Number of simulations from
  the simulator at each step of the MCMC algorithm
  (\code{integer}).} \item{burn}{Number of initial MCMC
  iterations that are discarded (\code{integer}).}
  \item{priorFun}{Function that takes a vector of
  parameters as input and the log-density of the prior as
  output. If the output is not finite the proposed point
  will be discarded. (\code{function}).  The function needs
  to have signature \code{fun(x, ...)}, where \code{x}
  represents the input parameters (\code{function}).}
  \item{propCov}{Matrix representing the covariance matrix
  to be used to perturb the parameters at each step of the
  MCMC chain (\code{matrix}).} \item{targetRate}{Target
  rate for the adaptive MCMC sampler. Should be in (0, 1),
  default is NULL (no adaptation). The adaptation uses the
  approach of Vihola (2011). (\code{numeric})}
  \item{recompute}{If TRUE the synthetic likelihood will be
  evaluated at the current and proposed positions in the
  parameter space (thus doubling the computational effort).
  If FALSE the likelihood of the current position won't be
  re-estimated (\code{logical}).} \item{multicore}{If TRUE
  the \code{object@simulator} and \code{object@summaries}
  functions will be executed in parallel. That is the nsim
  simulations will be divided in multiple cores
  (\code{logical}).} \item{ncores}{Number of cores to use
  if multicore == TRUE (\code{integer}).}
  \item{accRate}{Acceptance rate of the MCMC chain, between
  0 and 1 (\code{numeric}).} \item{chains}{Matrix of size
  niter by length(initPar) where the i-th row contains the
  position of the MCMC algorithm in the parameter space at
  the i-th (\code{matrix}).} \item{llkChain}{Vector of
  niter elements where the i-th element is contains the
  estimate of the synthetic likelihood at the i-th
  iteration (\code{numeric}).} \item{control}{Control
  parameters used by the MCMC sampler: \itemize{
  \item{\code{theta} = controls the speed of adaption.
  Should be between 0.5 and 1.  A lower gamma leads to
  faster adaption.} \item{\code{adaptStart} = iteration
  where the adaption starts. Default 0.}
  \item{\code{adaptStop} = iteration where the adaption
  stops. Default \code{burn + niter}} \item{\code{saveFile}
  = path to the file where the intermediate results will be
  stored (ex: "~/Res.RData").} \item{\code{saveFreq} =
  frequency with which the intermediate results will be
  saved on \code{saveFile}.  Default 100.}
  \item{\code{verbose} = if \code{TRUE} intermediate
  posterior means will be printer.} \item{\code{verbFreq} =
  frequency with which the intermediate posterior means
  will be printer. Default 500.} } } }
}
\examples{
# Load "synlik" object
data(ricker_sl)

plot(ricker_sl)

# MCMC estimation
set.seed(4235)
ricker_sl <- smcmc(ricker_sl,
                   initPar = c(3.2, -1, 2.6),
                   niter = 50,
                   burn = 3,
                   priorFun = function(input, ...) 1,
                   propCov = diag( c(0.1, 0.1, 0.1) )^2,
                   nsim = 200,
                   multicore = FALSE)

# Continue with additional 50 iterations
ricker_sl <- continue(ricker_sl, niter = 50)

plot(ricker_sl)
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}
\references{
  Vihola, M. (2011) Robust adaptive Metropolis algorithm
  with coerced acceptance rate.  Statistics and Computing.
}

