\name{continue.smcmc}
\alias{continue.smcmc}
\alias{plot,smcmc,missing-method}
\title{Continuing estimation of an \code{smcmc} object.}
\usage{
  continue.smcmc(object, niter = object@niter,
    nsim = object@nsim, propCov = object@propCov,
    targetRate = object@targetRate,
    recompute = object@recompute,
    multicore = object@multicore, ncores = object@ncores,
    cluster = NULL, control = object@control, ...)
}
\arguments{
  \item{object}{An object representing the results of an
  estimation procedure which we wish to continue.  For
  example it might represents and MCMC chain.}

  \item{niter}{Number of additional MCMC iterations
  requested. The output chain will have length object}

  \item{nsim}{see \code{\link{smcmc-class}}.}

  \item{propCov}{see \code{\link{smcmc-class}}.}

  \item{targetRate}{see \code{\link{smcmc-class}}.}

  \item{recompute}{see \code{\link{smcmc-class}}.}

  \item{multicore}{see \code{\link{smcmc-class}}.}

  \item{ncores}{see \code{\link{smcmc-class}}.}

  \item{cluster}{see \code{\link{smcmc}}.}

  \item{control}{see \code{\link{smcmc-class}}.}

  \item{...}{additional arguments to be passed to
  \code{slik()}, see \code{\link{slik}}.}
}
\value{
  An object of the same class as \code{object}, where the
  results of the estimation have been updated.
}
\description{
  Continues MCMC estimation of an object of class
  \code{smcmc}. All input parameters are defaulted to the
  corresponding slots in the input object, with the
  exception of cluster. The chain restarts were it ended,
  burn-in is set to zero, the same prior (if any) is used.
}
\examples{
# For an example see help("smcmc-class").
}
\seealso{
  \code{\link{smcmc-class}}, \code{\link{smcmc}},
  \code{\link{continue}}.
}

