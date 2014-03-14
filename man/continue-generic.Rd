\name{continue}
\alias{continue}
\title{Continuing estimation.}
\arguments{
  \item{object}{An object representing the results of an
  estimation procedure which we wish to continue.  For
  example it might represents an MCMC chain.}

  \item{...}{Additional arguments that might be needed to
  continue the estimation procedure.}
}
\value{
  An object of the same class as \code{object}, where the
  results of the estimation have been updated.
}
\description{
  Generic function, that given the results of an estimation
  procedure (ex. MCMC or maximum likelihood optimization)
  continues the procedure for some more iterations.
}
\seealso{
  For examples, see \code{\link{smcmc-class}}.
}

