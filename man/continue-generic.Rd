\docType{methods}
\name{continue}
\alias{continue}
\alias{continue,smcmc-method}
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
\details{
  When \code{is("smcmc", object) == TRUE} continues MCMC
  estimation of an object of class \code{smcmc}. All input
  parameters are defaulted to the corresponding slots in
  the input object, with the exception of cluster. The
  chain restarts were it ended, burn-in is set to zero, the
  same prior (if any) is used.
}
\seealso{
  For examples, see \code{\link{smcmc-class}}.
}

