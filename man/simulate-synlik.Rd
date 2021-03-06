\name{simulate}
\alias{simulate}
\alias{simulate,synlik-method}
\title{Simulate data or statistics from an object of class \code{synlik}.}
\arguments{
  \item{object}{An object of class \code{synlik}.}

  \item{nsim}{Number of simulations from the model.}

  \item{seed}{Random seed to be used. It is not passed to
  the simulator, but simply passed to \code{set.seed()}
  from within \code{simulate.synlik}.}

  \item{param}{Vector of parameters passed to
  \code{object@simulator}.}

  \item{stats}{If \code{TRUE} the function trasforms the
  simulated data into statistics using
  \code{object@summaries}.}

  \item{clean}{If \code{TRUE} the function tries to clean
  the statistics from NaNs or non-finite values.  Given
  that \code{object@summaries} has to returns a numeric
  vector or a matrix where each row is a simulation, rows
  containing non-finite values will be discarded.}

  \item{verbose}{If \code{TRUE} the function will complain
  if, for instance, the simulations contain lots of
  non-finite values.}

  \item{...}{additional arguments to be passed to
  \code{object@simulator} and \code{object@summaries}.  In
  general I would avoid using it and including
  \code{object@extraArgs} everything they need.}
}
\value{
  If \code{stats == FALSE} the output will that of
  \code{object@simulator}, which depends on the simulator
  used by the user.  If \code{stats == TRUE} the output
  will be a matrix where each row is vector of simulated
  summary statistics.
}
\description{
  Simulate data or statistics from an object of class
  \code{synlik}.
}
\examples{
data(ricker_sl)

# Simulate data
simulate(ricker_sl, nsim = 2)

# Simulate statistics
simulate(ricker_sl, nsim = 2, stats = TRUE)
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}
\seealso{
  \code{\link{synlik-class}}, \code{\link{simulate}}.
}

