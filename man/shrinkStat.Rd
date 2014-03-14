\name{shrinkStat}
\alias{shrinkStat}
\title{Creates a summaries function using ridge regression or lasso}
\usage{
  shrinkStat(object, nsim, mu, sigma, type = "ridge",
    constr = list(), clean = TRUE, verbose = TRUE, ...)
}
\arguments{
  \item{object}{An object of class "synlik".}

  \item{nsim}{Number of summary statistics to be
  simulated.}

  \item{mu}{mean around which the parameters will be
  simulated.}

  \item{sigma}{covariance matrix used to simulate the
  parameters.}

  \item{...}{even more additional arguments to be passed to
  object@simulator and object@summaries.}
}
\value{
  A function that transforms simulated data in summary
  statistics. This function can by used as object@summaries
  function in an object of class "synlik".
}
\description{
  Fearnhead and Prangle (2012) proposed to use posterior
  expected values of the parameters as summary statistics.
  Here this goal is approximated by using fitted values of
  the parameters (obtained by ridge regression) as summary
  statistics.
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>.
}

