\name{getRidgeStat}
\alias{getRidgeStat}
\title{Creates a summaries function by ridge regression}
\usage{
  getRidgeStat(object, nsim, upper, lower, ...)
}
\arguments{
  \item{object}{An object of class "synlik".}

  \item{nsim}{Number of summary statistics to be
  simulated.}

  \item{upper}{Vector giving the upper limit of the
  hypercube where the statistics will be simulated
  uniformly.}

  \item{lower}{Vector giving the upper limit of the
  hypercube where the statistics will be simulated
  uniformly.}

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

