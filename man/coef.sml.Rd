\name{coef.sml}
\alias{coef.sml}
\title{Parameters estimates for "sml" objects}
\usage{
  \method{coef}{sml} (object, lag = 10, ...)
}
\arguments{
  \item{object}{("sml") object.}

  \item{lag}{(integer) final estimate of the parameter is
  the mean of the last "lag" iterations. Hence lag <
  object@niter.}
}
\description{
  Parameters estimates for "sml" objects
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

