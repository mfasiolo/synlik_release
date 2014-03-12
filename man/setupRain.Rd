\name{setupRain}
\alias{setupRain}
\title{Set up the integrated thresholded rainfall as in Bhadhra et al, 2011.}
\usage{
  setupRain(rainfall, rainYear, nSteps, obsInterval = 1/12,
    lag = 5, v = 200)
}
\arguments{
  \item{rainfall}{Vector time series of monthly observed
  precipitations.}

  \item{rainYear}{Time of precipitation [0, 20)}

  \item{nSteps}{Number of multinomial-euler step in which
  each month is divided.}

  \item{obsInterval}{Interval of time between 2
  observations (typically 1/12 for monthly data).}

  \item{lag}{Number of month used as lag (see Bhadhra et al
  2011).}

  \item{v}{Threshold used for the rainfall (see Bhadhra et
  al 2011).}
}
\value{
  vector of length( length(rainfall) * nSteps ) containing
  integrated rainfall.
}
\description{
  Set up the integrated thresholded rainfall as in Bhadhra
  et al, 2011.
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

