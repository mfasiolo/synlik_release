\name{confint.sml}
\alias{confint.sml}
\title{Parameters covariance matrix for "sml" objects}
\usage{
  \method{confint}{sml} (object, parm, level = 0.95,
    lag = 10, ...)
}
\arguments{
  \item{object}{a fitted model object of class "sml".}

  \item{parm}{a specification of which parameters are to be
  given confidence intervals, either a vector of numbers or
  a vector of names.  If missing, all parameters are
  considered.}

  \item{level}{the confidence level required.}

  \item{lag}{(integer) final estimate of the parameter is
  the mean of the last "lag" iterations. Hence lag <=
  object@niter.}
}
\description{
  Computes confidence intervals for one or more parameters
  in a fitted model of class "sml".
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

