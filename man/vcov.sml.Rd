\name{vcov.sml}
\alias{vcov.sml}
\title{Parameters covariance matrix for "sml" objects}
\usage{
  \method{vcov}{sml} (object, nreps = 1000, boot = TRUE,
    ...)
}
\arguments{
  \item{object}{("sml") object.}

  \item{nreps}{(integer) number of simulations used to tilt
  the negative Hessian (-H) toward positive definiteness.
  Used only if -H is not PD.}

  \item{boot}{(logical) relevant only if -H is not PD.  If
  TRUE hessians will be simulated by resampling parameters
  and likelihoods.  If FALSE hessians will be simulated the
  asymptotic distribution of the regression coefficients.}
}
\description{
  Parameters covariance matrix for "sml" objects
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

