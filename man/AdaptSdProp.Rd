\name{AdaptSdProp}
\alias{AdaptSdProp}
\title{adapt the standard deviation of the proposal}
\usage{
  AdaptSdProp(sdProp, accepts, lowAccept = 0.15,
    highAccept = 0.4, lag = 20)
}
\arguments{
  \item{sdProp}{current standard deviation for the
  proposal}

  \item{accepts}{binary vector (TRUE/FALSE or 0/1) with 1
  or TRUE for accepted values}

  \item{lowAccept}{minimum acceptance rate}

  \item{highAccept}{maximum acceptance rate}

  \item{lag}{number of past iterations took into account}
}
\value{
  a scalar: the new standard deviation for proposals, with
  an attributes(sdProp)$noupdate set to TRUE if no update
  needed
}
\description{
  From a binary vector indicating acceptance, adapt the
  sampling deviation so that acceptance rate fall by
  default within 0.15 and 0.4 (see Gelman1996).  Be careful
  not to begin the sampling before that)
}
\seealso{
  MCMC
}

