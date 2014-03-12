\name{checkNorm}
\alias{checkNorm}
\title{Checking the multivariate normality approximation.}
\usage{
  checkNorm(object, observed = NULL, param = object@param,
    nsim = 1000, cex.axis = 1, cex.lab = 1, ...)
}
\arguments{
  \item{object}{An object of class "synlik" or a matrix
  where each row is a random vector.}

  \item{param}{A vector of model's parameters at which the
  summary statistics will be simulated.}

  \item{nsim}{number of summary statistics to be simulated
  if object is of class "synlik", otherwise it is not
  used.}

  \item{cex.axis}{Axis scale expansion factor.}

  \item{cex.lab}{Axis label expansion factor.}

  \item{...}{additional arguments to be passed to
  object@simulator and object@summaries.  In general I
  would avoid using it and including in object@extraArgs
  everything they need.}
}
\value{
  Mainly produces plots and prints output. Also an array
  indicating proportion of simulated statistics smaller
  than observed.
}
\description{
  Given an object of class "synlik" this routine provides a
  graphical check of whether the distribution of the random
  summary statistics is multivariate normal.
}
\details{
  The method is from section 7.5 of Krzanowski (1988). The
  replicate vectors of summary statistic \code{S} are
  transformed to variables which should be univariate chi
  squared r.v.s with DoF given by the number of rows of
  \code{S}.  An appropriate QQ-plot is produced, and the
  proportion of the data differing substantially from the
  ideal line is reported. Deviations at the right hand end
  of the plot indicate that the tail behaviour of the
  Normal approximation is poor: in the context of synthetic
  likelihood this is of little consequence.  Secondly,
  \code{s} is transformed to a vector which should be
  i.i.d. N(0,1) under multivariate normality, and a QQ plot
  is produced. Unfortunately this approach is not very
  useful unless the dimension of \code{s} is rather large.
  In simulations, perfectly MVN data produce highly
  variable results, so that the approach lacks any real
  power.
}
\author{
  Simon N. Wood, maintained by Matteo Fasiolo
  <matteo.fasiolo@gmail.com>.
}
\references{
  Krzanowski, W.J. (1988) Principles of Multivariate
  Analysis. Oxford.
}

