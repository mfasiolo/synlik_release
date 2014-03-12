\name{robust.vcov}
\alias{robust.vcov}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Robust covariance matrix estimation}
\description{
Obtains a robust estimate of the covariance matrix of a sample of multivariate data, using Campbell's (1980) method
as described on p231-235 of Krzanowski (1988). 
} %- end description

\usage{
robust.vcov(sY,alpha=2,beta=1.25)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{sY}{a matrix, where each column is a replicate observation on a multivariate r.v.}

\item{alpha}{tuning parameter, see details}

\item{beta}{tuning parameter, see details.}

}

\details{ 
Campbell (1980) suggests an estimator of the covariance matrix which downweights observations at more than some Mahalanobis distance \code{d.0} from the mean. \code{d.0} is \code{sqrt(nrow(sY))+alpha/sqrt(2)}. Weights are one for observations with Mahalanobis distance, \code{d}, less than \code{d.0}. Otherwise weights are are \code{d.0*exp(-.5*(d-d.0)^2/beta)/d}. The defaults are as recommended by Campbell.

This routine also uses pre-conditioning to ensure good scaling and stable numerical calculations.

}
\value{ A list. \code{E} is a square root of the inverse covariance matrix. i.e. the inverse cov matrix is \code{t(E)\%*\%E}.
\code{half.ldet.V} is half the log of the determinant of the covariance matrix. \code{mY} is the estimated mean, and 
\code{sd} the estimated standard deviations of each statistic.
}
\references{
Krzanowski, W.J. (1988) Principles of Multivariate Analysis. Oxford.

Campbell, N.A. (1980) Robust procedures in multivariate analysis I: robust covariance estimation. JRSSC 29, 231-237. 

\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  
get.trans, MVN.check
}

\examples{
  library(sl)
  p <- 5;n <- 100
  Y <- matrix(runif(p*n),p,n)
  robust.vcov(Y)
  
}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





