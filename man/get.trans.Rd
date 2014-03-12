\name{get.trans}
\alias{get.trans}
\alias{trans.stat}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Obtain or apply a simple transformation to normality}
\description{
Given a matrix where each row is a set off iid draws from some distribution, finds a piecewise linear transformation of each row to acheive approximate normality.   
} %- end description

\usage{
get.trans(S)
trans.stat(S,trans)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{S}{a matrix, where each row contains i.i.d. observations of some random variable}

\item{trans}{a list, with an item for each row of \code{S}, which can be used to transform that row to approximate normality.}

}

\details{ 
Plots the quantiles of a standard normal against the ordered observations in each row of \code{S}. \code{\link{locator}} is then used to produce a piecewise linear approximation to this plot. Left clicking with the mouse addes a node to the piecwise approximation. Right clicking finishes (or indicates that no transformation is required).

The piecewise linear segments at each end are extended, each by 100 times the range of the data. The list returned by \code{get.trans} can then be used by \code{trans.stat} to apply the transformation. 
}
\value{ \code{get.trans} returns a list, with an item for each row of \code{S}. A \code{NULL} item implies that no transform should be applied, otherwise the item is a list with elements \code{x} and \code{y} describing a piecewise linear transform of the varaible.

\code{trans.stat} returns the transformed version of \code{S}.
}
\references{

\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  
MVN.check, robust.vcov
}

\examples{
\dontrun{
  library(sl)
  S <- matrix(0,2,100)
  S[1,] <- rpois(100,7)
  S[2,] <- rgamma(100,2)
  tr <- get.trans(S)
  S1 <- trans.stat(S,tr)
  par(mfrow=c(2,2))
  qqnorm(S[1,]);qqline(S[1,],col=2)
  qqnorm(S[2,]);qqline(S[2,],col=2)
  qqnorm(S1[1,]);abline(0,1,col=2)
  qqnorm(S1[2,]);abline(0,1,col=2)
}
}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





