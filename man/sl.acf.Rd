\name{sl.acf}
\alias{sl.acf}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Turn time series into stats}
\description{
Takes replicate time series and summarizes them, quickly.
} %- end description

\usage{
sl.acf(x,max.lag=10)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{x}{a matrix. Each column contains a replicate series. }
\item{max.lag}{How many lags to use}
}

\details{ 
 Experimental

}
\value{Currently a matrix of the acf per series.
}
\references{


\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  
nlar, order.dist
}

\examples{
  library(sl)
  set.seed(10)
  x <- matrix(runif(1000),100,10)
  acf <- sl.acf(x)
}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





