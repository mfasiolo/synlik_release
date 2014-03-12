\name{order.dist}
\alias{order.dist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summarize marginal distribution of (differenced) series}
\description{
Summarizes (difference) distribution of replicate series, by regressing ordered differenced series on a reference series.
} %- end description

\usage{
order.dist(x,z,np,diff)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{x}{a matrix. Each column contains a replicate series. }

\item{z}{array of lags, for rhs terms.}

\item{np}{maximum power on rhs of regression}

\item{diff}{order of differencing (zero for none)}

}

\details{ 
 Experimental

}
\value{Currently a matrix of the coefs. per series.
}
\references{


\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  
sl.acf, nlar
}

\examples{
  library(sl)
  set.seed(10)
  n <- 100;nr <- 3
  x <- matrix(runif(n*nr),n,nr)
  z <- runif(n)
  beta <- order.dist(x,z,np=3,diff=1)
  
  zd <- z;xd <- x[,3]
  zd <- diff(zd,1);xd <- diff(xd,1)
  zd <- sort(zd);zd <- zd - mean(zd)
  xd <- sort(xd);xd <- xd - mean(xd)
  lm(xd~zd+I(zd^2)+I(zd^3)-1)

}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





