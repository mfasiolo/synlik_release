\name{nlar}
\alias{nlar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Turn time series into stats}
\description{
Polynomial autoregression for multiple reps.
} %- end description

\usage{
nlar(x,lag,power)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{x}{a matrix. Each column contains a replicate series. }

\item{lag}{array of lags, for rhs terms.}

\item{power}{array of powers, for rhs terms.}

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
order.dist, sl.acf
}

\examples{
  library(sl)
  set.seed(10)
  x <- matrix(runif(200),100,2)
  beta <- nlar(x,lag=c(1,1),power=c(1,2))
  y <- x[,1]
  y <- y - mean(y)
  z <- y[1:99];y <- y[2:100]
  lm(y~z+I(z^2)-1)
  beta

  ## NA testing
  x[5,1] <- x[45,2] <- NA
  beta <- nlar(x,lag=c(1,1),power=c(1,2))
  y <- x[,1]
  y <- y - mean(y,na.rm=TRUE)
  z <- y[1:99];y <- y[2:100]
  lm(y~z+I(z^2)-1)
  beta

  ## higher order...
  set.seed(10)
  x <- matrix(runif(100),100,2)
  beta <- nlar(x,lag=c(6,6,6,1,1),power=c(1,2,3,1,2))
  k <- 2
  y <- x[,k]
  y <- y - mean(y)
  ind <- (1+6):100
  y6 <- y[ind-6];y1 <- y[ind-1];y <- y[ind]
  beta0 <- coef(lm(y~y6+I(y6^2)+I(y6^3)+y1+I(y1^2)-1))
  as.numeric(beta[,k]);beta0;beta0-as.numeric(beta[,k])
  
}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





