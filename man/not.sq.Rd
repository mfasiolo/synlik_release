\name{not.sq}
\alias{not.sq}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Attenuated square function}
\description{
The normal pdf features a squared norm term. To make normal based likelihoods more robust, it can help
to replace the square with a function that grows much more slowly than the square, after some point.

} %- end description

\usage{
not.sq(x,alpha= .5,d0 = 5)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{x}{vector or scalar argument}

\item{alpha}{scalar determining exponent to replace 2, for \code{x>d0}.}

\item{d0}{scalar threshold. For \code{x} values below this value \code{x^2} is returned.}
}

\value{ For \code{x <= d0} returns \code{x^2}. For \code{x>d0} returns \code{k*x^alpha + c} where \code{c} and \code{k} are chosen to ensure continuity up to first derivative at \code{d0}. \code{x} is set to \code{abs(x)} at start, so that function is symmetric about zero.
}
\details{
The main use of this function is to replace the square function in the squared norm term in the normal likelihood. \code{d0} is 
set with reference to the appropriate chi-squared distribution. The substitution only has effect when the data are highly improbable under the model. For well fitting data the likelihood is unchanged.   
}

\author{ Simon N. Wood \email{simon.wood@r-project.org}}



\examples{
library(sl)
x <- seq(0,20,length=1000)
plot(x,not.sq(x,d0=5,alpha=.1),type="l")
lines(x,x^2,col=2)

}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





