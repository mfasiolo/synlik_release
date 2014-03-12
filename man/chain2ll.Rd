\name{chain2ll}
\alias{chain2ll}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Get likelihood asypmtotics information form MCMC output}
\description{
Takes output from an MCMC run using a synthetic likelihood, and estimates quadratic 
approximation to the log likelihood, near the maximum. Returns corresponding MLE, covariance matrix and 
standard errors.
} %- end description

\usage{
chain2ll(th,para=NULL,ll="ll",start=2000)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{th}{a matrix, each row of which contains MCMC output for one parameter, or the 
synthetic log likelihood itself. rows must be named. }
\item{para}{Optional list of parameter names, allowing \code{th} to contain more varaibles than
will be used for the quadratic approximation process. If \code{NULL} then the parameter names are 
obtained from the row names of \code{th}.}
\item{ll}{Name of the row of \code{th} containing the synthetic log likelihood values.}
\item{start}{Which column of \code{th} to start at. Columns before this one are discarded.}
}

\details{ 
Uses least squares to fir a model to the data in which synthetic log likelihood is the response, and the predictor is quadratic in the model parameters. From this fitted model the MLE of the the parameters and the approximate hessian of the (denoised) synthetic log likelihood can be obtained.

}
\value{list containing the MLE standard errors, Hessian and it's negative generalized inverse --- 
the approximate covariance matrix.
}
\references{


\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  
ricker
}

\examples{
## see ?ricker
}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





