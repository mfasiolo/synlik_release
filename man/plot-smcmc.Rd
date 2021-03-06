\name{plot-smcmc}
\alias{plot,smcmc,missing-method}
\alias{plot-smcmc}
\title{Method for plotting an object of class \code{smcmc}.}
\arguments{
  \item{x}{An object of class \code{smcmc}.}

  \item{trans}{Name list or vector containing names of
  transforms for some parameters (ex: \code{list("par1" =
  "exp", "par2" = "log")}).  The transformations will be
  applied before plotting.}

  \item{addPlot1}{Name of additional plotting function that
  will be call after the MCMC chain have been plotted. It
  has to have prototype \code{fun(nam, ...)} where
  \code{nam} will be the parameter name. See "examples".}

  \item{addPlot2}{Name of additional plotting function that
  will be call after the histograms have been plotted. It
  has to have prototype \code{fun(nam, ...)} where
  \code{nam} will be the parameter name. See "examples".}

  \item{...}{additional arguments to be passed to the
  plotting functions.}
}
\value{
  NULL
}
\description{
  Method for plotting an object of class \code{smcmc}.
}
\examples{
data(ricker_smcmc)

# Functions for additional annotations (true parameters)
addline1 <- function(parNam, ...){
               abline(h = exp(ricker_smcmc@param[parNam]), lwd = 2, lty = 2, col = 3)
               }
addline2 <- function(parNam, ...){
               abline(v = exp(ricker_smcmc@param[parNam]), lwd = 2, lty = 2, col = 3)
               }

# Transformations (exponentials)
trans <- rep("exp", 3)
names(trans) <- names(ricker_smcmc@param)

plot(ricker_smcmc,
     trans = trans,
     addPlot1 = "addline1",
     addPlot2 = "addline2")
}
\seealso{
  \code{\link{smcmc-class}}, \code{\link{plot}}.
}

