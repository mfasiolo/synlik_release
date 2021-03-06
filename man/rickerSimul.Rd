\name{rickerSimul}
\alias{rickerSimul}
\title{Simulates from the ricker model}
\usage{
  rickerSimul(param, nsim, extraArgs, ...)
}
\arguments{
  \item{param}{vector of log-parameters: logR, logSigma,
  logPhi. Alternatively a matrix \code{nsim} by 3 were each
  row is a different parameter vector.}

  \item{nsim}{Number of simulations from the model.}

  \item{extraArgs}{A named list of additional arguments:
  \itemize{ \item{\code{nObs} = Length of each simulated
  time series.} \item{\code{nBurn} = Number of initial
  steps to be discarded before saving the following
  \code{nObs} steps.} \item{\code{randInit} = if
  \code{TRUE} (default) the initial state N0 is
  \code{runif(0, 1)}, otherwise it is equal to
  \code{extraArgs$initVal}.} \item{\code{initVal} = initial
  value N0, used only if \code{extraArgs$randInit ==
  TRUE}.} }}

  \item{...}{Need for compatibility with \code{synlik}, but
  not used.}
}
\value{
  A matrix \code{nsim} by \code{nObs}, where each row is a
  simulated path.
}
\description{
  Simulator for the stochastic Ricker model, as described
  by Wood (2010). The observations are Y_t ~ Pois(Phi *
  N_t), and the dynamics of the hidden state are given by
  N_t = r * N_\{t-1\} * exp( -N_\{t-1\} + e_t ), where e_t
  ~ N(0, Sigma^2).
}
\examples{
tmp <- rickerSimul(c(3.8, -1.2, 2.3), nsim = 2, extraArgs = list("nObs" = 50, "nBurn" = 200))
matplot(t(tmp), type = 'l', ylab = "Y", xlab = "Time")

parMat <- rbind(c(3.8, -1.2, 2.3),  # Chaotic
                c(2.5, -1.2, 2.3))  # Not Chaotic

par(mfrow = c(2, 1))
tmp <- rickerSimul(parMat, nsim = 2, extraArgs = list("nObs" = 50, "nBurn" = 200))
plot(tmp[1, ], type = 'l', ylab = "Y", xlab = "Time")
plot(tmp[2, ], type = 'l', ylab = "Y", xlab = "Time")
}
\author{
  Simon Wood and Matteo Fasiolo <matteo.fasiolo@gmail.com>.
}
\references{
  Simon N Wood. Statistical inference for noisy nonlinear
  ecological dynamic systems. Nature, 466(7310):1102--1104,
  2010. \cr \cr
}
\seealso{
  \link{ricker_sl}
}

