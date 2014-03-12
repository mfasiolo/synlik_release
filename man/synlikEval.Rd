\name{synlikEval}
\alias{synlikEval}
\title{Evaluate the synthetic likelihood.}
\usage{
  synlikEval(object, param, nsim, saddle = FALSE,
    decay = 0.5, multicore = FALSE,
    ncores = detectCores() - 1, cluster = NULL, ...)
}
\arguments{
  \item{object}{An object of class "synlik".}

  \item{param}{Vector of parameters at which the synthetic
  likelihood will be evaluated.}

  \item{nsim}{Number of simulation from the model.}

  \item{saddle}{If FALSE the distribution of the summary
  statistics will be approximated by a normal, if TRUE an
  Empirical Saddlepoint approximation will be used.}

  \item{decay}{Useful only if \code{saddle} == TRUE. It is
  the rate at which the saddlepoint density falls back to
  the normal density as the distance between observed and
  simulated statistics increases.}

  \item{multicore}{(logical) if TRUE the object@simulator
  and object@summaries functions will be executed in
  parallel. That is the nsim simulations will be divided in
  multiple cores.}

  \item{ncores}{(integer) number of cores to use if
  multicore == TRUE.}

  \item{cluster}{an object of class c("SOCKcluster",
  "cluster"). This allowes the user to pass her own
  cluster, which will be used if multicore == TRUE. The
  user has to remember to stop the cluster.}

  \item{...}{additional arguments to be passed to
  object@simulator and object@summaries.  In general I
  would avoid using it and including in those two function
  everything they need.}
}
\value{
  A list where "logLik" is the log of the estimated
  synthetic likelihood and "mix" is the share of
  saddlepoint used (0 means no saddlepoint (only normal)
  and 1 means only saddlepoint (no normal) ).
}
\description{
  Evaluate the synthetic likelihood.
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com>
}

