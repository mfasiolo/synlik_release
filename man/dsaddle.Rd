\name{dsaddle}
\alias{dsaddle}
\title{Empirical saddlepoint density}
\usage{
  dsaddle(y, X, tol = 1e-06, decay = 0.5, deriv = FALSE,
    log = FALSE)
}
\arguments{
  \item{y}{Point at which the SPA is evaluated (d
  dimensional vector).}

  \item{X}{n by d matrix containing the data.}

  \item{tol}{Tolerance used to assess the convergence of
  the rootfinding routine used to fit the saddlepoint
  density. Default value is 1e-6.}

  \item{decay}{Rate at which the SPA falls back on a normal
  density. Should be a positive number, by default set to
  0.5.}

  \item{deriv}{If TRUE also the gradient of the
  log-saddlepoint density is returned.}

  \item{log}{If TRUE the log of the saddlepoint density is
  returned.}
}
\value{
  A list with entries: \itemize{ \item{ \code{llk} }{The
  value of the empirical saddlepoint at y;} \item{
  \code{mix} }{The mixture of normal-saddlepoint used (1
  means only saddlepoint);} \item{ \code{grad} }{The
  gradient of the log-density at y (optional);} }
}
\description{
  Gives a pointwise evaluation of the empirical saddlepoint
  and optionally of its gradient at position y
}
\author{
  Matteo Fasiolo <matteo.fasiolo@gmail.com> and Simon Wood.
}

