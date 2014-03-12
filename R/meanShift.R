#####################
### Mean-shift algorithm
#####################
#####
#' Mean-shift mode seeking algorithm
#' @description Given a sample from a d-dimensinal distribution, an initialization point and a bandwith
#'              the algorithm finds the nearest mode of the corresponding Gaussian kernel density.
#' @param X n by d matrix containing the data.
#' @param init d-dimensional vector containing the initing point for the optimization. By default
#'             it is equal to colMeans(X).
#' @param H Positive definite bandwidth matrix representing the covariance of 
#'          each component of the Gaussian kernel density.
#' @param tol Tolerance used to assess the convergence of the algorithm, which is stopped if the absolute values
#'            of increments along all the dimensions are smaller then tol at any iteration. Default value is 1e-6.
#' @param traj If FALSE only the latest iteration is returned, if TRUE the function will return a matrix where
#'             the i-th row is the position of the algorithms at the i-th iteration.
#' @return If traj == FALSE only the latest iteration is returned, if traj == TRUE the function will return a matrix where
#'             the i-th row is the position of the algorithms at the i-th iteration.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @export
#'

meanShift <- function(X, init, H, tol = 1e-6, traj = FALSE)
{
  if(is.matrix(X) == FALSE) X <- matrix(X, length(X), 1)
  if( !is.matrix(H) ) H <- diag(H, ncol(X))
  
  d <- length(init)
  N <- nrow(X)
  oldPos <- currPos <- trajectory <- init
  delta <- rep(2, d) * tol
  
  weights <- numeric(N)
  cholDec <- chol(H)
  
  while( any( delta > tol ) )
  {
    oldPos <- currPos
    weights <- dmvnFast(X = X, mu = oldPos, sigma = cholDec, isChol = TRUE)
    currPos <- weights %*% X / sum(weights) 
    delta <- abs(currPos - oldPos)
    
    if(traj) trajectory <- rbind(trajectory, currPos)
  }
  
  list("estim" = currPos, "traj" = trajectory)
}