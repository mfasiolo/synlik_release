#########
#### Fast computation of multivariate normal pdf
#########
#' Fast computation of the multivariate normal pdf
#'
#' @param X matrix n by d where each row is a d dimensional random vector. Alternatively
#'          X can be a d dimensional vector.
#' @param mu vector of length d, representing the mean the distribution.
#' @param sigma covariance matrix (d x d). Alternatively is can be the cholesky decomposition
#'              of the covariance. In that case isChol should be set to TRUE.
#' @param isChol boolean set to true is sigma is the cholesky decomposition
#'               of the covariance.
#' @param log boolean set to true the logarithm of the pdf is required.
#' @return a vector of length n where the i-the entry contains the pdf of the i-th random vector.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> 
#' @export dmvnFast

dmvnFast <- function(X, mu, sigma, log = FALSE, isChol = FALSE, verbose = TRUE){
  
  if( !is.matrix(X) ) X <- matrix(X, 1, length(X))
  
  # Checking if there are zeros on the diagonal
  fix <- which( diag(sigma) == 0 )
  anyFix <- ( length(fix) > 0 )
  
  if( anyFix )
  {
    if( isChol ) stop("The cholesky decomposition has zeros on the diagonal!")
    if( verbose ) warning("The covariance has zeros on the diagonal, those variables are not taken into account.")
    X <- X[ , -fix]
    mu <- mu[ -fix ]
    sigma <- sigma[-fix, -fix]
  }
  
  .Call( "dmvnCpp", 
         X_ = X, 
         mu_ = mu, 
         sigma_ = sigma, 
         log_ = log, 
         isChol_ = isChol, 
         PACKAGE = "synlik" )
}



###
# Test
###
# 
# library(mvtnorm)
# library(microbenchmark)
# 
# N <- 100
# d <- 5
# mu <- 1:d
# X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
# tmp <- matrix(rnorm(d^2), d, d)
# mcov <- tcrossprod(tmp, tmp)
# myChol <- chol(mcov)
# 
# 
# a <- cbind(
#       drop(dmvnFast(X, mu, mcov)),
#       drop(dmvnFast(X, mu, myChol, isChol = TRUE)),
#       dmvnorm(X, mu, mcov))
# a[ , 1] / a[, 3]
# a[ , 2] / a[, 3]
# 
# microbenchmark(dmvnFast(X, mu, 
#                          myChol, 
#                           isChol = TRUE), 
#                dmvnFast(X, mu, mcov), 
#                dmvnorm(X, mu, mcov))
