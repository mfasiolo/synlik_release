#####
#' Empirical cumulant generating function
#' @description Calculates the empirical cumulant generating function (CGF) and its derivatives
#'               given a sample of n d-dimentional vectors
#'
#' @param lambda Point at which the CGF is evaluated (d-dimensional vector).
#' @param X (n by d) matrix containing the data.
#' @param kum1 Mean vector of the data.
#' @param kum2 Covariance matrix of the data.
#' @param grad If grad == 0 only the value of the CGF is returned, 
#'             if grad == 1 also its first derivative wrt lambda 
#'             and if grad == 2 also the second derivarive wrt lambda.
#' @param deriv If TRUE the gradient of the empitical CGF wrt y (and at y) is returned.
#'              Otherwise the values of the empirical CGF (and possibly of its derivatives wrt
#'              lambda) at lambda is returned.
#' @param mix Mixture of empirical and normal CGF to use (if 1 only empirical CGF is used).
#' @param addList = list of additional (optional) arguments: 
#'         \itemize{
#'         \item{ \code{invCOV} }{The inverse of kum2;}
#'         \item{ \code{y} }{The point at which the underlying empirical saddlepoint is evaluated;}
#'         \item{ \code{grad} }{The decay rate of the saddlepoint. See ?dsaddle for details;}
#'         }
#' @return If deriv == FALSE a list with entries:
#'         \itemize{
#'         \item{ \code{K} }{The value of the empirical CGF at lambda;}
#'         \item{ \code{dK} }{The value of the gradient empirical CGF wrt lambda at lambda;}
#'         \item{ \code{d2K} }{The value of the hessian of the empirical CGF wrt lambda at lambda;}
#'         }
#'         If deriv == TRUE the gradient of the empitical CGF wrt y (and at y) is returned. 
#'         This is used to calculate the gradient of the underlying empirical saddlepoint density
#'         at y. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon Wood.
#'
.ecgf <- cmpfun( function(lambda, X, kum1, kum2, grad = 0, mix = 1, deriv = FALSE,  addList = list(NaN) ) {
  ## X[i,j] is ith rep of jth variable. Evaluate observed KGF 
  ## and its derivs w.r.t. lambda, without overflow...
  
  if(grad > 2 && is.nan(addList[[1]])) stop("if you want the gradient you have to specify the additional arguments")
  
  if(!is.vector(lambda)) lambda <- as.vector(lambda) 
  if (!is.matrix(X)) X <- matrix(X, length(X), 1)
  n <- nrow(X)
  d <- ncol(X)
  lx <- drop( colSums( lambda*t(X) ) )
  alpha <- max(lx)       # constant for preventing overflow
  elx <- exp(lx - alpha) # exp(lambda'x_i - alpha) vector
  
  tmp.K <- log( sum(elx) / n )+alpha
  K <-  mix * ( tmp.K ) + (1 - mix) * ( crossprod(kum1, lambda) + 0.5 * crossprod(lambda, kum2%*%lambda) )
  ret <- list("K" = K)
  
  if(grad > 0)
  {
    tmp.dK <- colSums( elx*X ) / sum( elx )
    ret$dK <- mix * tmp.dK + ( 1-mix ) * ( kum1 + kum2%*%lambda ) 
  }
  
  if(grad > 1)
  {
    tmp.d2K <- crossprod(X, elx * X) / sum(elx) - tcrossprod(tmp.dK, tmp.dK)
    ret$d2K <-  as.matrix( mix * tmp.d2K  + (1-mix) * kum2, drop = FALSE) 
  }
  
  if(deriv)
  {
    dGam <- drop( (mix * addList$decay / d) * addList$invSigma %*% (kum1 - addList$y) ) 
    
    D <- diag( diag(ret$d2K) ^ -0.5, nrow = d, ncol = d)
    d2KQR <- qr( D %*% ret$d2K %*% D, tol = 0 )
    
    dLambda <- D %*% qr.solve(d2KQR, D %*% (diag(1, d) - tcrossprod(tmp.dK - kum1 - kum2%*%lambda, dGam)), tol = 0)
    
    d3K <- array(NA, c(d, d, d) )
    
    for(ff in 1:d)
    {
      A <- crossprod(X, (elx*X)*X[ , ff] ) / sum(elx)
      B <- -( tmp.d2K + tcrossprod(tmp.dK, tmp.dK) ) * tmp.dK[ff] 
      C <- - tcrossprod(tmp.d2K[ , ff], tmp.dK) - t( tcrossprod(tmp.d2K[ , ff], tmp.dK) )
      d3K[ , , ff] <- A + B + C 
    }  
    
    d2Kdy <- matrix(0, d, d)
    spaGrad <- numeric(d)
    for( ff in 1:d )
    {
      d2Kdy <- d2Kdy * 0 
      for(zz in 1:d)
      {
        d2Kdy <- d2Kdy + d3K[ , , zz] * dLambda[zz, ff] # Should it be [ff, zz] ??
      }
      
      spaGrad[ff] <- - mix * 0.5 * .Trace( D %*% qr.solve(d2KQR, D %*% d2Kdy, tol = 0) ) 
    }
    
    spaGrad <- spaGrad - 
      lambda + 
      dLambda %*% (ret$dK - addList$y) -  
      0.5 * .Trace( D %*% qr.solve(d2KQR, D %*% ( tmp.d2K - kum2 ), tol = 0) ) * dGam + 
      dGam * drop( tmp.K - crossprod(kum1, lambda) - 0.5 * crossprod(lambda, kum2%*%lambda) )
    
    return(spaGrad)
    
  } else return(ret)
  
})