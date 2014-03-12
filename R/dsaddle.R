#####
#' Empirical saddlepoint density 
#' @description Gives a pointwise evaluation of the empirical saddlepoint and optionally of
#'              its gradient at position y
#'
#' @param y Point at which the SPA is evaluated (d dimensional vector).
#' @param X n by d matrix containing the data.
#' @param tol Tolerance used to assess the convergence of the rootfinding routine used to fit
#'            the saddlepoint density. Default value is 1e-6.
#' @param decay Rate at which the SPA falls back on a normal density. Should be a positive number,
#'              by default set to 0.5.
#' @param deriv If TRUE also the gradient of the log-saddlepoint density is returned.
#' @param log If TRUE the log of the saddlepoint density is returned.
#' @return A list with entries:
#'         \itemize{
#'         \item{ \code{llk} }{The value of the empirical saddlepoint at y;}
#'         \item{ \code{mix} }{The mixture of normal-saddlepoint used (1 means only saddlepoint);}
#'         \item{ \code{grad} }{The gradient of the log-density at y (optional);}
#'         }
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon Wood.
#' @export
#'
dsaddle <- cmpfun(function(y, X, tol=1e-6, decay = 0.5, deriv = FALSE, log = FALSE) {
  ## X[i,j] is ith rep of jth variable; y is vector of variables.
  ## evaluate saddle point approximation based on empirical CGF
  
  if( !is.vector(y) ) y <- as.vector(y)
  
  m <- length(y)
  if(!is.matrix(X)){if(m > 1){stop("Error: simulated data must be entered in matrix form")
  }else{X <- matrix(X, length(X), 1)}}
  
  er <- .robCov(t(X), alpha2 = 4, beta2 = 1.25)
  # If there are some statistics with zero variance we remove them
  if( length(er$lowVar) )
  {
    y <- y[-er$lowVar]
    m <- m - length(er$lowVar)
    X <- X[ , -er$lowVar, drop = FALSE]
  }
  
  rss <- sum((er$E %*% (y - er$mY))^2)
  mix <- as.numeric( exp(-decay*rss / (2*m)) )
  
  # Weighting the statistics for a downweight outliers
  X <- er$weights * X
  
  # Initial guess of the root is the solution to the Gaussian case
  # the gain is one step less of Newton on average.
  lambda <- drop( crossprod(er$E, er$E %*% (y - er$mY)) )
  
  b <- .ecgf(lambda, X, kum1 = er$mY, kum2 = er$COV, mix = mix, grad = 2)
  
  kk <- 1
  while(TRUE) ## Newton loop to minimize K(lambda) - t(lambda)%*%y
  { 
    kk <- kk + 1
    
    D <- diag(diag(b$d2K)^-0.5, nrow = m, ncol = m)
    d.lambda <- -drop( D %*% qr.solve(D%*%b$d2K%*%D, D%*%(b$dK-y), tol = 0) )
    lambda1 <- lambda + d.lambda ## trial lambda
    
    b1 <- .ecgf(lambda1, X, kum1 = er$mY, kum2 = er$COV, mix = mix, grad = 2)
    if ( sum( abs(b1$d2K) ) == 0 ) return(NA) ## spa breakdown (possibly too tight)
    
    jj <- 1
    c1 <- 10^-4
    alpha <- 1
    rho <- 0.5
    ## Line seach checking Arminjo condition and step halving
    while( ( b1$K - crossprod(lambda1, y) ) > 
             ( b$K - crossprod(lambda, y) ) + c1 * alpha * crossprod(d.lambda, drop(b$dK-y)) && jj < 50)  
    {
      jj <- jj + 1
      alpha <- alpha * rho
      d.lambda <- d.lambda * alpha
      lambda1 <- lambda + d.lambda
      b1 <- .ecgf(lambda1, X, kum1 = er$mY, kum2 = er$COV, mix = mix, grad = 2)
    }
    
    ## now update lambda, K, dK, d2K
    lambda <- lambda1 
    b <- b1
    
    ## Converence test, +1 there just in case b$K = 0 at the solution
    if ( all( abs(b1$dK-y) < tol ) || kk > 100) break
    
  } ## end of Newton loop
  ## lambda is the SPA lambda...
  
  if(kk > 50 || jj > 20) warning(paste("The convergence of the saddlepoint root-finding is quite slow! Inner =", jj,"Outer=", kk))
  
  spa <- b$K - crossprod(lambda, y) - 0.5 * log( 2*pi ) * m - 0.5 * as.numeric( determinant(b$d2K, logarithm=TRUE)$modulus )
  if( !log ) spa <- exp(spa)
  
  outList <- list("llk" = drop(spa), "mix" = mix)
  
  if(deriv)
  {
    grad <- .ecgf(lambda = lambda, X = X, kum1 = er$mY, kum2 = er$COV, mix = mix, grad = 2, 
                 deriv = TRUE, addList = list("y" = y, "decay" = decay, "invSigma" = crossprod(er$E, er$E)) )
    outList$grad <- grad
  }
  
  return(outList)
  
})
