########
#' Parameters covariance matrix for "sml" objects
#'
#' @param llk     (numeric) vector containing the log-likelihood values.
#' @param parMat  (matrix) where the i-th row is the parameter vector corresponding to the i-th log-likelihood.
#' @param nreps   (integer) number of simulations used to tilt the negative Hessian (-H)
#'                toward positive definiteness. Used only if -H is not PD.
#' @param boot    (logical) relevant only if -H is not PD. 
#'                If TRUE hessians will be simulated by resampling parameters and likelihoods. 
#'                If FALSE hessians will be simulated the asymptotic distribution of the regression
#'                coefficients.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @export
#' 

.vcov.ml <- function(llk, parMat, nreps = 1000, boot = TRUE, ...)
{ 
  npar <- ncol(parMat)
  linPar <- 1 : (npar+1)
  
  # Non finite loglikelihoods are excluded
  good <- which( is.finite(llk) )
  parMat <- parMat[good, ]
  llk <-    llk[ good ]
  nval <- length(llk)
  
  # Creating model matrix and weights
  X <- .quadModMat(parMat)
  W <- exp(llk - max(llk))
  
  # Regressing estimated log-likelihoods of parameters
  fit <- lm(llk ~ -1 + X, weights = W)
  meanCoef <- coef(fit)[ -linPar ]
  ncoef <- length(meanCoef)
  
  # Extract hessian and its eigen-values
  hess <- - .extractHessian(meanCoef, npar)
  eig <- eigen(hess, only.values = TRUE)$values
  
  ###### Start of PD correction
  # If hessian is not PD, bootstrap or simulate using asymptotic covariance to get it PD (hopefully).
  # Our final estimate is the positive definite hessian with the lowest conditioning number.
  if( any(eig < 0) )
  {
    message("The estimated Hessian of the log-lik is not negative definite, I will try to bootstrap it.")
    
    # Resampled coefficients will be stored here by row
    coefMat <- matrix(NA, nreps, ncoef)
    
    # Using boostrap
    if(boot)
    {
      for(ii in 1:nreps)
      {
        index <- sample(1:nval, nval, replace = TRUE)
        
        tmpX <- X[index, ]
        tmpllk  <- llk[index]
        tmpW <- exp( tmpllk - max(tmpllk) )
        
        coefMat[ii, ] <- lm.wfit(x = tmpX, y = tmpllk, w = tmpW)$coefficients[ -linPar ]
      }
    } else {
      # Using asymptotic covariance
      coefMat <- .rmvn(nreps, mu = meanCoef, sigma = vcov(fit)[ -linPar, -linPar ]) 
    }
    
    # Extract simulated hessian and their conditioning number
    coefMat <- split(t(coefMat), rep(1:nreps, each = ncoef))
    hess <- lapply(coefMat, function(inMat) - .extractHessian(inMat, npar))
    cond <- sapply(hess, 
                   function(inHess){ 
                     eig <- eigen(inHess, only.values = TRUE)$values
                     if( all(eig > 0) ){ return( max(eig) / min(eig) ) } else { return( NA ) } 
                   })
    
    if( all(is.na(cond)) ) stop("Cannot get a positive definite negative Hessian.")
    
    # Use hessian with lowest conditioning number
    hess <- hess[[ which.min(cond) ]]
  } 
  ###### End of PD correction
  
  # Getting covariance, standard errors and confidence intervals
  covar <- .qrInverse(hess)
  rownames(covar) <- colnames(covar) <- colnames(parMat)
  
  return(covar) 
}

