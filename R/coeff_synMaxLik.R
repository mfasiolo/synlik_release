##########################
#' Get the final parameters' estimates from a synMaxLik object
#' 
#' @description This function extract the final estimate for each parameter from object@@resultPar 
#'              and calculates confidence intervals (if the estimated Hessian is positive definite).
#' 
#' @param object An object of class "synMaxlik".
#' @param Lag Number of values over which the Hessian is averaged to give it's final estimate. Should 
#'            be a number between 1 and object@@nIter.
#' @param tol Tolerance used to determine whether the estimated Hessian is positive definite (PD). If 
#'            smallest_eigen > largest_eigen * tol than the Hessian is considered to be PD.
#' @return A matrix where each row contains a 95% confidence interval for each parameter. If the estimated Hessian
#'         is not PD then the boundaries of the CI will be set to NA.       
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>                       
#' #@export coef.synMaxlik
#' #@method coef.synMaxlik synMaxlik
#' #@S3method coef.synMaxlik synMaxlik

coef.synMaxlik <- function(object, Lag = 10, tol = 1e-12, ...)
{
  if(!is(object, "synMaxlik")) stop("object has to be of class \"synMaxlik\" ")
  
  nIter <- nrow(object@resultPar)
  if(nIter == 0) stop("There are no stored iterations in the object (object@resultPar is empty)")
  
  parEstim <- object@resultPar
  nPar <- length(object@param)
  
  # The hessian is an average of the last "Lag" hessians
  Lag <- min(Lag, nIter)
  hess <- object@resultHess[ (nIter - Lag + 1):nIter ]
  hess <- Reduce("+", hess) / Lag
  
  # If the hessian is not positive definite we remove some bad parameters until it is
  tmp <- .cleanHessian(hess, tol = tol)
  COV <- qr.solve(tmp$hessian, tol = 0)
  badParam <- tmp$badParam
  
  if( length(badParam) > 0 ) 
    warning(paste("The estimated Hessian is not positive definite: cannot calculate confidence intervals for",  length(badParam), "parameters."))
    
  iGood <- 1
  output <- matrix(NA, nPar, 3)
  for(iPar in 1:nPar)
  {
    parVals <- parEstim[ , iPar]
    finalPar <- tail(parVals, 1)
    
    # Calculate confidence intervals the parameter has not been removed from the Hessian
    if( !(iPar %in% badParam) ){
      stdErr <- sqrt( COV[iGood, iGood] )
      confInt <- c(finalPar - 1.96 * stdErr, finalPar + 1.96 * stdErr)
      iGood <- iGood + 1
    }
    
    rownames(output) <- names(object@param)
    colnames(output) <- c("Lower C.I.", "Estimate", "Upper C.I.")
    
    output[iPar, 2] <- finalPar
    output[iPar, 1] <- ifelse( !(iPar %in% badParam), confInt[1], NA)
    output[iPar, 3] <- ifelse( !(iPar %in% badParam), confInt[2], NA)
    
  }
  
  if( length(badParam) )
  {
    dimnames(COV) <- list(rownames(output)[-badParam], rownames(output)[-badParam])
  } else{
    dimnames(COV) <- list(rownames(output), rownames(output))
  }
  
  return( list("par" = output, "cov" = COV) )
}


setMethod("coef", 
          signature = signature(object = "synMaxlik"), 
          definition = coef.synMaxlik)
