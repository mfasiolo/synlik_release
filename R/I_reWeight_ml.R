#######
# Recycling
#######

.reWeight.ml <- function(storage, currMean, currCov, minW = 0.5, maxW = 1.0, verbose = FALSE)
{
  mixNam <- ls( storage )
  
  if( all( diag(currCov) != 0 ) )
  {
    sigma <- chol(currCov)
    isChol <- TRUE
  } else{
    sigma <- currCov
    isChol <- FALSE
  }
  
  # Go through all the components of the mixture currently saved
  for(iter in mixNam)
  {
    # Calculate new log-weights
    logW <- pmin( dmvnFast(X = storage[[iter]]$X, mu = currMean, sigma = sigma, 
                           log = F, isChol = isChol, verbose = FALSE) - storage[[iter]]$dens, log(maxW) )
    logW[ logW < log(minW) ] <- NA
    
    storage[[iter]]$logW <- logW + storage[[iter]]$llk + storage[[iter]]$logprior
  }
  
  return( storage )
}
