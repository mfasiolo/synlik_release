##########################
#' Creates a summaries function using ridge regression or lasso
#' 
#' @description Fearnhead and Prangle (2012) proposed to use posterior
#'              expected values of the parameters as summary statistics. Here this goal is approximated
#'              by using fitted values of the parameters (obtained by ridge regression) as summary
#'              statistics. 
#' 
#' @param object An object of class "synlik".
#' @param nsim Number of summary statistics to be simulated.
#' @param mu   mean around which the parameters will be simulated.
#' @param sigma covariance matrix used to simulate the parameters.
#' @param ... even more additional arguments to be passed to object@@simulator and object@@summaries.
#' @return A function that transforms simulated data in summary statistics. This function can by used as
#'         object@@summaries function in an object of class "synlik".
#'            
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.                         
#' @export shrinkStat         

shrinkStat <- function(object, nsim, mu, sigma, type = "ridge", constr = list(), clean = TRUE, verbose = TRUE, ...) 
{
  
  regrCoef <- unname( shrinkCoef(object = object, nsim = nsim, mu = mu, sigma = sigma, type = type, constr = constr, clean = clean, verbose = verbose, ...)$regrCoef )
  
  origStats <- force( object@summaries ) # Forcing evaluation just in case
  
  # Create new summaries function
  shrinkSummaries <- function(x, extraArgs, ...)
  {
    simulStats <- origStats(x = x, extraArgs = extraArgs, ...)
    nsim <- nrow(simulStats)
    
    # Cleaning the stats from NANs
    clean <- .cleanStats(simulStats)
    if(clean$nBanned > 0) simulStats <- clean$cleanX
    if(nrow(simulStats) < nsim / 3) warning(paste(nsim - nrow(simulStats), "out of", nsim, "statistics vectors", "contain NAs and will not be used"))
        
    return( cbind(1, simulStats) %*% t(regrCoef) ) 
    
  }
  
  return(shrinkSummaries)
}


