#######
### Method to simulate from the model
#######

# If nsim == 1 returns a "synlik" with the simulation in @data, if nsim > 1
# it retuns a matrix nsim simulations (which could be anything if stats == FALSE)
simulate.synlik  <- function(object,
                             nsim, 
                             seed = NULL,
                             param = object@param, 
                             stats = FALSE,
                             clean = FALSE,
                             verbose = TRUE,
                             ...)
{
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
  
  # Reduce the object to "synlik" so that I avoid moving around all the additional slots of the "synMaxlik" class
  if( !class(object)[[1]] != "synlik" ) object <- as(object, "synlik")
  
  if(is.null(seed) == FALSE) set.seed(seed)
  
  # I copy these function so I can mtrace() them
  simulator <- object@simulator
  summaries <- object@summaries
  extraArgs <- object@extraArgs
  
  simul <- simulator(param = param, nsim = nsim, extraArgs = extraArgs, ...)
  
  # Cleaning the stats from nans ONLY if simul is a matrix
  if( clean ) simul <- .clean(X = simul, verbose = verbose)$cleanX
  
  if( (stats == TRUE) && !is.null(summaries) ) simul <- summaries(x = simul, extraArgs = extraArgs, ...)
  
  # Cleaning the stats from NANs
  if( clean ) simul <- .clean(X = simul, verbose = verbose)$cleanX
  
  return( simul )
  
}

setMethod("simulate", 
          signature = signature(object = "synlik"), 
          definition = simulate.synlik)