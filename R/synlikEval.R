######################
###### Evaluates the likelihood function
######################
#' Evaluate the synthetic likelihood.
#' 
#' @param object An object of class "synlik".
#' @param param Vector of parameters at which the synthetic likelihood will be evaluated.
#' @param nsim  Number of simulation from the model.
#' @param saddle If FALSE the distribution of the summary statistics will be approximated by a normal, if TRUE an
#'               Empirical Saddlepoint approximation will be used.
#' @param decay Useful only if \code{saddle} == TRUE. It is the rate at which the saddlepoint density falls back to the 
#'              normal density as the distance between observed and simulated statistics increases.               
#' @param multicore  (logical) if TRUE the object@@simulator and object@@summaries functions will
#'                    be executed in parallel. That is the nsim simulations will be divided in multiple cores.
#' @param ncores  (integer) number of cores to use if multicore == TRUE.
#' @param cluster an object of class c("SOCKcluster", "cluster"). This allowes the user to pass her own cluster,
#'                which will be used if multicore == TRUE. The user has to remember to stop the cluster. 
#' @param ... additional arguments to be passed to object@@simulator and object@@summaries.
#'            In general I would avoid using it and including in those two function everything they need.
#' @return A list where "logLik" is the log of the estimated synthetic likelihood and "mix" is the share of saddlepoint used
#'         (0 means no saddlepoint (only normal) and 1 means only saddlepoint (no normal) ).
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>                         
#' @export
#' 
synlikEval <- function(object, param, nsim, saddle = FALSE, decay = 0.5, multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, ...) 
{
  
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
  
  # Reduce the object to "synlik" so that I avoid moving around all the additional slots of the "synMaxlik" class
  if( !class(object)[[1]] != "synlik" ) object <- as(object, "synlik")
  
  if( !is.vector(param) ) stop("param should be a numeric vector.")
    
  # Simulating from model
  if(multicore)
  {
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
       
    coresSchedule <- c( rep(floor(nsim / ncores), ncores - 1), floor(nsim / ncores) + nsim %% ncores)
    
    simulData <- clusterApply(
      cluster, 
      coresSchedule, 
      function(input, ...)
      {
        simulate.synlik(object, param = param, nsim = input, stats = TRUE, clean = TRUE, verbose = FALSE, ...)
      } 
      , ...
    )
    
    simulData <- do.call(rbind, simulData)
    
    if(clusterCreated) cluster <- stopCluster(cluster)
    
  } else {
    simulData <- simulate.synlik(object, param = param, nsim = nsim, stats = TRUE, clean = TRUE, verbose = FALSE, ...) 
  }
  
  if(nrow(simulData) < nsim / 3) warning(paste(nsim - nrow(simulData), "out of", nsim, "statistics vectors", "contain NAs and will not be used"))
  
  # Transforming the observation into stats
  summaries <- object@summaries
  obsStats <- if( !is.null(summaries) ) summaries(x = object@data, extraArgs = object@extraArgs, ...) else object@data
    
  # Calculating log-likelihood
  return( .empDens(y = obsStats, X = simulData, saddle = saddle, decay = decay, tol = 1e-6, log = TRUE) )
  
}





