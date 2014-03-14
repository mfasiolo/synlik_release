#'
#' Estimated coverage of the synMaxlik optimization rountine
#'
#' @param object  ("synlik") object.
#' @param nrep    Number of time the whole simulation is repeated in order ot obtain the coverage.
#' @param niter   (integer) numer of iterations.
#' @param nsim    (integer) numer of simulations from the model at each step.
#' @param initCov  (matrix) initial covariance matrix used to simulate the paramters at each step.
#' @param initPar  (numeric) vector of initial values of the parameters.
#' @param addRegr  (logical) 
#'           if FALSE the statistics calculated by object@@summaries will be used (SL approach)
#'           if TRUE the simulated parameters will be regressed on the statistics and the 
#'           fitted values of the paramaters given the _observed_ statistics will be used as statistics
#'            (SL+ approach)
#' @param constr (named list) of 3 elements:
#'           [["indexes"]] = (numeric integers) indexes of the elements to check;
#'           [["upper"]]  = (numeric) upper bounds for the elements in "indexes";
#'           [["lower"]]  = (numeric) lower bounds for the elements in "indexes".
#' @param control Named list of control setting for the optimization routine.
#' @param multicore  (logical) if TRUE the object@@simulator and object@@summaries functions will
#'                    be executed in parallel. That is the nsim simulations will be divided in multiple cores.
#' @param ncores  (integer) number of cores to use if multicore == TRUE.
#' @param cluster an object of class c("SOCKcluster", "cluster"). This allowes the user to pass her own cluster,
#'                which will be used if multicore == TRUE. The user has to remember to stop the cluster. 
#' @param verbose  (logical) if TRUE lots of things will be printed.
#' @param ...  additional arguments to be passed to object@@simulator and object@@summaries.
#'             In general I would avoid using it and including in those two function everything they need.
#' @return object 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @export checkCoverage.synMaxlik

checkCoverage.synMaxlik <- function(object, nrep, niter, nsim, initCov, initPar = object@param, 
                                    addRegr = TRUE, constr = list(), control = list(),
                                    multicore = FALSE, 
                                    ncores = min(detectCores() - 1, nrep), 
                                    cluster = NULL, 
                                    verbose = FALSE,  ...)
{
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
    
  trueParam <- object@param
  
  # Function that simulates new data and estimates the parameters
  funToApply <- function(input, ...)
  {
    object@data <- simulate(object, nsim = 1, stats = FALSE, ...)
    object <- synMaxlik(object = object, 
                        niter = niter, 
                        nsim = nsim, 
                        initCov = initCov, 
                        initPar = initPar, 
                        addRegr = addRegr, 
                        constr = constr, 
                        control = control,
                        multicore = FALSE,
                        verbose = verbose,  ...)
    return( suppressWarnings( coef.synMaxlik(object)$par ) )
  }
  
  # Running synMaxlik optimization nrep times, possibly on multiple cores
  if(multicore){
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik", exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    
    # I put the environment of the function to apply to .Global to avoid that exports all the enviroment at every round of clusterApply.
    environment(funToApply) <- .GlobalEnv
    
    estimates <- clusterApply(cluster, seq.int(1, nrep), funToApply, ...)
    
    if(clusterCreated) cluster <- stopCluster(cluster)
  } else {
    
    estimates <- lapply(seq.int(1, nrep), funToApply, ...)
    
  }
  
  
  # Calculating the coverage
  tmp <- lapply(estimates, 
                function(estim)
                {
                  tmp <- sweep(estim, 1, trueParam)
                  inCover <- (tmp[ , 1] < 0) & (tmp[ , 3] > 0)
                  return( inCover )
                })
  
  cover <- numeric(length(trueParam))
  nGood <- 0
  for(ii in 1:nrep)
  {
    if( !any(is.na(tmp[[ii]])) ) 
      {
      cover <- cover + tmp[[ii]]
      nGood <- nGood + 1
      }
  }
  
  if(nGood < nrep) warning(paste(nrep-nGood, "out of", nrep, "confidence intervals couldn't be calculated because the Hessian is not positive definite."))
  
  cover <- matrix(cover/nGood * 100, 1, length(trueParam))
  rownames(cover) <- "Coverage %"
  colnames(cover) <- names(trueParam)
                              
  print(cover)
                            
  return( invisible( list("Coverage" = cover, "Estimates" = estimates) ) )
}
  
