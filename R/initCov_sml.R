
# Initialization for the covariance of "sml" objects.

initCov.sml <- function(object, initPar, initCov, np, nsim, priorFun = NULL, 
                        multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                        constr = list(), verbose = FALSE, ...)
{
  if( is.null(names(initPar)) ) names(initPar) <- names(object@param)
  
  # Force evaluation of everything in the environment, so it will available to likfun on cluster
  if( multicore ) .forceEval(ALL = TRUE)
  
  # Function that will be used by sapply() or clusterApply to evaluate the likelihood
  likFun <- function(param, temper, ...)
  {
    slik(object, param, nsim, multicore = FALSE, cluster = NULL, temper = temper, ...)
  }
  
  # Calling general maximum likelihood method
  covar <- initCov.ml(likFun = likFun, 
                      initPar = initPar, 
                      initCov = initCov, 
                      np = np, 
                      nsim = nsim, 
                      priorFun = priorFun,
                      multicore = multicore, 
                      ncores = ncores, 
                      cluster = cluster, 
                      constr = constr,
                      verbose = verbose,
                      ...)
  
  rownames(covar) <- colnames(covar) <- names( object@param[ diag(initCov) > 0 ] )
  
  return( covar )
}