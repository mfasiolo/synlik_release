
# Initialization for the covariance of "sml" objects.

initCov.sml <- function(object, initpar, initcov, np, nsim, priorfun = NULL, 
                        multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                        constr = list(), verbose = FALSE, ...)
{
  if( is.null(names(initpar)) ) names(initpar) <- names(object@param)
  
  # Force evaluation of everything in the environment, so it will available to likfun on cluster
  if( multicore ) .forceEval(ALL = TRUE)
  
  # Function that will be used by sapply() or clusterApply to evaluate the likelihood
  likfun <- function(param, temper, ...)
  {
    synlikEval(object, param, nsim, multicore = FALSE, cluster = NULL, temper = temper, ...)$logLik
  }
  
  # Calling general maximum likelihood method
  covar <- initCov.ml(likfun = likfun, 
                      initpar = initpar, 
                      initcov = initcov, 
                      np = np, 
                      nsim = nsim, 
                      priorfun = priorfun,
                      multicore = multicore, 
                      ncores = ncores, 
                      cluster = cluster, 
                      constr = constr,
                      verbose = verbose,
                      ...)
  
  rownames(covar) <- colnames(covar) <- names( object@param[ diag(initcov) > 0 ] )
  
  return( covar )
}