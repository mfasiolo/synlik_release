
initCov.ml <- function(likfun, initpar, initcov, np, priorfun = NULL, 
                       multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                       constr = list(), verbose = FALSE, ...)
{
  # Including only parameters that were estimated
  varPar <- diag(initcov) > 0
  
  # Simulating parameter vectors
  simpar <- .paramsSimulator(theMean = initpar, covar = initcov, nsim = np, constr = constr)
  
  # Evaluating the likelihoods
  llk <- .funEval(parMat = simpar, fun = likfun, multicore = multicore, ncores = ncores, cluster = cluster, 
                  libraries = "synlik", ...) # NB "libraries" is being passed to .clusterExport()
  
  covar <- .vcov.ml(llk = llk, parMat = simpar[ , varPar], ...)
    
  return(covar)
  
}