
initCov.ml <- function(likFun, initPar, initCov, np, priorFun = NULL, 
                       multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                       constr = list(), verbose = FALSE, ...)
{
  # Including only parameters that were estimated
  varPar <- diag(initCov) > 0
  
  # Simulating parameter vectors
  simPar <- .paramsSimulator(theMean = initPar, covar = initCov, nsim = np, constr = constr)
  
  # Evaluating the likelihoods
  llk <- .funEval(parMat = simPar, fun = likFun, multicore = multicore, ncores = ncores, cluster = cluster, 
                  libraries = "synlik", ...) # NB "libraries" is being passed to .clusterExport()
  
  covar <- .vcov.ml(llk = llk, parMat = simPar[ , varPar], ...)
    
  return(covar)
  
}