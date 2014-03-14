
sml <- function(object, initPar, initCov, np, nsim, niter, alpha = 0.95,
                priorFun = NULL, temper = NULL, recycle = FALSE,
                multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                constr = list(), verbose = FALSE, ...)
{
  if( is.null(names(initPar)) ) names(initPar) <- names(object@param)
  
  message("remember to check why the recycling isn't working!!")
  
  # Force evaluation of everything in the environment, so it will available to likfun on cluster
  if( multicore ) .forceEval(ALL = TRUE)
  
  # Function that will be used by sapply() or clusterApply to evaluate the likelihood
  likFun <- function(param, temper, ...)
  {
    slik(object, param, nsim, multicore = FALSE, cluster = NULL, temper = temper, ...)
  }
  
  # Calling general maximum likelihood method
  tmp <- ml(likFun = likFun, 
            initPar = initPar, 
            initCov = initCov, 
            np = np, 
            niter = niter,
            priorFun = priorFun,
            alpha = alpha,
            temper = temper,
            recycle = recycle,
            multicore = multicore, 
            ncores = ncores, 
            cluster = cluster, 
            constr = constr,
            verbose = verbose,
            ...)
  
  return( .sml(object, 
               initPar = initPar, 
               initCov = initCov, 
               np = np, 
               nsim = nsim, 
               niter = niter,
               priorFun = priorFun,
               alpha = alpha,
               temper = temper,
               recycle = recycle,
               multicore = multicore, 
               ncores = ncores, 
               constr = constr,
               
               estim = tmp$estim,
               simLogLik = tmp$simLogLik,
               simLogPrior = tmp$simLogPrior,
               simPar = tmp$simPar)  )
}









########
# Adaptation
########
#     if( adapt )
#     {
#       propUpdate <- matrix(0, length(parMean), length(parMean))
#       for(kk in 1:np)
#       {
#         propUpdate <- propUpdate + w[kk] * tcrossprod(simPar[kk, ] - parMean, simPar[kk, ] - parMean)
#       }
#       
#       ESS <- 1 / sum(w ^ 2)
#       beta <- min(ESS / length(parMean), 0.25)
#       
#       parCov <- alpha * ( (1 - beta) * parCov + beta * propUpdate )
#     }