
sml <- function(object, initpar, initcov, np, nsim, niter, alpha = 0.95,
                priorfun = NULL, temper = NULL, recycle = FALSE,
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
  tmp <- ml(likfun = likfun, 
            initpar = initpar, 
            initcov = initcov, 
            np = np, 
            niter = niter,
            priorfun = priorfun,
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
               initpar = initpar, 
               initcov = initcov, 
               np = np, 
               nsim = nsim, 
               niter = niter,
               priorfun = priorfun,
               alpha = alpha,
               temper = temper,
               recycle = recycle,
               multicore = multicore, 
               ncores = ncores, 
               constr = constr,
               
               estim = tmp$estim,
               simloglik = tmp$simloglik,
               simlogprior = tmp$simlogprior,
               simpar = tmp$simpar)  )
}









########
# Adaptation
########
#     if( adapt )
#     {
#       propUpdate <- matrix(0, length(parMean), length(parMean))
#       for(kk in 1:np)
#       {
#         propUpdate <- propUpdate + w[kk] * tcrossprod(simpar[kk, ] - parMean, simpar[kk, ] - parMean)
#       }
#       
#       ESS <- 1 / sum(w ^ 2)
#       beta <- min(ESS / length(parMean), 0.25)
#       
#       parCov <- alpha * ( (1 - beta) * parCov + beta * propUpdate )
#     }