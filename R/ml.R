
# Maximimum likelihood

ml <- function(likfun, initpar, initcov, np, niter, priorfun = NULL, alpha = 0.95, temper = rep(1, niter), recycle = FALSE,
               multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
               constr = list(), verbose = FALSE, ...)
{  
  parCov <- initcov
  npar <- length(initpar)
  parMean <- initpar
  
  stopifnot( is.vector(temper), length(temper) == niter )
  
  outEstim <- matrix(NA, niter, npar)
  
  # Setting up an environment where recycled likelihoods, parameters and densities will be stored
  if( recycle ) {
    storage <- new.env( parent = emptyenv() )
  } else {
    storage <- NULL
  }
  
  # Set up the cluster
  if(multicore)
  {
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik") 
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
  }
  
  # Main loop 
  for(ii in 1:niter)
  { 
    # Simulating parameter vectors
    simpar <- .paramsSimulator(theMean = parMean, covar = parCov, nsim = np, constr = constr)
    
    # Evaluating the likelihoods
    llk <- .funEval(parMat = simpar, fun = likfun, multicore = multicore, cluster = cluster, ...)
    
    # Evaluate log-prior of each simulated set of parameters
    logprior <- if( !is.null(priorfun) ){ apply(simpar, 1, priorfun, ...) } else { rep(0, np) }
    
    # Adding the latest component to the mixture
    storage[[ as.character(ii) ]] <- list("X" = simpar, 
                                          "llk" = llk,
                                          "logprior" = logprior,
                                          "logW"   = rep(0, np),
                                          "dens" = dmvnFast(X = simpar, mu = parMean, sigma = parCov, log = F, verbose = FALSE)
    )
    
    w <- llk + logprior
    
    # Recycling old likelihood estimates
    if( recycle ){
      storage <- .reWeight.ml(storage, currMean = parMean, currCov = parCov, verbose = FALSE)
      
      w <- do.call("c", lapply(storage, "[[", "logW") )
      simpar <- do.call("rbind", lapply(storage, "[[", "X"))
    }
    
    # Discarding parameter with bad weights
    good <- !is.na(w)
    w <- w[ good ] 
    simpar <- simpar[good, ]
    
    # Tempering, exponential trick and normalize the weights
    w <- exp( temper[ii] * (w - max(w)) )
    w <- w / sum(w)
    
    # Update the parameters. If all the weights are zero we don't move.
    if( all(w == 0) ){
      warning("All weights equal to zero! Maybe use some tempering?")
      outEstim[ii, ] <- parMean 
    } else{
      parMean <- outEstim[ii, ] <- parMean + colSums( t(t(simpar) - parMean) * w  )      
    }
    
    parCov <- alpha * parCov
  }
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  # Extract all estimated likelihood and simulated parameters 
  timeOrder <- as.character(sort(as.numeric(ls(storage))))
  outloglik <- do.call("c", lapply(storage, "[[", "llk")[ timeOrder ] )
  outlogprior <- do.call("c", lapply(storage, "[[", "logprior")[ timeOrder ] )
  outPar <- do.call("rbind", lapply(storage, "[[", "X")[ timeOrder ] )
  
  colnames(outEstim) <- colnames(outPar) <- names(initpar)
  
  return( list("estim" = outEstim, "simpar" =  outPar, "simloglik" = outloglik, "simlogprior" = outlogprior  ) )
}
