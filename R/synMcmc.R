#######
### Method to estimate the parameters through MCMC
#######
#' MCMC parameter estimation for "synlik" objects.
#' 
#' @param object An object of class "synlik".
#' @param initPar Vector of initial parameters where the MCMC chain will start.
#' @param nIter Number of MCMC iterations.
#' @param burnIn Number of initial MCMC iterations which will be discarded.
#' @param priorFun Function that takes a vector of parameters as input and gives either 0 or 1 as output (uniform prior).
#' @param propCov Matrix representing the covariance matrix to be used to perturb the parameters at each step of the MCMC
#'                chain.
#' @param nsim  Number of simulated statistics at each step.
#' @param sameSeed If TRUE the synthetic likelihood will be evaluated at the current and proposed positions in the parameter
#'                 space using the same seed (thus doubling the computational effort). If FALSE the likelihood of the current
#'                 position won't be re-estimated.             
#' @param multicore  (logical) if TRUE the object@@simulator and object@@summaries functions will
#'                    be executed in parallel. That is the nsim simulations will be divided in multiple cores.
#' @param ncores  (integer) number of cores to use if multicore == TRUE.
#' @param cluster an object of class c("SOCKcluster", "cluster"). This allowes the user to pass her own cluster,
#'                which will be used if multicore == TRUE. The user has to remember to stop the cluster. 
#' @param ... additional arguments to be passed to synlikEval, see \code{\link{synlikEval}}.
#' @return An object of class "synMcmc".
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>                         
#' @export
#' 
synMcmc <- function(object, 
                    initPar, 
                    nIter, 
                    nsim,
                    propCov, 
                    burnIn = 0,
                    priorFun = function(param){ 1 },
                    targetRate = NULL,
                    recompute = FALSE,
                    multicore = !is.null(cluster),
                    cluster = NULL,
                    control = list(),
                    ...)
{
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
  
  # Reduce the object to "synlik" so that I avoid moving around all the additional slots of the "synMaxlik" class
  if( !class(object)[[1]] != "synlik" ) object <- as(object, "synlik")
  
  y <- object@data
  totalIter <- nIter + burnIn
  
  # Control list which will be used internally
  ctrl <- list( "ncores" = detectCores() - 1, 
                "theta" = 0.5,
                "adaptStart" = 0,
                "adaptStop" = totalIter,
                "saveFile" = NULL,
                "saveFreq" = 100,
                "verbFreq" = 500, 
                "verbose"  = FALSE )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  # Safety checks on ctrl
  if(ctrl$theta > 1 || ctrl$theta < 0.5) stop("control$theta should be between 0.5 and 1")
  if( !is.null(ctrl$saveFile) && !is.character(ctrl$saveFile) ) stop("\"ctrl$saveFile\" should be a character vector")
  stopifnot( ctrl$adaptStart <= ctrl$adaptStop, 
             ctrl$adaptStart >= 0,  
             ctrl$adaptStop <= totalIter) 
  
  # Check other arguments
  if(is.matrix(propCov) == FALSE) propCov <- diag(propCov)
  if(nrow(propCov) != ncol(propCov)) stop("propCov should be a square matrix")
  if(nrow(propCov) != length(initPar)) stop("nrow(propCov) != length(initPar)")
  
  # If a parameter has variance 0 in the proposal we save it's index in "fixPar"
  # we modidy the covariance and we save the initial covariance
  fixPar <- which( diag(propCov) == 0 )
  anyFix <- ( length(fixPar) > 0 )
  savedCov <- propCov
  if(anyFix) diag(propCov)[fixPar] <- 0.1

  cholFact <- t( chol( unname(propCov) ) )
  
  nPar <- length(initPar);
  currPar <- unname( initPar );
  propPar <- numeric(nPar);
  
  mcmcSample <- matrix(NA, nIter, nPar);
  LogLikChain <- numeric(nIter);
  
  currLogLik <- propLogLik <- tmpLik <- -10^99;
  currPrior <- priorFun(initPar, ...)
  
  accept <- 0
  
  unifVar <- runif(totalIter)
  
  if(multicore){
    tmp <- .clusterSetUp(cluster = cluster, ncores = ctrl$ncores, libraries = "synlik")
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
  }

  # Mcmc main loop
  storeIndex <- 1
  for(ii in 1:totalIter){
    
    # Propose new parameters
    pert <- rnorm(nPar)
    propPar <- currPar + as.vector( cholFact %*% pert )
    
    # Fix some parameters (if necessary) and check prior
    if( anyFix ) propPar[fixPar] <- currPar[fixPar]
    propPrior <- priorFun(propPar, ...)
    
    if( propPrior != -Inf )
    { 
      # Compute likelihood of proposed param
      propLogLik <- try( synlikEval(object, param = propPar, nsim = nsim, multicore = multicore, ncores = ctrl$ncores, cluster = cluster, ...)$logLik )
      if( !is.numeric(propLogLik) || !is.finite(propLogLik) ) propLogLik <- -Inf
      
      # (Optionally) recompute likelihood at old parameters
      if(recompute){ 
        tmpLik <- try( synlikEval(object, param = currPar, nsim = nsim, multicore = multicore, ncores = ctrl$ncores, cluster = cluster, ...)$logLik )
        if( is.numeric(tmpLik) && is.finite(tmpLik) ) currLogLik <- tmpLik
      }
      
      # Compute acceptance ratio
      likDiff <- propLogLik  - currLogLik + propPrior - currPrior
      alpha <- min(1, exp(likDiff))
      if( !is.finite(alpha) ) alpha <- ifelse( likDiff >= 0, 1, 0) 
      
      # Accept/Reject
      if ( likDiff > log(unifVar[ii]) ) {
        currPar <- propPar
        currPrior <- propPrior
        currLogLik <- propLogLik
        if(ii > burnIn) accept <- accept + 1
      }
      
    } else { alpha <- 0 }
    
    # Store iteration if iteration > burn-in
    if(ii > burnIn) {
      mcmcSample[storeIndex, ] <- currPar;
      LogLikChain[storeIndex] <- currLogLik;
      storeIndex <- storeIndex + 1
    }
    
    # (Optionally) adapt the proposal distribution, by updatint the transpose of its Cholesky factor
    if( !is.null(targetRate) && (ii >= ctrl$adaptStart) && (ii <= ctrl$adaptStop) )
    {
      cholFact <- .adaptChol(nPar = nPar, iter = ii, S = cholFact, U = pert, 
                             gamma = ctrl$theta, alpha = alpha, accRate = targetRate)
    }
    
    # (Optionally) save the object to file
    if( !is.null(ctrl$saveFile) && !(ii %% ctrl$saveFreq) ){ 
      save(file = ctrl$saveFile, 
           .synMcmc(object,
                    initPar = initPar,
                    nIter = nIter,
                    nsim = nsim, 
                    propCov = propCov,
                    burnIn = burnIn,
                    priorFun = priorFun,
                    targetRate = targetRate,
                    recompute = recompute,
                    multicore = multicore,
                    control = control,
                    
                    acceptRate = accept/nIter,
                    mcmcChain = mcmcSample,
                    LogLikChain = LogLikChain))
    }
    
    # (Optionally) print out intermediate results
    if( ctrl$verbose == TRUE && (ii > burnIn) && !(ii %% ctrl$verbFreq) )
    {
      tmp <- colMeans(mcmcSample[1:ii, ])
      names(tmp) <- names(object@param)
      cat(paste("Empirical posterior means at iteration", ii - burnIn, "\n"))
      print(tmp)
    }
    
  }
    
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  return( .synMcmc(object,
                   initPar = initPar,
                   nIter = nIter,
                   nsim = nsim, 
                   propCov = propCov,
                   burnIn = burnIn,
                   priorFun = priorFun,
                   targetRate = targetRate,
                   recompute = recompute,
                   multicore = multicore,
                   control = control,
                   
                   acceptRate = accept/nIter,
                   mcmcChain = mcmcSample,
                   LogLikChain = LogLikChain)  )
    
}

