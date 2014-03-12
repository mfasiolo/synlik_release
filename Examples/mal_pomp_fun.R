########################################################################################################################
################           POMP MODELS FOR       V_O_L_E_S      M_O_D_E_L
########################################################################################################################

#######
# Observed density evaluator
#######

malObsDens <- function (y, x, t, params, log, ...) {
  
  dnbinom(x = y["Y"], mu = exp(params["rho"]) * x["H_Cases"], size = 1 / exp(params["psiSquare"]), log = log)
  
}


#######
# Observed process simulator
#######

malObsSim <- function (x, t, params, ...) {
  
  X <- x["H_Cases"]
  
  return( drop( unname( rnbinom(n = length(X), mu = exp(params["rho"]) * X, size = 1 / exp(params["psiSquare"]) ) ) ) )
          
}


#######
# Pomp model for voles
#######

malariaPomp <- pomp(
  data=data.frame(
    time=1:240,
    Y=NA
  ),
  times="time",
  rmeasure = malObsSim,
  dmeasure = malObsDens,
  t0=0
)


#########
# Process simulator
#########


malSirPomCreator <- function(population, rainfall, nSteps, obsInterval, splineBasis)
{
  
  # Malaria process function where the step function has been enclosed
  malProc <- function (xstart, times, params, ...)
  {
    nrep <- ncol(xstart)      # number of realizations
    ntimes <- length(times)   # number of timepoints
    
    if(nrow(params) != 26) stop("We need 25 parameters")
    
    params <- t(params)
    
    ## x is the array of states to be returned: it must have rownames
    x <- array(0, dim = c(8, nrep, ntimes) )
    rownames(x) <- rownames(xstart)
    
    # If this is the first iteration we have to untrasform the initial values
    if(times[1] == 0)  ###########!!!!!!!!#############
{
  tmp <- logistic(xstart[1:5, , drop = FALSE])
  xstart[1:5, ] <- round( tmp / colSums(tmp) * population[1] )
  
  # Exponenting lam1 and lam2
  xstart[6:7, ] <- exp(xstart[6:7, , drop = FALSE])
}

## xnow holds the current state values
x[ , , 1] <- xstart
rownames(x) <- c("S1", "S2", "E", "I1", "I2", "LAM1", "LAM2", "H_Cases")

for (kk in seq.int(2, ntimes, by = 1)) {
  
  tmp <- try(
    .Call("malStepPompCpp",
          params_ = params[ , 1:10, drop = FALSE], 
          splineParams_ = params[ , 13:18, drop = FALSE],
          splineBasis_ = splineBasis,
          initStates_ = t(x[ , , kk - 1]),
          rainfall_ = rainfall,
          Pop_      = population[times[kk]], #population[kk-1],
          nSteps_   = nSteps, 
          nSimul_   = nrep, 
          obsInterval_ = obsInterval,
          currIteration_ = times[kk-1], # -1 because in C++ we start from zero 
          PACKAGE = "synlik")
  )
  
  x[ , , kk] <- t(tmp) 
}

return(x)
  } 

return(malProc)

}








#######
### Method to estimate the parameters through MCMC
#######
#' MCMC parameter estimation for "synlik" objects.
#' 
#' @param object An object of class "pomp".
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
malariaPMCMC <- function(object, 
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
  ctrl <- synlik:::.ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
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
  if(anyFix) 
    for(ii in fixPar) propCov[ii, ii] <- 0.1
  
  cholFact <- t(chol(propCov))
  
  nPar <- length(initPar);
  currPar <- initPar;
  propPar <- numeric(nPar);
  
  mcmcSample <- matrix(NA, nIter, nPar);
  LogLikChain <- numeric(nIter);
  
  currLogLik <- propLogLik <- tmpLik <- -10^99;
  currPrior <- priorFun(initPar)
  
  accept <- 0
  
  unifVar <- runif(totalIter)
  
  # Mcmc main loop
  storeIndex <- 1
  for(ii in 1:totalIter){
    
    # Propose new parameters
    pert <- rnorm(nPar)
    propPar <- currPar + as.vector( cholFact %*% pert )
    
    # Fix some parameters (if necessary) and check prior
    if( anyFix ) propPar[fixPar] <- currPar[fixPar]
    propPrior <- priorFun(propPar)
    
    # Fix some parameters (if necessary) and check prior
    if( anyFix ) propPar[fixPar] <- currPar[fixPar]
    propPrior <- priorFun(propPar)
    
    if( propPrior != -Inf )
    { 
      # Compute likelihood of proposed param
      propLogLik <- try( pfilter(object = object, params = propPar, Np = nsim, ...)@loglik )
      if( !is.numeric(propLogLik) || !is.finite(propLogLik) ) propLogLik <- -Inf
      
      # (Optionally) recompute likelihood at old parameters
      if(recompute){ 
        tmpLik <- try( pfilter(object = object, params = currPar, Np = nsim, ...)@loglik )
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
      cholFact <- synlik:::.adaptChol(nPar = nPar, iter = ii, S = cholFact, U = pert, 
                                      gamma = ctrl$theta, alpha = alpha, accRate = targetRate)
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
  
  return(list("mcmcSample" = mcmcSample, "LogLikChain" = LogLikChain, "AccRate" = accept/totalIter));
  
}


