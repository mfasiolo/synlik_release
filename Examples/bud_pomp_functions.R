budModCreator <- function (which) {
  if (missing(which)) {
    datasets <- c("food","para1","para2","tri")
    cat("available datasets:",sQuote(datasets),"\n")
    invisible(datasets)
  } else {
    which <- as.character(substitute(which))
      pomp(
        data=data.frame(
          time=seq(from=0,to=60,by=1),
          Qobs=NA,Nobs=NA,Sobs=NA
        ),
        time="time",
        t0=-1,
        params=switch(
          which,
          tri=c(
            alpha=0.5, sig.alpha=0.1, gam=50, lambda=22,
            sig.lambda=0.25, g=0.08, delta=10,
            a=1.7, sig.a=0.1, w=0.15, beta0=0, beta1=35, u=0.9,
            sigQobs=0.03, sigNobs=0.5, sigSobs=0.1,
            Q.0=0.96, N.0=0.02, S.0=0.22
          ),
          food=c(
            alpha=0.5, sig.alpha=0.1, gam=20, lambda=5,
            sig.lambda=0.25, g=0.02, delta=10,
            a=1, sig.a=0.1, w=0, beta0=0, beta1=35, u=0.9,
            sigQobs=0.03, sigNobs=0.5, sigSobs=0.1,
            Q.0=0.96, N.0=0.02, S.0=0.22
          ),
          para1=c(
            alpha=0.5, sig.alpha=0.1, gam=50, lambda=22,
            sig.lambda=0.25, g=0.08, delta=0.5,
            a=1.7, sig.a=0.1, w=0.15, beta0=0, beta1=35, u=0.9,
            sigQobs=0.03, sigNobs=0.5, sigSobs=0.1,
            Q.0=0.96, N.0=0.02, S.0=0.22
          ),
          para2=c(
            alpha=0.5, sig.alpha=0.1, gam=50, lambda=10,
            sig.lambda=5, g=0.08, delta=0.5,
            a=1.7, sig.a=1, w=0.15, beta0=0, beta1=35, u=0.9,
            sigQobs=0.03, sigNobs=0.5, sigSobs=0.1,
            Q.0=0.96, N.0=0.02, S.0=0.22
          ),
          stop("unrecognized dataset ",sQuote(which),call.=FALSE)
        ),
        rprocess=euler.sim(
          step.fun="budmoth_map",
          delta.t=1,
          PACKAGE="pompExamples"
        ),
        dprocess=onestep.dens(
          dens.fun="budmoth_density",
          PACKAGE="pompExamples"
        ),
        rmeasure="budmoth_rmeasure",
        dmeasure="budmoth_dmeasure",
        skeleton.type="map",
        skeleton="budmoth_skeleton",
        PACKAGE="pompExamples",
        paramnames=c(
          "alpha","sig.alpha","gam","lambda","sig.lambda",
          "g","delta","a","sig.a",
          "w","beta0","beta1","u",
          "sigQobs","sigNobs","sigSobs"
        ),
        statenames=c(
          "Alpha","Lambda","A","Q","N","S"
        ),
        obsnames=c("Qobs","Nobs","Sobs"),
        initializer=function (params, t0, ...) {
          x <- c(params[c("Q.0","N.0","S.0")],c(0,0,0))
          names(x) <- c("Q","N","S","Alpha","Lambda","A")
          x
        },
        logitvar = switch(
          which,
          tri = c("alpha","Q.0","S.0","u"),
          food = c("alpha","Q.0","S.0","u"),
          para1 = c("alpha","Q.0","S.0","u"),
          para2 = c("alpha","Q.0","S.0","u"),
          stop("unrecognized dataset ",sQuote(which),call.=FALSE)
        ),
        logvar = switch(
          which,
          tri = c(
            "sig.alpha","gam","lambda","sig.lambda",
            "g","delta","a","w","sig.a","beta1","sigQobs",
            "sigNobs", "sigSobs","N.0"),
          food = c(
            "sig.alpha","gam","lambda","sig.lambda",
            "g","delta","a","sig.a","beta1","sigQobs",
            "sigNobs", "sigSobs","N.0"),
          para1 = c(
            "sig.alpha","gam","lambda","sig.lambda",
            "g","delta","a","w","sig.a","beta1","sigQobs",
            "sigNobs", "sigSobs","N.0"),
          para2 = c(
            "sig.alpha","gam","lambda","sig.lambda",
            "g","delta","a","w","sig.a","beta1","sigQobs",
            "sigNobs", "sigSobs","N.0"),
          stop("unrecognized dataset ",sQuote(which),call.=FALSE)
        ),
        
        parameter.transform=function (params, logitvar,
                                      logvar, ...) {
          params[logitvar] <- plogis(params[logitvar])
          params[logvar] <- exp(params[logvar])
          params
        },
        parameter.inv.transform=function (params, logitvar,
                                          logvar, ...) {
          params[logitvar] <- qlogis(params[logitvar])
          params[logvar] <- log(params[logvar])
          params
        }
      )
  }
}



budWrap <- function(object, params, Np, trans = FALSE, ...)
{
  
  if( trans )
  {
    params[object@userdata$logitvar] <- logistic(params[object@userdata$logitvar])
    params[object@userdata$logvar] <-   exp(params[object@userdata$logvar])
  } 
  
  pfilter(object = object, params = params, Np = Np, ...)@loglik
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
budPMCMC <- function(object, 
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
  if(anyFix) propCov[fixPar, fixPar] <- 0.1
  
  cholFact <- t(chol(propCov))
  
  nPar <- length(initPar);
  currPar <- initPar;
  propPar <- numeric(nPar);
  
  mcmcSample <- matrix(NA, nIter, nPar);
  LogLikChain <- numeric(nIter);
  
  currLogLik <- propLogLik <- tmpLik <- -10^99;
  currPrior <- priorFun(initPar, ...)
  
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
    propPrior <- priorFun(propPar, ...)
    
    if( propPrior != -Inf )
    { 
      # Compute likelihood of proposed param
      propLogLik <- try( budWrap(object = object, params = propPar, Np = nsim, ...) )
      if( !is.numeric(propLogLik) || !is.finite(propLogLik) ) propLogLik <- -Inf
      
      # (Optionally) recompute likelihood at old parameters
      if(recompute){ 
        tmpLik <- try( budWrap(object = object, params = currPar, Np = nsim, ...) )
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
