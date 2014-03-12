#########
### Method to continue MCMC estimation of a synMCMC object
#########

continue.synMcmc <- function(object, 
                             nIter = object@nIter,
                             nsim = object@nsim,
                             propCov = object@propCov, 
                             targetRate = object@targetRate,
                             recompute = object@recompute,
                             multicore = object@multicore,
                             cluster = NULL,
                             control = object@control,
                             ...)
{
  if(!is(object, "synMcmc")) stop("To use mcmc you need an object of class \"synMcmc\"")
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = object@control, outerCtrl = control)
  
  # For initPar and burnIn unless they have been specified by the user, we put
  # initPar to the final mcmc points in "object" and we don't do any more burn in.
  tmpObject <- synMcmc(object, 
                       initPar = as.vector(tail(object@mcmcChain, 1)),
                       nIter = nIter,
                       nsim = nsim,
                       propCov = propCov, 
                       burnIn = 0,
                       priorFun = object@priorFun,
                       targetRate = targetRate,
                       recompute = recompute,
                       multicore = multicore,
                       cluster = cluster,
                       control = ctrl,
                       ...)
  
  avgAcceptRate <- (object@acceptRate*object@nIter + tmpObject@acceptRate*tmpObject@nIter) / (object@nIter + tmpObject@nIter)
  
  return(.synMcmc(tmpObject,
                  
                  initPar = object@initPar,   # Resetting the values of these two param, so we don't lose information  
                  burnIn = object@burnIn,
                  
                  nIter = object@nIter + nIter,
                  
                  acceptRate = avgAcceptRate, 
                  mcmcChain = rbind(object@mcmcChain, tmpObject@mcmcChain),
                  LogLikChain = append(object@LogLikChain, tmpObject@LogLikChain)
  ))
}

setMethod("continue", 
          signature = signature(object = "synMcmc"), 
          definition = continue.synMcmc)