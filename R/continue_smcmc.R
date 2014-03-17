#########
### Method to continue MCMC estimation of a synMCMC object
#########

##############################################################
#' Continuing estimation of an \code{smcmc} object.
#'
#' Continues MCMC estimation of an object of class \code{smcmc}. All input parameters are defaulted to the corresponding
#' slots in the input object, with the exception of cluster. The chain restarts were it ended, burn-in is set to zero, the
#' same prior (if any) is used. 
#'
#' @param object An object representing the results of an estimation procedure which we wish to continue.
#'               For example it might represents and MCMC chain.
#' @param niter  Number of additional MCMC iterations requested. The output chain will have length object@niter + niter.
#' @param nsim  see \code{\link{smcmc-class}}.
#' @param propCov see \code{\link{smcmc-class}}.
#' @param targetRate see \code{\link{smcmc-class}}.
#' @param recompute see \code{\link{smcmc-class}}.             
#' @param multicore  see \code{\link{smcmc-class}}.
#' @param ncores   see \code{\link{smcmc-class}}.
#' @param cluster see \code{\link{smcmc}}. 
#' @param control see \code{\link{smcmc-class}}. 
#' @param ... additional arguments to be passed to \code{slik()}, see \code{\link{slik}}.
#'
#' @return An object of the same class as \code{object}, where the results of the estimation have been updated.
#' 
#' @seealso \code{\link{smcmc-class}}, \code{\link{smcmc}}, \code{\link{continue}}.
#' 
#' @aliases continue,smcmc,missing-method
#' @examples
#' # For an example see help("smcmc-class").
#' @rdname continue-smcmc

continue.smcmc <- function(object, 
                           niter = object@niter,
                           nsim = object@nsim,
                           propCov = object@propCov, 
                           targetRate = object@targetRate,
                           recompute = object@recompute,
                           multicore = object@multicore,
                           ncores = object@ncores,
                           cluster = NULL,
                           control = object@control,
                           ...)
{
  if(!is(object, "smcmc")) stop("To use mcmc you need an object of class \"smcmc\"")
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = object@control, outerCtrl = control)
  
  # For initPar and burn unless they have been specified by the user, we put
  # initPar to the final mcmc points in "object" and we don't do any more burn in.
  tmpObject <- smcmc(object = object, 
                     initPar = drop( tail(object@chains, 1) ),
                     niter = niter,
                     nsim = nsim,
                     propCov = propCov, 
                     burn = 0,
                     priorFun = object@priorFun,
                     targetRate = targetRate,
                     recompute = recompute,
                     multicore = multicore,
                     ncores = ncores,
                     cluster = cluster,
                     control = ctrl,
                     ...)
  
  avgAcceptRate <- (object@accRate*object@niter + tmpObject@accRate*tmpObject@niter) / (object@niter + tmpObject@niter)
  
  return(new(   "smcmc",
                tmpObject,
                
                initPar = object@initPar,   # Resetting the values of these two param, so we don't lose information  
                burn = as.integer(object@burn),
                
                niter = as.integer(object@niter + niter),
                
                accRate = avgAcceptRate, 
                chains = rbind(object@chains, tmpObject@chains),
                llkChain = append(object@llkChain, tmpObject@llkChain)
  ))
}

setMethod("continue", 
          signature = signature(object = "smcmc"), 
          definition = continue.smcmc)