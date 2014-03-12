#####################################################################################################################

########################################## CLASS DEFINITIONS AND VALIDITY CHECKS ####################################

#####################################################################################################################

setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("ANYOrNULL", c("ANY", "NULL"))
setClassUnion("functionOrNULL", c("function", "NULL"))

##################################
######### synlik: the base class
##################################

### Validity check

.check.synlik <- function(object)
{
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
  
  errors <- character()
  
  if(length(object@param) == 0) errors <- c(errors, "length(param) == 0")
  if(is.null(names(object@param)) || any("" %in% names(object@param)) )
    errors <- c(errors, "param has to be a named vector")
  
  simulArgs <- names(as.list(args(object@simulator)))
  if( length(simulArgs) < 5 || !identical(simulArgs[1:3], c("param", "nsim", "extraArgs")) || simulArgs[length(simulArgs) - 1] != "...") 
    stop("The first 3 arguments of the simulator should be \"param\", \"nsim\" and \"extraArgs\" (in that order) and the last should be \"...\"")
  
  if( !is.null(object@summaries) )
  {
    statsArgs <- names(as.list(args(object@summaries)))
    if( length(statsArgs) < 4 || (statsArgs[1] != "x") || statsArgs[length(statsArgs) - 1] != "...") 
      stop("The first 2 argument of the \"summaries\" function should be \"x\" and \"extraArgs\" (in that order) and the last should be \"...\"")
  }
  
  if(length(errors) == 0) TRUE else errors
}


### Class Definition
#' \code{synlik-class}
#'
#' \describe{
#'    \item{param}{Vector of parameters used by object@@simulator (\code{numeric}).}
#'    \item{simulator}{Function that simulates from the model (\code{function}).}
#'    \item{summaries}{Function that transforms simulated data into summary statistics (\code{function}).}
#'    \item{data}{Vector containing the observed data (\code{numeric}).}
#'    \item{extraArgs}{List containing all the extra arguments to be passed to object@@simulator and object@@summaries (\code{list}).}
#'  }
#' @name synlik-class
#' @rdname synlik-class
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @exportClass synlik
.synlik <- setClass( "synlik",
                     representation( param = "numeric",
                                     simulator = "function",
                                     summaries = "functionOrNULL",
                                     statTrans = "functionOrNULL",
                                     data = "ANY",
                                     extraArgs = "list"
                     ),
                       prototype = prototype(
                       param = numeric(),
                       simulator = function() NULL,
                       summaries = NULL,
                       statTrans = NULL,
                       data = NULL, 
                       extraArgs = list()
                       ),
                     
                     validity = .check.synlik
)



##################################
######### synMcmc: a synlik object after MCMC
##################################

### Validity check

.check.synMcmc <- function(object)
{
  if(!is(object, "synMcmc")) stop("object has to be of class \"synMcmc\" ")

  errors <- character()
  
  if(length(object@initPar) == 0) errors <- c(errors, "length(initPar) should be > 0")
  
  if(length(errors) == 0) TRUE else errors
}



### Class Definition
#' \code{synMcmc-class}
#'
#' \describe{
#'    \item{initPar}{Vector of initial parameters where the MCMC chain will start (\code{numeric}).}
#'    \item{nIter}{Number of MCMC iterations (\code{numeric}).}
#'    \item{burnIn}{Number of initial MCMC iterations which are discarded (\code{numeric}).}
#'    \item{priorFun}{Function that takes a vector of parameters as input and gives 
#'                    either 0 or 1 as output (uniform prior) (\code{function}).}
#'    \item{propCov}{Matrix representing the covariance matrix to be used to perturb the 
#'                   parameters at each step of the MCMC chain (\code{matrix}).}
#'    \item{nsim}{Number of simulations from the simulator at each step of the MCMC algorithm (\code{numeric}).}
#'    \item{sameSeed}{If TRUE the synthetic likelihood will be evaluated at the current and proposed positions in the parameter
#'                    space using the same seed (thus doubling the computational effort). If FALSE the likelihood of the current
#'                    position won't be re-estimated (\code{logical}).}
#'    \item{multicore}{If TRUE the object@@simulator and object@@summaries functions will
#'                     be executed in parallel. That is the nsim simulations will be divided in multiple cores (\code{logical}).}
#'    \item{ncores}{Number of cores to use if multicore == TRUE (\code{numeric}).}
#'    \item{acceptRate}{Acceptance rate of the MCMC chain, between 0 and 1 (\code{numeric}).}
#'    \item{mcmcChain}{Matrix of size nIter by length(initPar) where the i-th row contains the position of the MCMC algorithm
#'                      in the parameter space at the i-th (\code{matrix}).}
#'    \item{LoglikChain}{Vector of nIter elements where the i-th element is contains the estimate of the 
#'                       synthetic likelihood at the i-th iteration (\code{numeric}).}
#'  }
#'  
#' @name synMcmc-class
#' @rdname synMcmc-class
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @exportClass synMcmc
#'
.synMcmc <-setClass("synMcmc",
                    representation( initPar = "numeric",
                                    nIter = "numeric",
                                    nsim = "numeric",
                                    propCov = "matrix",
                                    burnIn = "numeric",
                                    priorFun = "function",
                                    
                                    targetRate = "numericOrNULL",
                                    recompute = "logical",
                                    multicore = "logical",
                                    control = "list",
                                    
                                    acceptRate = "numeric",
                                    mcmcChain = "matrix",
                                    LogLikChain = "numeric"
                    ),
                    prototype = prototype(initPar = numeric(),
                                          nIter = numeric(),
                                          nsim = numeric(), 
                                          propCov = matrix( , 0, 0),
                                          burnIn = 0,
                                          priorFun = function(param) 1,
                                          
                                          targetRate = 0.5,
                                          recompute = FALSE,
                                          multicore = FALSE,
                                          control = list(),
                                          
                                          acceptRate = numeric(),
                                          mcmcChain = matrix( , 0, 0),
                                          LogLikChain = numeric()),
                    contains = "synlik",
                    validity = .check.synMcmc
)



##################################
######### synMaxlik: a synlik object after stochOptim
##################################

### Validity check

.check.synMaxlik <- function(object)
{
  if(!is(object, "synMaxlik")) stop("object has to be of class \"synMaxlik\" ")
  
  errors <- character()
  
  if(length(object@initPar) == 0) errors <- c(errors, "length(initPar) should be > 0")
  
  if(length(errors) == 0) TRUE else errors
}

### Class Definition
#' \code{synMaxlik-class}
#'
#' \describe{
#'    \item{initPar}{Vector of initial parameters from which the optimization starts (\code{numeric}).}
#'    \item{nIter}{Number of iterations of the optimization routine (\code{numeric}).}
#'    \item{nsim}{Number of simulations from the simulator at each step of the optimization routine (\code{numeric}).}
#'    \item{initCov}{Covariance matrix used to simulate the parameter at each step of the optimization routine (\code{matrix}).}
#'    \item{addRegr}{If FALSE the statistics calculated by object@@summaries will be used (SL approach). If TRUE the simulated 
#'                   parameters will be regressed on the statistics and the fitted values of the paramaters given the _observed_
#'                   statistics will be used as statistics (this is called SL+ approach) (\code{logical}).}
#'    \item{constr}{Named list of contraints on the parameters, it has 3 elements: 
#'                  [["indexes"]] = (numeric integers) indexes of the elements to check;
#'                  [["upper"]]  = (numeric) upper bounds for the elements in "indexes";
#'                  [["lower"]]  = (numeric) lower bounds for the elements in "indexes" (\code{list}).}
#'    \item{control}{Named list of control setting for the optimization routine (\code{list}).}
#'    \item{multicore}{If TRUE the object@@simulator and object@@summaries functions will
#'                     be executed in parallel. That is the nsim simulations will be divided in multiple cores (\code{logical}).}
#'    \item{ncores}{Number of cores to use if multicore == TRUE (\code{numeric}).}
#'    \item{resultPar}{Matrix typically nIter by length(initPar) where the i-th row contains the estimate of the 
#'                     parameters at the i-th iteration (\code{matrix}).}
#'    \item{resultGrad}{Matrix typically nIter by length(initPar) where the i-th row contains the estimate of the 
#'                     gradient of the synthetic likelihood at the i-th iteration (\code{matrix}).}
#'    \item{resultHess}{List of nIter elements where the i-th element is contains the estimate of the 
#'                     Hessian of the synthetic likelihood at the i-th iteration (\code{list}).}
#'    \item{resultCovar}{List of nIter elements where the i-th element is contains the estimate of the 
#'                     covariance matrix of the parameters at the i-th iteration (\code{list}).}
#'    \item{resultLoglik}{Vector of nIter elements where the i-th element is contains the estimate of the 
#'                         synthetic likelihood at the i-th iteration (\code{numeric}).}
#'  }
#' @name synMaxlik-class
#' @rdname synMaxlik-class
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @exportClass synMaxlik
.synMaxlik <-setClass("synMaxlik",
                      representation( initPar = "numeric",
                                      nIter   = "numeric",
                                      nsim    = "numeric",
                                      initCov   = "matrix",
                                      addRegr  = "logical",
                                      constr = "list",
                                      control = "list",
                                      continueCtrl = "list",
                                      multicore = "logical",
                                      ncores = "numeric",
                                    
                                      resultPar = "matrix",
                                      resultGrad = "matrix",
                                      resultHess  = "list",
                                      resultCovar = "list",
                                      resultLoglik = "numeric"
                    ),
                    prototype = prototype(initPar = numeric(),
                                          nIter   = numeric(),
                                          nsim    = numeric(),
                                          initCov   = matrix( , 0, 0),
                                          addRegr  = TRUE,
                                          constr = list(),
                                          control = list(),
                                          continueCtrl = list(),
                                          multicore = FALSE,
                                          ncores = detectCores() - 1,
                                          
                                          resultPar = matrix( , 0, 0),
                                          resultGrad = matrix( , 0, 0),
                                          resultHess = list(),
                                          resultCovar = list(),
                                          resultLoglik = numeric()
                                          ),
                    contains = "synlik",
                    validity = .check.synMaxlik)




### Validity check

.check.sml <- function(object)
{
  if(!is(object, "sml")) stop("object has to be of class \"sml\" ")
  
  errors <- character()
  
  if(length(object@initpar) == 0) errors <- c(errors, "length(initpar) should be > 0")
  
  if(length(errors) == 0) TRUE else errors
}


### Class Definition
#' \code{sml-class}
#'
#' \describe{
#'    \item{initPar}{Vector of initial parameters from which the optimization starts (\code{numeric}).}
#'    \item{nIter}{Number of iterations of the optimization routine (\code{numeric}).}
#'    \item{nsim}{Number of simulations from the simulator at each step of the optimization routine (\code{numeric}).}
#'    \item{initCov}{Covariance matrix used to simulate the parameter at each step of the optimization routine (\code{matrix}).}
#'    \item{constr}{Named list of contraints on the parameters, it has 3 elements: 
#'                  [["indexes"]] = (numeric integers) indexes of the elements to check;
#'                  [["upper"]]  = (numeric) upper bounds for the elements in "indexes";
#'                  [["lower"]]  = (numeric) lower bounds for the elements in "indexes" (\code{list}).}
#'    \item{multicore}{If TRUE the object@@simulator and object@@summaries functions will
#'                     be executed in parallel. That is the nsim simulations will be divided in multiple cores (\code{logical}).}
#'    \item{ncores}{Number of cores to use if multicore == TRUE (\code{numeric}).}
#'    \item{estim}{Matrix typically nIter by length(initPar) where the i-th row contains the estimate of the 
#'                 parameters at the i-th iteration (\code{matrix}).}
#'    \item{simPar}{Matrix typically nIter * nP by length(initPar) where the i-th row contains the estimate of the 
#'                 parameters at the i-th iteration (\code{matrix}).}
#'    \item{resultHess}{List of nIter elements where the i-th element is contains the estimate of the 
#'                     Hessian of the synthetic likelihood at the i-th iteration (\code{list}).}
#'    \item{resultCovar}{List of nIter elements where the i-th element is contains the estimate of the 
#'                     covariance matrix of the parameters at the i-th iteration (\code{list}).}
#'    \item{resultLoglik}{Vector of nIter elements where the i-th element is contains the estimate of the 
#'                         synthetic likelihood at the i-th iteration (\code{numeric}).}
#'  }
#' @name sml-class
#' @rdname sml-class
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @exportClass sml
.sml <-setClass("sml",
                      representation( initpar = "numeric",
                                      initcov   = "matrix",
                                      niter   = "numeric",
                                      nsim    = "numeric",
                                      np      = "numeric",
                                      priorfun = "functionOrNULL",
                                      alpha = "numeric",
                                      constr = "list",
                                      temper = "numericOrNULL",
                                      recycle = "logical",
                                      multicore = "logical",
                                      ncores = "numeric",
                                      
                                      estim = "matrix",
                                      simloglik = "numeric",
                                      simlogprior = "numeric",
                                      simpar  = "matrix"
                      ),
                      prototype = prototype(initpar = numeric(),
                                            initcov   = matrix( , 0, 0),
                                            niter   = numeric(),
                                            nsim    = numeric(),
                                            np      = numeric(),
                                            priorfun = NULL,
                                            alpha = 0.95,
                                            constr = list(),
                                            temper = NULL,
                                            recycle = FALSE,
                                            multicore = FALSE,
                                            ncores = detectCores() - 1,
                                            
                                            estim  = matrix( , 0, 0),
                                            simloglik = numeric(),
                                            simlogprior = numeric(),
                                            simpar = matrix( , 0, 0)
                                            
                      ),
                      contains = "synlik",
                      validity = .check.sml)
