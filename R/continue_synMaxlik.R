#########
### Method to continue synMaxlik estimation of a synMaxlik object
#########

# For initPar, initCov amd control unless they have been specified by the user, we use the final 
# values contained in object (ex: initPar is the final value in "object")
continue.synMaxlik <- function(object, 
                               nIter = object@nIter, 
                               nsim  = object@nsim,
                               initCov = NULL, 
                               addRegr = object@addRegr, 
                               constr = object@constr, 
                               control = list(), 
                               multicore = object@multicore,
                               ncores = object@ncores,
                               cluster = NULL,
                               ...)
{
  if(!is(object, "synMaxlik")) stop("To use mcmc you need an object of class \"synMaxlik\"")
  
  ### Setting up control
  # If the user provided some inputs in "control" include them in tmpControl and object@control
  tmpControl <- .ctrlSetup(innerCtrl = object@continueCtrl, outerCtrl = control)
  object@control <- .ctrlSetup(innerCtrl = object@control, outerCtrl = control)  
  
  tmpObject <- synMaxlik(object, 
                         initPar = as.vector( tail(object@resultPar, 1) ), 
                         nIter = nIter, 
                         nsim = nsim,
                         initCov = if( is.null(initCov) ) tail(object@resultCovar, 1)[[1]] else initCov,  
                         addRegr = addRegr,  
                         constr = constr, 
                         control = tmpControl, 
                         multicore = multicore,
                         ncores = ncores,
                         cluster = cluster,
                         ...)
  
  
  return(.synMaxlik(tmpObject,
                    
                    initPar = object@initPar,  # Resetting the values of these three param, so we don't lose information
                    control = object@control, 
                    initCov = object@initCov,
                    
                    nIter = object@nIter + nIter,
                    
                    resultPar = rbind(object@resultPar, tmpObject@resultPar),
                    resultGrad = rbind(object@resultGrad, tmpObject@resultGrad),
                    resultHess  = append(object@resultHess, tmpObject@resultHess),
                    resultCovar = append(object@resultCovar, tmpObject@resultCovar),
                    resultLoglik = append(object@resultLoglik, tmpObject@resultLoglik)
  ))
}


setMethod("continue", 
          signature = signature(object = "synMaxlik"), 
          definition = continue.synMaxlik)