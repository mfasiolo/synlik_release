########
#' Parameters covariance matrix for "sml" objects
#'
#' @param object  ("sml") object.
#' @param nreps   (integer) number of simulations used to tilt the negative Hessian (-H)
#'                toward positive definiteness. Used only if -H is not PD.
#' @param boot    (logical) relevant only if -H is not PD. 
#'                If TRUE hessians will be simulated by resampling parameters and likelihoods. 
#'                If FALSE hessians will be simulated the asymptotic distribution of the regression
#'                coefficients.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @method vcov sml
#' @export
#' 

vcov.sml <- function(object, nreps = 1000, boot = TRUE, ...)
{ 
  # Including only parameters that were estimated
  varPar <- diag(object@initCov) > 0
  
  covar <- .vcov.ml(llk = object@simLogLik, 
                    parMat = object@simPar[ , varPar],
                    nreps = nreps, 
                    boot = boot, ...)
  
  rownames(covar) <- colnames(covar) <- names(object@param)[ varPar ]
  
  return(covar) 
}

setMethod("vcov",
          signature = signature(object = "sml"),
          definition = vcov.sml
)

