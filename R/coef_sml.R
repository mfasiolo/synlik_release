#########
#' Parameters estimates for "sml" objects
#'
#' @param object  ("sml") object.
#' @param lag     (integer) final estimate of the parameter is the mean of the 
#'                last "lag" iterations. Hence lag < object@@niter.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @method coef sml
#' @export
#' 

coef.sml <- function(object, lag = 10, ...)
{
  est <- colMeans( tail(object@estim, min(object@niter, lag)) )
  names(est) <- names(object@param)
  
  return(est)
}

setMethod("coef",
          signature = signature(object = "sml"),
          definition = coef.sml
)