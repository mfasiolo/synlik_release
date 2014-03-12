########
#' Parameters covariance matrix for "sml" objects
#'
#' @param object    a fitted model object of class "sml".
#' @param parm      a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. 
#'                  If missing, all parameters are considered.
#' @param level     the confidence level required.
#' @param lag     (integer) final estimate of the parameter is the mean of the 
#'                last "lag" iterations. Hence lag <= object@niter.
#' @description Computes confidence intervals for one or more parameters in a fitted model of class "sml".
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @export
#' 

confint.sml <- function(object, parm, level = 0.95, lag = 10, ...)
{
  # Code copies from stats::confint.default()
  cf <- coef.sml(object, lag = lag, ...)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(vcov.sml(object, ...)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}

setMethod("confint",
          signature = signature(object = "sml"),
          definition = confint.sml
)









