#####
#' Finding the mode of the empirical saddelpoint density
#' @description Given a sample from a d-dimensinal distribution the routine
#'              find the mode of the corresponding empirical saddlepoint density.
#' @param X n by d matrix containing the data.
#' @param init d-dimensional vector containing the starting point for the optimization. By default
#'             it is equal to colMeans(X).
#' @param decay Rate at which the SPA falls back on a normal density. Should be a positive number,
#'              by default set to 0.5.
#' @param method Optimization method used by stats::optim(), see ?optim for details. By default it is "BFGS".
#' @param sadTol Tolerance used to assess the convergence of the rootfinding routine used to fit
#'               the saddlepoint density. Default value is 1e-6.
#' @param ... Extra arguments to be passed to the optimization routine stats::optim. 
#' @return A list where \code{mode} is the location of mode of the empirical saddlepoint density,
#'         while the other entries are the same as for stats::optim.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @export
#'

findMode <- function(X, init = NULL, decay = 0.5, method = "BFGS", 
                     sadTol = 1e-6, marginal = FALSE, draw = FALSE, ...)
{
  switch(class(X),
         "matrix"  = theData <- X,
         "numeric" = theData <- matrix(X, length(X), 1),
         stop("X should be either of class \"vector\" or a matrix")
  )
  
  nDims <- ncol(theData)
  nObs <- nrow(theData)
  
  if( nDims > nObs )
  {
    warning("The input data should be seems to have more dimensions (columns) than data points (rows): I'm transposing it!")
    theData <- t(theData)
  }
  
  if(is.null(init)){ 
    init <- colMeans(theData) 
  } else { 
    stopifnot( is.vector(init), length(init) == ncol(theData) ) 
  }
  
  objFun  <- function(x) -dsaddle(y = as.numeric(x), X = theData, tol = sadTol, decay = decay, log = TRUE)$llk
  objGrad <- function(x) -dsaddle(y = as.numeric(x), X = theData, tol = sadTol, decay = decay, deriv = TRUE)$grad
  
  optOut <- optim(par = init, fn = objFun, gr = objGrad, method = method, ...)
  
  names(optOut)[c(1, 2, 3)] <- c("mode", "log-density", "number of evaluations")
  
  if(draw == TRUE)
  {
    nPoints <- 100
    dens  <- evalPoints <- numeric(nPoints)
    for(ii in 1:nDims)
    {
      evalPoints <- seq(min(theData[ , ii]), max(theData[ , ii]), length.out = nPoints)
      densFun <- function(x) dsaddle(y = as.numeric(x), X = theData, tol = sadTol, decay = decay)$llk
      
      for(kk in 1:nPoints)
      {
        x <- optOut$mode
        x[ii] <- evalPoints[kk]
        dens[kk] <- densFun(x = x)
      }
      readline(prompt = "Press <Enter> to see the next density...")
      
      par(mfrow = c(2, 1))
      hist(theData[ , ii], 
           main = if(class(input) == "synMcmc") paste("Marginal posterior sample for parameter", names(input@params)[ii]) else paste("Marginal data for variable", ii),
           xlab = if(class(input) == "synMcmc") names(input@params)[ii] else "Variable value")
      
      plot(evalPoints, dens, type = 'l',
           main = if(class(input) == "synMcmc") paste("Transect of log-saddlepoint density wtr parameter", names(input@params)[ii]) else paste("Transect of log-saddlepoint density wtr variable", ii), 
           xlab = if(class(input) == "synMcmc") names(input@params)[ii] else "Variable value", 
           ylab = "Log-saddlepoint density")
      abline(v = optOut$mode[ii], lty = 2)
    }
  }
  
  return(optOut)
}