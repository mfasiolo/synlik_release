##########
#### Method to plot the object
##########
#
# ARGS:
# x = object of class "synMaxlik"
# trans = (list) that control the trasformations of the parameters to be plotted
#         trans$parIndex = (numeric integers) indexes of the parameters to be transformed (ex: c(1, 3))
#         trans$transform = (character vector) _names_ of the functions to be applied to each paramters (ex: c("exp", "log")) 

plot.synMaxlik <- function(x, trans = NULL, Lag = 10, tol = 1e-8, verbose = TRUE, ...)
{
  if(!is(x, "synMaxlik")) stop("object has to be of class \"synMaxlik\" ")
  
  if(!is.null(trans)) stopifnot( identical(names(trans), c("parIndex", "transform")), all(is.function(get(trans$transform))),
                                 length(trans$parIndex) == length(trans$transform) )
  
  parEstim <- x@resultPar
  nIter <- nrow(parEstim)
  
  if(nIter == 0) {
    warning("object@resultPar is empty (i.e. there are no results to plot)")
    invisible( callNextMethod(x) )
  }
  
  ####
  # Plot of the convergence of the parameter's estimates
  ####
  
  nPar <- ncol(parEstim)
  panelDim <- min( ceiling(sqrt(nPar)), 3 )
  par(mfrow = c(panelDim, panelDim))
  
  readline(prompt = "Press <Enter> to see the next plot...")
  cat("Convergence of parameters' estimates")
  
  # The hessian is an average of the last "Lag" hessians
  tmp <- coef.synMaxlik(object = x, Lag = Lag, tol = tol)
  confInt <- tmp$par
  covariance <- tmp$cov
  
  jj <- 1
  for(iPar in 1:nPar){
    
    parVals <- parEstim[ , iPar]
    finalPar <- tail(parVals, 1)
    
    valid <- !is.na( confInt[iPar, 1] )
    
    # Getting the parameters transformation, the default it the identity
    if(iPar %in% trans$parIndex)
    {
      transFun <- get(trans$transform[jj])
      jj <- jj + 1
    } else {
      transFun <- function(input) input
    }
    
    # Calculate confidence intervals if the stderror is positive
    if(valid)
    {
      stdErr <- diff( confInt[iPar, 1:2] ) / 1.96
      yVals <- seq(finalPar - 3*stdErr, finalPar + 3*stdErr, length.out = 100)
      llk <- dnorm(yVals,  mean = finalPar, sd = stdErr, log = TRUE)
      llk <- llk - max(llk)
      ylim <- c(min(parVals, finalPar - 3 * stdErr), max(parVals, finalPar + 3 * stdErr))
    }
    
    plot(1:nIter, 
         y = transFun(parVals), 
         type = 'l', 
         ylim = if(valid) transFun(ylim) else NULL, 
         main = names(x@param)[iPar],
         ylab = names(x@param)[iPar], 
         xlab = "Iteration")
    
    # Add confidence intervals if the stderror is positive
    if(valid)
    { 
      lines(nIter + llk, transFun(yVals))
      abline(h = transFun( confInt[iPar, 1] ), lty = 2)
      abline(h = transFun( confInt[iPar, 3] ), lty = 2)
    }
    
    if( !(iPar %% (panelDim^2)) && (iPar != nPar) ) readline(prompt = "Press <Enter> to see the next plot...") 
    
  }
  
  ####
  # Plot of correlation matrix of the parameters estimates
  ####
  
  readline(prompt = "Press <Enter> to see the next plot...")
  cat("Estimated asymptotic correlation matrix of the estimated parameters")
    
  par(mfrow = c(1, 1))
  saveOpt <- par()$mar
  finalCorr <- extractCorr(covariance)
  .plotMatrix(finalCorr, 
             title = "Estimated correlation matrix", 
             xLabels = colnames(finalCorr), 
             yLabels = colnames(finalCorr), 
             scaleLab = "Correlation", correl = TRUE)
  suppressWarnings( par(mar = saveOpt) )

  ####
  # Plot of evolution of the diagonal entries of the Hessian
  ####
  
  if(verbose)
  {
    readline(prompt = "Press <Enter> to see the next plot...")
    cat("Convergence of the diagonal entries of the Hessian matrix")
    
    par(mfrow = c(panelDim, panelDim))
    saveOpt <- par()$mar
    par(mar = c(5.1, 4.5, 4.1, 2.1))
    hessEstim <- x@resultHess
    for(iPar in 1:nPar){
      tmp <- sapply(hessEstim, function(input) input[iPar, iPar])
      plot(1:nrow(parEstim), tmp, main = names(x@param)[iPar],
           ylab = substitute("d" ^ 2 * "log-lik / d " * theta^2, list(theta = names(x@param)[iPar])), xlab = "Iteration")
      lines( .rollMean(vett = tmp, Lag = min(10, nIter)), col = 2 )
      
      if( !(iPar %% (panelDim^2)) && (iPar != nPar) ) readline(prompt = "Press <Enter> to see the next plot...") 
      
    }
    suppressWarnings( par(mar = saveOpt) )
  }
  
  ####
  # Plot of the convergence of the log-likelihood (objective function)
  ####
  
  readline(prompt = "Press <Enter> to see the next plot...")
  cat("Plotting the evolution of the log-likelihood (objective function)")
  
  par(mfrow = c(1, 1))
  plot(1:nrow(parEstim), x@resultLoglik, xlab = "Iteration", ylab = "Log-likelihood", main = "Log-likelihood evolution", type = 'l')
  
  
  ####
  # Plot of the convergence of the covariance matrix of the local model
  ####
  
  if(verbose)
  {
    readline(prompt = "Press <Enter> to see the next plot...")
    cat("Plotting the convergence of the log-variances used for the local model")
    
    par(mfrow = c(panelDim, panelDim))
    covarEstim <- x@resultCovar
    for(iPar in 1:nPar){
      tmp <- sapply(covarEstim, function(input) input[iPar, iPar])
      plot(1:nrow(parEstim), log(tmp), type = 'l', main = names(x@param)[iPar],
           ylab = paste("log-variance of", names(x@param)[iPar]), xlab = "Iteration")
      abline(h = log(x@control$limCov$upper[iPar]), col = "green")
      abline(h = log(x@control$limCov$lower[iPar]), col = "red")
      if( !(iPar %% (panelDim^2)) && (iPar != nPar) ) readline(prompt = "Press <Enter> to see the next plot...") 
    }
  }
  
  ####
  # Plot of the convergence of the correlation matrix" of the local model
  ####
  
  if(verbose)
  {
    readline(prompt = "Press <Enter> to see the next plot...")
    cat("Plotting final correlation matrix of the local model")
    
    par(mfrow = c(1, 1))
    saveOpt <- par()$mar
    finalCovar <- tail(covarEstim, 1)[[1]]
    finalCorr <- diag(sqrt(diag(finalCovar))^-1, nPar)%*%finalCovar%*%diag(sqrt(diag(finalCovar))^-1, nPar)
    .plotMatrix(finalCorr, title = "Final correlation matrix of the local model", xLabels = names(x@param), yLabels = names(x@param), 
               scaleLab = "Correlation", correl = TRUE)
    suppressWarnings( par(mar = saveOpt) )
    
  }
  
  readline(prompt = "Press <Enter> to see the next plot...")
  invisible( callNextMethod(x) )
  
}

setMethod("plot",
          signature = signature(x = "synMaxlik", y = "missing"),
          definition = plot.synMaxlik
)


