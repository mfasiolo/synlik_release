#########
#### Method to plot the object
##########
#
# ARGS:
# x = object of class "synMcmc"
# trans = (list) of transforms for each parameters (ex: list("par1" = "exp", "par2" = "log")) 

plot.synMcmc <- function(x, trans = NULL, addplot1 = NULL, addplot2 = NULL, ...)
{
  if(!is(x, "synMcmc")) stop("object has to be of class \"synMcmc\" ")
  chains <- x@mcmcChain
  colnames(chains) <- names(x@param)
  
  if(nrow(chains) > 0 )
  {
    print("Plotting the MCMC chains")
    .plotIter(chains, trans, type = "line", addplot = addplot1, ...)
    
    readline(prompt = "Press <Enter> to see the next plot...")
    
    print("The posterior densities")
    .plotIter(chains, trans, type = "hist", addplot = addplot2, ...)
    
    readline(prompt = "Press <Enter> to see the next plot...")
    
    print("Plotting the log-likelihood chain")
    par(mfrow = c(1, 1))
    plot(1:nrow(chains), x@LogLikChain, xlab = "Iteration", ylab = "Log-likelihood", main = "Log-likelihood chain", type = 'l', ...)
  }
  
  readline(prompt = "Press <Enter> to see the next plot...")
  
  print("Plotting correlation structure of the posterior sample")
  plotMatrix(cor(chains), title = "Correlations", xLabels = names(x@param), yLabels = names(x@param), 
             scaleLab = "Correlation", correl = TRUE, ...)
  
  readline(prompt = "Press <Enter> to see the next plot...")
  
  invisible( callNextMethod(x) )
}


setMethod("plot",
          signature = signature(x = "synMcmc", y = "missing"),
          definition = plot.synMcmc
)
