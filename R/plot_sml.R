#########
#### Method to plot the object
##########
#
# ARGS:
# x = object of class "synMcmc"
# trans = (list) that control the trasformations of the parameters to be plotted
#         trans$parIndex = (numeric integers) indexes of the parameters to be transformed (ex: c(1, 3))
#         trans$transform = (character vector) _names_ of the functions to be applied to each paramters (ex: c("exp", "log")) 

plot.sml <- function(x, trans = NULL, addPlot = NULL, ...)
{
  if(!is(x, "sml")) stop("object has to be of class \"synMcmc\" ")
  estim <- x@estim
  addPoints <- list("x" = rep(1:x@niter, each = x@np), "y" = x@simPar)
  
  print("Plotting the convergence trajectories")
  invisible( .plotIter(estim, trans, type = "line", addPoints = addPoints, addPlot = addPlot, ...) )
}


setMethod("plot",
          signature = signature(x = "sml", y = "missing"),
          definition = plot.sml
)
