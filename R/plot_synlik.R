##########
#### Method to plot the object
##########

plot.synlik <- function(x)
{
  if(!is(x, "synlik")) stop("object has to be of class \"synlik\" ")
  
  if(is.numeric(x@data))
  {
    if(length(x@data) > 0){
      cat("The data")
      plot(1:length(x@data), x@data, type = 'l', ylab = "Data", xlab = "Observation", main = "The data")
    }
  } else {
    warnings("plot.synlik can plot only data of class \"numeric\"")
  }
  
}

setMethod("plot",
          signature = signature(x = "synlik", y = "missing"),
          definition = plot.synlik
)