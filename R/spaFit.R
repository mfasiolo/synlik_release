###########################
###### Univariate Spa fit
###########################

spaFit <- function(x, theData, decay = 0.5, draw = TRUE, log = FALSE, ...)
{
  stopifnot( is.vector(x) )
  if( !is.vector(theData) ) theData <- as.vector(theData)
  
  llk <- x * 0
  for(ii in 1:length(x)) llk[ii] <- dsaddle(y = x[ii], X = theData, decay = decay, log = log, ...)$llk
  
  plot(x, llk, type = 'l', main = c("Black = saddle, red = normal"))
  lines(x, dnorm(x, mean(theData), sd(theData)), col = 2)
  
  invisible( llk )
}