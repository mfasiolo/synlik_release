#########
### Function used to set up the integrated thresholded rainfall as in Bhadhra et al 2011
#########
#' Set up the integrated thresholded rainfall as in Bhadhra et al, 2011.
#' 
#' @param rainfall Vector time series of monthly observed precipitations.
#' @param rainYear Time of precipitation [0, 20)
#' @param nSteps Number of multinomial-euler step in which each month is divided.
#' @param obsInterval Interval of time between 2 observations (typically 1/12 for monthly data).
#' @param lag Number of month used as lag (see Bhadhra et al 2011).
#' @param v Threshold used for the rainfall (see Bhadhra et al 2011).
#' 
#' @return vector of length( length(rainfall) * nSteps ) containing integrated rainfall.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>                         
#' @export setupRain

setupRain <- function(rainfall, rainYear, nSteps, obsInterval = 1/12, lag = 5, v = 200)
{
  par(mfrow = c(2, 2))
  plot(rainYear, rainfall, type = 'l', main = "Input rainfall")
  
  # Interpolating the rainfall
  interpolator <- approxfun(x = rainYear, y = rainfall, method = "linear")
  newx <- seq(rainYear[1], tail(rainYear, n = 1), length.out = length(rainfall) * nSteps)
  inRain <- interpolator(newx)
  plot(newx, inRain, type = 'l', main = "Interpolated rainfall" )
  
  # Setting lag month to be equal to the average rainfall
  nMonths <- 12
  totLen <- length(inRain)
  addRain <- rep(mean(inRain), lag * nSteps)
  inRain <- c(addRain, inRain)
  dt <- ( lag / nMonths ) / nSteps
  
  # Calculated integrated thresholded rainfall
  outRain <-  totLen
  jj = 1
  for(ii in (length(addRain)+1) : length(inRain))
  {
    start <- ii - lag*nSteps
    outRain[jj] <- max( sum(inRain[start:ii]*dt) - v, 0)
    jj <- jj + 1
  } 
  plot(newx, outRain, type = 'l', main = "Integrated thresholded rainfall")
  
  # Standardizing output
  outRain <- (outRain - mean(outRain)) / sd(outRain)
  
  plot(newx, outRain, type = 'l', main = "Output standardized integrated thresholded rainfall")
  
  return(outRain)
}