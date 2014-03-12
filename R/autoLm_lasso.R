# Modified version of lars::lars() that automatically finds the penalty
# that minimized GCV by root finding. 
# ARGS:
# The arguments are the same as for lars() and cv.lars(),  but I removed some options.
# For details see ?lars::lars() and ?lars::cv.lars
.autoLm.lasso <- function(x, y, verbose = TRUE, maxSteps = 8 * min(ncol(x), nrow(x-1)))
{
  object <- cv.lars(x, y, plot.it = verbose, max.steps = maxSteps)
  
  cv <- object$cv
  errors <- object$cv.error
  index <- object$index
  
  best <- which.min(cv)
  upper <- cv[best] + errors[best]
  
  # Select the most parsimonious model within one sd error from the best cross validation
  best <- which( cv[1:best] < upper )[1]
  abline(v = index[best], lty = 2, col = 2)
  
  # Fit again to get to get a "lars" object
  object <- lars(x, y, type = "lasso", max.steps = maxSteps)
  
  # Get the coefficients for best fraction, calculate intercept and get fitted values
  tmpCoef <- coef(object, s = index[best], mode = "fraction")
  inter <- drop( object$mu - t(tmpCoef) %*% object$meanx )
  fit <- inter + drop( x%*%tmpCoef )
  
  list("coef" = c(inter, tmpCoef), "fit" = fit, "fraction" = index[best] )
}



