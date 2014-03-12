#########
### Simulator
#########

cholSimul  <- function(param, nsim, extraArgs, trans = FALSE, ...)
{
  if( is.null(extraArgs$pompModel)  ) stop("extraArgs$pompModel is missing")
  if( !is.loaded("pompExamples") ) library(pompExamples)
  
  pompModel <- extraArgs$pompModel
  
  names(param) <- names(pompModel@params)
  
  tmp <- simulate(pompModel, params = param, nsim = nsim)
  
  if(nsim > 1)
  {
    tmp <- t( sapply(tmp, function(input) input@data[1, ]) )
    
  } else {
    tmp <- tmp@data[1, ]
  }
  
  # Sometimes the simulator goes below zero!!
  tmp[tmp < 0] <- 0
  
  return( tmp )
}


cholSimulMulti  <- function(param, nsim, extraArgs, trans = FALSE, ...)
{
  if( is.matrix(param) )
  {
    tmp <- drop(cholSimul(param = param[1, ], nsim = 1, extraArgs, trans = FALSE, ...))
    
    out <- matrix(NA, nsim, length(tmp))
    out[1, ] <- tmp
    for(ii in 2:nsim) out[ii, ] <- drop(cholSimul(param = param[ii, ], nsim = 1, extraArgs, trans = FALSE, ...))
  }

  out  
}


########
## Statistics
########

cholStats  <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  if( !is.loaded("synlik") ) library("synlik")
  if( !is.loaded("matrixStats") ) library("matrixStats")
  
  obsData <- extraArgs$obsData
  
  stopifnot(is.vector(obsData), length(obsData) != 0)
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  qX <- x ^ 0.2
  tqX <- t(qX)
  
  nObs <- length(obsData)
  
  mu <- rowMeans(x)
  
  X0 <- cbind(# Capturing dynamics
    t(nlar(tqX, lag=c(1, 1, 2, 3, 3), power=c(1, 2, 1, 1, 2))), 
    t(sl.acf(tx, max.lag = 5))[ , -1, drop = FALSE],     #the first element is the variance (now removed)
    
    # Capturing distribution of the differences
    t(order.dist(tx, obsData, np=3, diff=1)),
    
    # Capturing skewness in the marginal distribution of the data
    mu, 
    (mu - rowMedians(x)) / mu,
    ( rowOrderStats(x, round(nObs*0.75)) - rowOrderStats(x, round(nObs*0.25)) ) / mu, # 450 - 150 = IQR(x) approximately, but this is faster
    rowOrderStats(x, round(nObs*0.9)),
    rowOrderStats(x, round(nObs*0.1)) / mu, 
    
    #Number of turning points
    apply( t( abs( apply( sign( apply(x, 1, diff) ), 2, diff ) )), 1, sum)/2 )
  
  # Capturing mean, time trend and dominant frequencies
  tim <- (1:nObs) / 12
  modelMat <- cbind(1, tim, sin(0.12*2*pi*tim), cos(0.12*2*pi*tim),  
                            sin(1*2*pi*tim)   , cos(1*2*pi*tim),
                            sin(2*2*pi*tim)   , cos(2*2*pi*tim),
                            cos(3*2*pi*tim)   , sin(4*2*pi*tim) )
  tmp <- unname( lm.fit(x = modelMat, y = tqX)$coef )
  X0 <- cbind(unname(X0), t(tmp)[ , -1, drop = FALSE])
  
  return(X0)
}  ## end of ss.model.matrix function
