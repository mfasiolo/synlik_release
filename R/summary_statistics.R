############################################################################################################################################
##################################               S U M M A R Y    S T A T I S T I C S                     ##################################
############################################################################################################################################
############ All individual and grouped summary statistics are stored here



###########################################################################################################################################
##########################
######## INDIVIDUAL SUMMARY STATISTICS
##########################

sl.acf <- function(x, max.lag=10) {
  ## `x' is a matrix containing replicate simulations in its columns.
  ## sl.acf turns these into acf's
  
  NAcode <- -1e70
  x[is.na(x)] <- NAcode
  
  acf <- matrix(0,max.lag+1,ncol(x))
  oo <- .C("slacf",acf=as.double(acf),x=as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),
           as.integer(max.lag),as.double(NAcode),correlation=as.integer(0),PACKAGE="synlik")
  
  acf <- matrix(oo$acf,max.lag+1,ncol(x))
  acf[acf == NAcode] <- NA
  acf
  
} ## end of sl.acf

nlar <- function(x,lag,power) {
  ## relatively efficient polynomial autoregression for multiple reps.
  ## each column of `x' is a replicate. 
  ## `lag[i]' is the lag for term i on rhs of autoregression
  ## `power[i]' is the power for term i on rhs of autoregression 
  beta <- matrix(0,length(lag),ncol(x))
  
  NAcode <- -1e70
  x[is.na(x)] <- NAcode  
  
  oo <- .C("slnlar",beta = as.double(beta), x = as.double(x),
           n=as.integer(nrow(x)),n.reps=as.integer(ncol(x)),n.terms=as.integer(length(lag)),
           as.integer(lag),as.integer(power),as.double(NAcode),PACKAGE="synlik")
  
  beta <- matrix(oo$beta,length(lag),ncol(x))
  
  beta
} ## end of nlar


order.dist <- function(x,z,np=3,diff=1) {
  ## Routine to obtain coefficients summarizing distribution of (differenced) columns
  ## of x, by regression of sorted differenced columns of x on sorted differenced z's. 
  ## regression is with order `np' polynomial (no intercept as all centred). `diff'
  ## is order of differencing to apply.
  
  beta <- matrix(0,np,ncol(x))
  oo <- .C("order_reg",beta=as.double(beta), as.double(x),as.double(z),n=as.integer(nrow(x)),
           as.integer(ncol(x)),as.integer(np),as.integer(diff),PACKAGE="synlik")
  
  beta <- matrix(oo$beta,np,ncol(x))
  
} ## end of order.dist


# FP statistics
fpstat <- cmpfun(function(Y) {
  
  n <- dim(Y)[1]
  n.paths <- dim(Y)[2]
  res = matrix(NA,n.paths,3) 
  
  R <- log(Y[2:n, , drop = FALSE]) - log(Y[1:(n-1), , drop = FALSE])
  R[which(!is.finite(R))] <- NA
  X <- Y[1:(n-1), , drop = FALSE]
  X[which(is.na(R))] <- NA
  
  # Estimating r and phi
  xcm <- colMeans(X, na.rm=TRUE)
  X <- sweep(X, 2, xcm)
  beta <- (colSums(R*X,na.rm=TRUE)/colSums(X^2,na.rm=TRUE))
  r.hat <- colMeans(R,na.rm=TRUE) - xcm*beta
  phi.hat <- -1/beta
  beta
  
  # This is the estimator for sigma as in the discussion, but it doesn't work that well.
  # The sigmas estimated to be < 0 are resampled from those that are > 0.
  e <- sweep(R, 2, r.hat) + sweep(Y[1:(n-1), , drop = FALSE], 2, phi.hat, "/")
  v <- (1/Y[2:n, , drop = FALSE]) + sweep(-1/Y[1:(n-1), , drop = FALSE], 2, 1/phi.hat, "+")^2 * Y[1:(n-1), , drop = FALSE]
  e2 <- e^2 - v
  sig.hat <- colMeans(e2, na.rm=TRUE)
  badSigmas <- which(sig.hat < 0)
  nBad <- length(badSigmas)
  if(nBad > 0){ 
    if(nBad == n.paths) stop("Estimated sigma is always below negative")
    sig.hat[badSigmas] <- sample(sig.hat[-badSigmas], nBad, replace = TRUE)
  }
  
  sig.hat <- sig.hat^.5
  res <- unname(cbind(r.hat, log(pmax(sig.hat, 10^-3)), log(pmax(phi.hat, 10^-3))))
  
  res
})








###########################################################################################################################################
##########################
######## Grouped SUMMARY STATISTICS
##########################

########
### WOOD2010 original statistics
########

WOOD2010 <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- as.vector(extraArgs$obsData)
  stopifnot(length(obsData) != 0)
  
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(order.dist(tx, obsData, np=3, diff=1))        ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(tx^.3, lag=c(1,1), power=c(1,2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x), rowSums(x==0))         ## mean values of Y, # of 0's
  X0 <- cbind(X0, t(sl.acf(tx, max.lag=5)))                 ## autocovariances up to lag 5 (the first element is the variance)
  
  X0
}

########
### WOOD2010 + FP STATISTICS
########

woodPlusFP <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- as.vector(extraArgs$obsData)
  stopifnot(length(obsData) != 0)
  
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(order.dist(tx, obsData, np=3, diff=1))        ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(tx^.3, lag=c(1,1), power=c(1,2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x), rowSums(x==0))         ## mean values of Y, # of 0's
  X0 <- cbind(X0, t(sl.acf(tx, max.lag=5)))                 ## autocovariances up to lag 5 (the first element is the variance)
  X0 <- cbind(X0, fpstat(tx))         
  
  X0
}


#########################
####### Blowfly Statistics
#########################

blowStats <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- as.vector(extraArgs$obsData)
  stopifnot(length(obsData) != 0)
  
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(order.dist(t(x), obsData,np=3,diff=1)) ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(t(x),lag=c(6, 6, 6, 1, 1),power=c(1,2,3,1,2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x), rowMeans(x) - apply(x,1,median)) 
  X0 <- cbind(X0, t(sl.acf(t(x), max.lag=11)))  ## autocovariances up to lag 5 (the first element is the variance)
  
  X0 <- cbind(X0, apply( t( abs( apply( sign( apply(x,1,diff) ), 2, diff ) )), 1, sum)/2) #Number of turning points 
  
  X0
}



########
### Statistics for stable distribution
########

stableStats <- function(x, extraArgs, ...){
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  ## build a model matrix for linear regression
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  
  X0 <- t( apply(x, 1, quantile, probs = seq(0.1, 0.9, length.out = 15)) )  # Sample quantiles
  
  unname(X0)
}


####################
######## VOLES: statistics for the vole model (Turchin)
####################

volesStats  <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- extraArgs$obsData
  
  stopifnot(is.vector(obsData), length(obsData) != 0)
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(order.dist(tx, obsData, np=3, diff=1)) ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(t(x),lag=c(6, 6, 6, 1, 1), power=c(1, 2, 3, 1, 2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x)) # mean values of Y, # of 0's
  X0 <- cbind(X0, t(sl.acf(tx, max.lag = 5)))  ## autocovariances up to lag 5 (the first element is the variance)
  X0 <- cbind(X0, apply( t( abs( apply( sign( apply(x,1,diff) ), 2, diff ) )), 1, sum)/2) #Number of turning points 
  #X0 <- t( log(apply(x, 1, function(input) spec.pgram(input, plot = FALSE)$spec)) )[ , 1:20, drop = FALSE] # Spectrogram
  
  return(X0)
}  ## end of ss.model.matrix function


volesStatsLarge  <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- extraArgs$obsData
  
  stopifnot(is.vector(obsData), length(obsData) != 0)
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(order.dist(tx, obsData, np=3, diff=1)) ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(t(x),lag=c(15, 15, 10, 10, 6, 6, 6, 3, 3, 1, 1), power=c(1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x)) # mean values of Y, # of 0's
  X0 <- cbind(X0, t(sl.acf(tx, max.lag = 10)))  ## autocovariances up to lag 5 (the first element is the variance)
  X0 <- cbind(X0, apply( t( abs( apply( sign( apply(x, 1, diff) ), 2, diff ) )), 1, sum)/2) #Number of turning points 
  #X0 <- cbind(X0, fpstat(tx)[ , -2, drop = FALSE]) 
  X0 <- cbind(X0, unname(t(apply(x, 1, quantile, probs = seq(0.01, 0.99, by = 0.05)))))
  X0 <- cbind(X0, apply(x, 1, function(input) sd(diff(input, 1))), apply(x, 1, function(input) sd(diff(input, 3))))
  #X0 <- t( log(apply(x, 1, function(input) spec.pgram(input, plot = FALSE)$spec)) )[ , 1:20, drop = FALSE] # Spectrogram
  
  # From freanhead
  X0 <- cbind(X0, X0^2)
  #X0 <- cbind(X0, log(rowSums(x^2)), log(rowSums(x^3)), log(rowSums(x^4)), log(rowSums(x^5)), log(rowSums(x^6)))
  
  return(X0)
}  ## end of ss.model.matrix function


###############
#### Lokta statistics
##############


# It's just the squared distance between simulated and observed paths
lvNormStats <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- extraArgs$obsData
  
  stopifnot(is.vector(obsData), length(obsData) != 0)
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- matrix(rowMeans( (sweep(x, 2, obsData))^2 ), nrow(x), 1)
  
  X0
}

loktaStats <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- extraArgs$obsData
  
  stopifnot(is.vector(obsData), length(obsData) != 0)
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(order.dist(tx, obsData, np=3, diff=1)) ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(order.dist(tx, obsData, np = 1, diff = 20))) ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(t(x), lag=c(1, 1, 6, 6, 6, 10, 10, 15, 15, 25, 25), power=c(1, 2, 1, 2, 3, 1, 2, 1, 2, 1, 2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x), apply(x, 1, sd)) # mean values of Y, # of 0's
  X0 <- cbind(X0, unname(t(apply(x, 1, quantile, probs = seq(0.1, 0.9, by = 0.1)))))
  X0 <- cbind(X0, t(sl.acf(tx, max.lag = 5)))  ## autocovariances up to lag 5 (the first element is the variance)
  #X0 <- t( log(apply(x, 1, function(input) spec.pgram(input, plot = FALSE)$spec)) )[ , 1:20, drop = FALSE] # Spectrogram
  
  return(X0)
}  ## end of ss.model.matrix function




####################
######## Malaria: statistics for the malaria model
####################

malariaStats <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- extraArgs$obsData
  theRain <- extraArgs$theRain
  
  stopifnot(is.vector(obsData), length(obsData) != 0, is.vector(theRain), length(theRain) != 0)
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  ####
  # Targeting the CONDITIONAL distribution of the data 
  ####
  # One month ahead
  X0 <- t(order.dist(tx, obsData, np=3, diff=1)) # cubic regs coeff of ordered diffs
  
  # 3 months ahead
  mAhead <- 3
  X0 <- cbind(X0, unname( t( apply(x, 
                                   1,
                                   function(input)
                                   {
                                     quantile( diff(input, mAhead), c(0.1, 0.3, 0.5, 0.7, 0.9) )  
                                   }))))
  #####
  # Targeting the MARGINAL distribution of the data
  #####
  X0 <- cbind(X0, unname(t(apply(x, 1, quantile, probs = c(0.1, 0.3, 0.5, 0.7, 0.9)))))
  
  ####
  # Targeting the time dependence
  ####
  # Non linear auto-regressive
  X0 <- cbind(X0, t(nlar(tx,lag=c(20, 15, 10, 6, 6, 6, 3, 3, 1, 1), power=c(1, 1, 1, 1, 2, 3, 1, 2, 1, 2))))
  # Autocovariance (the first element is the variance)
  X0 <- cbind(X0, t(sl.acf(tx, max.lag = 20))[ , c(1, 3, 5, 7, 10, 15, 20), drop = FALSE]) 
  
  ####
  # Targeting the dependence on rainfall
  ####
  # Regressing number of cases_t on rain_t, rain_(t-1) and rain_(t-2)    
  myLag <- 3
  lenRain <- length(theRain)
  modelMat <- cbind(1, theRain[-(1:myLag)], 
                    theRain[myLag:(lenRain - 1)], 
                    theRain[(myLag-1):(lenRain - 2)], 
                    theRain[(myLag-2):(lenRain - 3)])
  tmp <- unname( lm.fit(x = modelMat, y = tx[-(1:myLag), ])$coef )
  
  if( !is.matrix(tmp) ) tmp <- matrix(tmp, length(tmp), 1)
  X0 <- cbind(X0, t(tmp[2:(myLag+2), ]))
  
#   # Look at the regression if you want
#   llag <- 3 # From 0 to myLag
#   plot(theRain[1:(240 - llag)], drop(malaria_sl@data)[(1+llag):240])
#   lines(theRain[1:(240 - llag)], tmp[1] + tmp[llag+2]*theRain[1:(240 - llag)])
  
  ####
  # Other stats
  ####
  #Number of turning points 
  X0 <- cbind(X0, apply( t( abs( apply( sign( apply(x, 1, diff) ), 2, diff ) )), 1, sum)/2)  

  return(X0)
}  ## end of ss.model.matrix function



##########################
###### Exponential Distribution Stats
##########################

expStats <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  
  X0 <- matrix( 1 / rowMeans(x), nrow(x), 1)
  
  return(X0) 
}