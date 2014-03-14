############################################################################################################################################
##################################               S U M M A R Y    S T A T I S T I C S                     ##################################
############################################################################################################################################
############ All individual and grouped summary statistics are stored here


##########################
######## INDIVIDUAL SUMMARY STATISTICS
##########################

slAcf <- function(x, max.lag=10) {
  ## `x' is a matrix containing replicate simulations in its columns.
  ## slAcf turns these into acf's
  
  NAcode <- -1e70
  x[is.na(x)] <- NAcode
  
  acf <- matrix(0,max.lag+1,ncol(x))
  oo <- .C("slacf",acf=as.double(acf),x=as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),
           as.integer(max.lag),as.double(NAcode),correlation=as.integer(0),PACKAGE="synlik")
  
  acf <- matrix(oo$acf,max.lag+1,ncol(x))
  acf[acf == NAcode] <- NA
  acf
  
} ## end of slAcf

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


orderDist <- function(x,z,np=3,diff=1) {
  ## Routine to obtain coefficients summarizing distribution of (differenced) columns
  ## of x, by regression of sorted differenced columns of x on sorted differenced z's. 
  ## regression is with order `np' polynomial (no intercept as all centred). `diff'
  ## is order of differencing to apply.
  
  beta <- matrix(0,np,ncol(x))
  oo <- .C("order_reg",beta=as.double(beta), as.double(x),as.double(z),n=as.integer(nrow(x)),
           as.integer(ncol(x)),as.integer(np),as.integer(diff),PACKAGE="synlik")
  
  beta <- matrix(oo$beta,np,ncol(x))
  
} ## end of orderDist




##########################
######## Grouped SUMMARY STATISTICS
##########################

########
### WOOD2010 original statistics
########

rickerStats <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- as.vector(extraArgs$obsData)
  stopifnot(length(obsData) != 0)
  
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(orderDist(tx, obsData, np=3, diff=1))        ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(tx^.3, lag=c(1,1), power=c(1,2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x), rowSums(x==0))         ## mean values of Y, # of 0's
  X0 <- cbind(X0, t(slAcf(tx, max.lag=5)))                 ## autocovariances up to lag 5 (the first element is the variance)
  
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
  
  X0 <- t(orderDist(t(x), obsData,np=3,diff=1)) ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(t(x),lag=c(6, 6, 6, 1, 1),power=c(1,2,3,1,2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x), rowMeans(x) - apply(x,1,median)) 
  X0 <- cbind(X0, t(slAcf(t(x), max.lag=11)))  ## autocovariances up to lag 5 (the first element is the variance)
  
  X0 <- cbind(X0, apply( t( abs( apply( sign( apply(x,1,diff) ), 2, diff ) )), 1, sum)/2) #Number of turning points 
  
  X0
}