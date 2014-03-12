#########
### Simulator
#########

budSimul <- function(param, nsim, extraArgs, trans = FALSE, ...)
{
  if( is.null(extraArgs$pompModel)  ) stop("extraArgs$pompModel is missing")
  if( !is.loaded("pompExamples") ) library(pompExamples)
  
  pompModel <- extraArgs$pompModel
  
  names(param) <- names(pompModel@params)
  
  if( trans )
  {
    param[pompModel@userdata$logitvar] <- logistic(param[pompModel@userdata$logitvar])
    param[pompModel@userdata$logvar] <-   exp(param[pompModel@userdata$logvar])
  } 
  
  tmp <- simulate(pompModel, params = param, nsim = nsim)
  
  if(nsim > 1)
  {
    Q <- t( sapply(tmp, function(input) input@data[1, ]) )
    N <- t( sapply(tmp, function(input) input@data[2, ]) )
    S <- t( sapply(tmp, function(input) input@data[3, ]) )
  } else {
    Q <- tmp@data[1, ]
    N <- tmp@data[2, ]
    S <- tmp@data[3, ] 
  }
  
  return(list("Qobs" = Q, "Nobs" = N, "Sobs" = S))
}



################################
############ Statistics
################################

# Calculates variance of increments: Var( X_t - X_(t-lag) )
# X is a matrix, where each row is a time series.
.lagDistr <- function(X, lag)
{
  n <- ncol(X)
  tmp <- X[ , 1:(n - lag)] - X[ , -(1:lag)]
  
  return( rowMeans(tmp^2) / ncol(X) )
}

budStats <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  # Extract data
  obsQ <- drop(extraArgs$obsData[["Qobs"]])
  obsBud <- drop(extraArgs$obsData[["Nobs"]])
  obsP <- drop(extraArgs$obsData[["Sobs"]])
  stopifnot(length(obsP) != 0, length(obsQ) != 0, length(obsBud) != 0)
    
  # Extract simulations
  Q <- x[[1]]
  if (!is.matrix(Q)) Q <- matrix(Q, 1, length(Q))
  tQ <- t(Q)
  
  Bud <- x[[2]]
  if (!is.matrix(Bud)) Bud <- matrix(Bud, 1, length(Bud))
  tBud <- t(Bud)
  
  P <- x[[3]]
  if (!is.matrix(P)) P <- matrix(P, 1, length(P))
  tP <- t(P)
  
  nsim <- nrow(Q)
  nObs <- ncol(Bud)
  
  # Constructing lagged simulations
  laggedB <- log( Bud[ , 2:(nObs-1), drop = FALSE] )
  laggedQ <- log( Q  [ , 2:(nObs-1), drop = FALSE] )
  laggedP <- log( P  [ , 2:(nObs-1), drop = FALSE] )
  
  laggedB2 <- log( Bud[ , 1:(nObs-2), drop = FALSE] )
  laggedQ2 <- log( Q  [ , 1:(nObs-2), drop = FALSE] )
  laggedP2 <- log( P  [ , 1:(nObs-2), drop = FALSE] )
    
  #####
  # Quality
  #####
  X0 <- t(order.dist(tQ, obsQ, np=3, diff=1))        
  #X0 <- cbind(X0, t( nlar(tQ, lag=c(1,3), power=c(1, 1, 2)) ))        
  ### X0 <- cbind(X0, t(sl.acf(tQ, max.lag=5)))
  # X0 <- apply( t( abs( apply( sign( apply(Q, 1, diff) ), 2, diff ) )), 1, sum)/2
  
  # Regressing Quality on lagged Budmoth
  tmp <- matrix(NA, nsim, 7)
  for(ii in 1:nsim)
  {  
     modelMat <- cbind(1, laggedB[ii, ], laggedB[ii, ]^2, laggedB[ii, ]^3, laggedB2[ii, ], laggedB2[ii, ]^2,
                          laggedQ[ii, ])
     fit <- fastLmPure(modelMat, log(Q[ii, -(1:2)]) )
     tmp[ii, 1:6] <- fit$coef[-1]
     tmp[ii, 7]   <- fit$stderr[1]
  }
  
  X0 <- cbind(X0, tmp)
    
  ######
  # Budmoth
  ######
  X0 <- cbind(X0, t(order.dist(tBud, obsBud, np=3, diff=1)))        
  #X0 <- cbind(X0, t(nlar(log(tBud), lag=c(1, 4), power=c(1, 1))))
  #X0 <- cbind(X0, rowMeans(Bud))         
  #X0 <- cbind(X0, t(sl.acf(tBud, max.lag=5))[ , c(1, 3, 5)]) 
  
  # Regressing Bud on lagged Quality and lagged Parasite 
  # Regressing Quality on lagged Budmoth
  tmp <- matrix(NA, nsim, 11)
  for(ii in 1:nsim)
  {  
    modelMat <- cbind(1, laggedB[ii, ], laggedB2[ii, ], laggedB2[ii, ]^2, # laggedB3[ii, ], laggedB3[ii, ]^2
                         laggedQ[ii, ], laggedQ2[ii, ],
                         laggedP[ii, ], laggedP[ii, ]^2, laggedP2[ii, ], laggedP2[ii, ]^2)
    fit <- fastLmPure(modelMat, log(Bud[ii, -(1:2)]) )
    tmp[ii, 1:10] <- fit$coef
    tmp[ii, 11]   <- fit$stderr[1]
  }
  
  X0 <- cbind(X0, tmp)
  
  ######
  # Parasite
  ######

  X0 <- cbind(X0, t(order.dist(tP, obsP, np=3, diff=1)))        
  #X0 <- cbind(X0, t(nlar(log(tP), lag=c(1), power=c(1)))) 
  #X0 <- cbind(X0, rowMeans(P)) 
  #X0 <- cbind(X0, t(sl.acf(tP, max.lag=3)))
  #X0 <- cbind(X0, apply( t( abs( apply( sign( apply(P, 1, diff) ), 2, diff ) )), 1, sum)/2)
  
  # Regressing parasite on lagged bud    
  tmp <- matrix(NA, nsim, 5)
  for(ii in 1:nsim)
  {  
    modelMat <- cbind(1, laggedB[ii, ], laggedB[ii, ]^2, # laggedB2[ii, ]^2,
                         laggedP[ii, ])
    fit <- fastLmPure(modelMat, log(P[ii, -(1:2)]) )
    tmp[ii, 1:4] <- fit$coef
    tmp[ii, 5]   <- fit$stderr[1]
  }
  
  X0 <- cbind(X0, tmp)
  
  X0
}
























