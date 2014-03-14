########
### Updates the parameters
########
## Does a Newton-Raphson update of the parameters
## maxJump is a scalar that is used to limit the size of the maximum move

.updateParams <- function(currPar, 
                          grad, 
                          currCov,
                          oldCov, 
                          gain = 1, 
                          constr = list(), 
                          quant = 0.1 )
{
  stopifnot( is.vector(currPar), is.vector(grad), is.matrix(currCov) )
  
  nPar <- length( grad )
  
  # Propose a jump and limit its length so that it doesn't go
  # behond a certain quantile of the local normal model
  delta <- drop( gain * currCov %*% grad )
  
  lim <- qchisq(quant, nPar)
  quad <- drop( crossprod(delta, qr.solve(oldCov, delta)) )
  delta <- delta * min( 1, sqrt(lim / quad) )
  
  newPar <- currPar - delta
  
  # Appling contraints to the new parameter vector
  if(length(constr))
  {
    newPar[constr$indexes] <- pmin(newPar[constr$indexes], constr$upper)
    newPar[constr$indexes] <- pmax(newPar[constr$indexes], constr$lower)
  }
  
  return(round(newPar, 14))
  
}



#######
#### Gets positive-definite covariace out of hessian 
#######
#
# ARGS:
# hessian = (matrix)
# uplim and uplim = (numeric) of length ncol(hessian) that limit the size of the diagonal
#                   elements of the covariance                    
# verbose = (logical) print whether the limits are violated

.getCovariance <- function(hessian, upLim, lowLim, verbose = FALSE) 
{
  stopifnot( is.matrix(hessian), is.vector(upLim), is.vector(lowLim), all(upLim > lowLim))
  
  nPar <- nrow(as.matrix(hessian))
  
  # Inverting Hessian to get the covariance
  covar <- .qrInverse(hessian, tolQR = 0, imposePD = TRUE, tilt = 1e-5)
  
  if( any(diag(covar) > upLim) || any(diag(covar) < lowLim) )
  {
    if(verbose == TRUE){ 
      print("Variance of proposal is too high or too low! (probably to high)") 
      print(diag(covar))
    }
    
    corr <- diag(sqrt(diag(covar))^-1, nPar)%*%covar%*%diag(sqrt(diag(covar))^-1, nPar)
    sdev <- sqrt(diag(covar))
    sdev <- pmax( pmin(sdev, sqrt(upLim)), sqrt(lowLim))
    covar  <- diag(sdev, nPar)%*%corr%*%diag(sdev, nPar)
  }
  
  round(covar, 14)
}

####

.mahaMovAv <- function(pos, Cov, oldPos, oldGrad, oldHess, redux = 1e-1, quant = 0.1)
{
  nPar <- length(pos)
  
  dist <- mahalanobis(oldPos, pos, cov = Cov * redux)
  w <- exp( -dist )
  w[ dist > qchisq(quant, nPar) ] <- 0.0
  w <- w / sum(w)

  grad <- colSums( oldGrad * w )

  hess <- Reduce("+", mapply("*", oldHess, w, SIMPLIFY = F) ) 
  
  list("grad" = grad, "hess" = hess)
}


#########
#### The stochastic optimization routine
#########
#' Synthetic likelihood maximization routine.
#'
#' @param object  ("synlik") object.
#' @param nIter   (integer) numer of iterations.
#' @param nsim    (integer) numer of simulations from the model at each step.
#' @param initCov  (matrix) initial covariance matrix used to simulate the paramters at each step.
#' @param initPar  (numeric) vector of initial values of the parameters.
#' @param addRegr  (logical) 
#'           if FALSE the statistics calculated by object@@summaries will be used (SL approach)
#'           if TRUE the simulated parameters will be regressed on the statistics and the 
#'           fitted values of the paramaters given the _observed_ statistics will be used as statistics
#'            (SL+ approach)
#' @param constr (named list) of 3 elements:
#'           [["indexes"]] = (numeric integers) indexes of the elements to check;
#'           [["upper"]]  = (numeric) upper bounds for the elements in "indexes";
#'           [["lower"]]  = (numeric) lower bounds for the elements in "indexes".
#' @param control Named list of control setting for the optimization routine.
#' @param multicore  (logical) if TRUE the object@@simulator and object@@summaries functions will
#'                    be executed in parallel. That is the nsim simulations will be divided in multiple cores.
#' @param ncores  (integer) number of cores to use if multicore == TRUE.
#' @param cluster an object of class c("SOCKcluster", "cluster"). This allowes the user to pass her own cluster,
#'                which will be used if multicore == TRUE. The user has to remember to stop the cluster. 
#' @param verbose  (logical) if TRUE lots of things will be printed.
#' @param ...  additional arguments to be passed to object@@simulator and object@@summaries.
#'             In general I would avoid using it and including in those two function everything they need.
#' @return object of class "synOptim": see info about that class.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>
#' @export
 
synMaxlik <- function(object, nIter, nsim, initCov, initPar = object@param, 
                       addRegr = TRUE, constr = list(), control = list(),
                       multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                       verbose = FALSE,  ...)
{
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
  
  # Reduce the object to "synlik" so that I avoid moving around all the additional slots of the "synMaxlik" class
  if( !class(object)[[1]] != "synlik" ) object <- as(object, "synlik")
  
  stopifnot(is.vector(initPar), is.matrix(initCov))
  
  nPar <- length(initPar)
  currCov <- oldCov <- initCov
  currPar <- oldPar <- initPar
  
  # Setting up control parameter
  ctrl <- list( "gain" = 1, 
                "gainExp" = 0.602, 
                "lag" = 1,
                "maxLag" = 30,
                "lagGrow" = 0.333,
                "limCov" = list("upper" = diag(currCov) * 10, "lower" = diag(currCov) * 0.1),
                "quant" = 0.9,
                "initGrad" = numeric(nPar),
                "initHess" = solve(currCov) )
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control)
  
  if(!identical(names(ctrl$limCov), c("upper", "lower"))) stop("limCov should have names \"upper\" and \"lower\" (in that order)")
  
  # Checking if the contraints are valid
  if( length(constr) )
    stopifnot( all(c("indexes", "upper", "lower") %in% names(constr)),
               length(constr$upper) == length(constr$lower), length(constr$lower) <= nPar, 
               length(constr$indexes) == length(constr$lower), all(constr$upper > constr$lower)
               ) 
    
  gainSeq <- ctrl$gain / ( (1:nIter) ^ ctrl$gainExp )
  
  gradient <- numeric(nPar)
  hessian_hat <- hessian_hat_hat <- matrix(NA, nPar, nPar)
  
  oldGradient <- ctrl$initGrad
  oldHessian <- ctrl$initHess
  
  # I will store the results here
  resultPar <- resultGrad <- matrix(NA, nIter, nPar)
  resultHess <- resultCovar <- list()
  resultLoglik <- numeric(nIter) 
  
  # Setting up a cluster if needed
  if(multicore) {
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik")
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
  }
  
  currLag <- ctrl$lag
  
  # The main loop of the optimization
  for(ii in 1:nIter)
  {
    
    currLag <- min(currLag + ctrl$lagGrow, ctrl$maxLag)
    
    # Calculate gradient and hessian                     
    tmp <- synGradient(object, param = currPar, nsim = nsim, covariance = currCov, addRegr, constr = constr, 
                       multicore = multicore, ncores = ncores, cluster = cluster, ...)  
    
    tmpGrad <- resultGrad[ii, ] <- -round(tmp$gradient, 14) 
    hessian_hat <- resultHess[[paste("Iter", ii, sep = "")]] <- -round(tmp$hessian, 14)
    resultLoglik[ii] <- tmp$llk
    
    # Gradient and Hessian are weighted averages of current and past values
    resultPar[ii, ] <- currPar
    tmp <- .mahaMovAv(pos = currPar, Cov = currCov, oldPos = resultPar[1:ii, , drop = F], 
                      oldGrad = resultGrad[1:ii, , drop = F], oldHess = resultHess[1:ii])
    gradient <- tmp$grad
    hessian_hat <- tmp$hess
#     gradient <- currLag/(1+currLag) * oldGradient + 1/(1+currLag) * tmpGrad
#     oldGradient  <- tmpGrad
#     
#     hessian_hat <- currLag/(1+currLag) * oldHessian + 1/(1+currLag) * hessian_hat
#     oldHessian <- hessian_hat
    
    # Get the covariance which will be used to simulate the parameters at the next iteration
    oldCov <- currCov
    currCov <- .getCovariance(hessian = hessian_hat, upLim = ctrl$limCov$upper, lowLim = ctrl$limCov$lower)
    
    oldPar <- currPar
    
    # Update paramters values
    currPar <- .updateParams(currPar = currPar, grad = gradient, 
                             currCov = currCov, oldCov = oldCov,
                             gain = gainSeq[ii], 
                             constr = constr)
        
    resultPar[ii, ] <- currPar
    resultCovar[[paste("Iter", ii, sep = "")]] <- currCov
    
    if(verbose) print( paste("Iteration", ii, ", parameters =", currPar ) )
    
  }
  
  if(multicore && clusterCreated) stopCluster(cluster)
  
  colnames(resultPar) <- names(object@param)
  for(ii in 1:nIter) dimnames(resultHess[[ii]]) <- dimnames(resultCovar[[ii]]) <- list(names(object@param),
                                                                                       names(object@param))
  
  # Setting up control list for "continue.synMaxlik" method
  toContinue <- ctrl
  toContinue$initGrad <- oldGradient
  toContinue$initHess <- oldHessian
  toContinue$lag <- currLag
  toContinue$gain <- tail(gainSeq, 1) 
  
  ### Class Definition
  .synMaxlik( object,
              initPar = initPar,
              nIter   = nIter,
              nsim    = nsim,
              initCov = initCov,
              addRegr  = addRegr,
              constr = constr,
              control = ctrl,
              continueCtrl = toContinue,
              multicore = multicore,
              ncores = ncores,
                                        
              resultPar = resultPar,
              resultGrad = resultGrad,
              resultHess = resultHess,
              resultCovar = resultCovar,
              resultLoglik = resultLoglik)
  
}


############################
##### Function to estimate gradient and hessian of synthetic likelihood
############################
# This is Simon's local regression to get gradient and hessian with one modification:
# the second regression (the one involving the residuals) is linear not quadratic. This is suggested 
# by the article "Local Polynomial Variance-Function Estimation" since the variance shouldn't need non linear
# terms locally. 
# 
#' Estimate gradient and hessian of the synthetic likelihood.
#' 
#' @param object an object of class "synlik".
#' @param param (numeric) vector containing the current values of the model's parameters.
#' @param nsim   (integer) number of simulated statistics.
#' @param covariance (matrix) covariance matrix used to simulate the parameters.
#' @param addRegr (logical). 
#'           If FALSE the statistics calculated by object@@summaries will be used (SL approach);
#'           if TRUE the simulated parameters will be regressed on the statistics and the 
#'           fitted values of the paramaters given the _observed_ statistics will be used as statistics
#'           (referred to as SL+ approach).
#' @param constr (named list) used to impose constraints on the parameters.
#'           Composed of 3 elements:
#'           [["indexes"]] = (numeric integers) indexes of the elements to check;
#'           [["upper"]]   = (numeric) upper bounds for the elements in "indexes";
#'           [["lower"]]   = (numeric) lower bounds for the elements in "indexes".
#' @param multicore  (logical) if TRUE the object@@simulator and object@@summaries functions will
#'                    be executed in parallel. That is the nsim simulations will be divided in multiple cores.
#' @param ncores  (integer) number of cores to use if multicore == TRUE.
#' @param cluster an object of class c("SOCKcluster", "cluster"). This allowes the user to pass her own cluster,
#'                which will be used if multicore == TRUE. The user has to remember to stop the cluster. 
#' @param ... additional arguments to be passed to object@@simulator and object@@summaries.
#'            In general I would avoid using it and including in those two function everything they need.
#'
#' @return  a list containing: ["gradient"] = (numeric) estimated gradient of the synthetic log-likelihood at currPar;
#'                          ["hessian"]  = (matrix) estimated hessian of the synthetic log-likelihood at currPar;  
#'                          ["llk"]      = (scalar) estimated value of the synthetic log-likelihood at currPar.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>                         
#' @export
#' 
synGradient <- cmpfun(function(object, param, nsim, covariance, addRegr = TRUE, constr = list(), 
                               multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, 
                               tolVar = 10 * .Machine$double.eps, ...)
{
  if(!is(object, "synlik")) stop("object has to be of class \"synlik\" ")
  
  # Reduce the object to "synlik" so that I avoid moving around all the additional slots of the "synMaxlik" class (for example)
  if( !class(object)[[1]] != "synlik" ) object <- as(object, "synlik")
  
  theData <- object@data
  
  nPar <- length(param)
  
  param <- unname(param)
  covariance <- unname(covariance)
  
  if(multicore)
  { 
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik") 
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    
    coresSchedule <- c( rep(floor(nsim / ncores), ncores - 1), floor(nsim / ncores) + nsim %% ncores)
    
    tmplist <- clusterApply( cluster, 
                             coresSchedule, 
                             function(input, ...)
                             {
                               # Simulate parameters and data
                               simulParams <- .paramsSimulator(theMean = param, covar = covariance, nsim = input, constr = constr)
                               
                               simulStats <- simulate.synlik(object, param = simulParams, nsim = input, stats = TRUE, clean = FALSE, verbose = FALSE, ...)
                               
                               return( list("simulParams" = simulParams, "simulStats" = simulStats) )
                             } 
                             , ... )
    
    simulParams <- do.call(rbind, lapply( tmplist, function(input) input$simulParams ) )
    simulStats  <- do.call(rbind, lapply( tmplist, function(input) input$simulStats ) )
    rm(tmplist)
    
    if(clusterCreated) cluster <- stopCluster(cluster)
    
  } else{
    
    # Simulate parameters and data
    simulParams <- .paramsSimulator(theMean = param, covar = covariance, nsim = nsim, constr = constr)
    simulStats  <- simulate.synlik(object, param = simulParams, nsim = nsim, stats = TRUE, clean = FALSE, verbose = FALSE, ...)  
    
  }
  
  clean <- .clean(X = simulStats, verbose = TRUE)
  if(clean$nBanned > 0){
    simulStats <- clean$cleanX
    simulParams <- simulParams[-clean$banned, ]
  }
  rm(clean)
  
  nGood <- nrow(simulStats)

  # if addRegr == TRUE we use SL+ and we regress parameters on statistics 
  if(addRegr){ 
    simulStats <- cbind(1, simulStats)
    nStats <- nPar
    
    qrx <- qr(simulStats, tol = 0)
    # matrix of linear regs coeffs for all summary statistics
    fearn_beta <- t( qr.coef(qrx, simulParams) )
    fearn_beta[is.na(fearn_beta)] <- 0 # Putting to zero the regression coefficients that are NAs
    
    simulStats <- tcrossprod(fearn_beta, simulStats)
    
    summaries <- object@summaries
    obserStats <- if( !is.null(summaries) ) drop( summaries(x = theData, extraArgs = object@extraArgs, ...) ) else drop( theData )
    obserStats <- fearn_beta %*% c(1, obserStats)
  } else{
    simulStats <- t(simulStats)
    nStats <- nrow(simulStats)
    
    summaries <- object@summaries
    obserStats <- if( !is.null(summaries) ) drop( summaries(x = theData, extraArgs = object@extraArgs, ...) ) else drop( theData )
  }
  
  ##### Regression of summary statistic (simulStats) on the parameters of each replicate (simulParams)
  
  X1 <- simulParams
  X1 <- sweep(X1, 2, param)    #Center around current position
  
  tmpModel <- cbind(1, X1, (X1^2)/2)
  
  if(nPar > 1){
    comb <- t( combn(nPar, 2) )
    for(jj in 1:nrow(comb)){tmpModel <- cbind(tmpModel, X1[ , comb[jj, 1]] * X1[ , comb[jj, 2]])}
  }
  
  X1 <- tmpModel;              # Model matrix for first regression                             
  X2 <- X1[ , 1:(nPar+1)]      # Model matrix for second regression
  rm(tmpModel)
  
  qrx1 <- qr(X1, tol = 0)
  qrx2 <- qr(X2, tol = 0) 
  
  # Calculate the betas of the 1st regression
  beta <- t( qr.coef(qrx1, t(simulStats)) )
  beta[is.na(beta)] <- 0
  
  ########################################################
  ######### Regressing entries of the covariance matrix on parameters 
  ########################################################
  
  # WARNING: it is better to simulated very close to the current value of parameters,
  # otherwise the quadratic model is biased, we get mean(res) != 0 and SIGMA_HAT is not invertible!
  
  # Obtaining residuals (res) in the two ways previously described.
  res <- matrix(NA, nStats, nGood)
  for(ii in 1:nStats) res[ii, ] <- simulStats[ii, ] - tcrossprod(beta[ii, ], X1)
  
  # Down-weighting the most extreme residuals
  fixCov <- robCov(res)
  res <- res * fixCov$weights
  
  # Array where the residual matrices for each path are stack one over the other.
  Dvett <- array(NA, c(nStats, nStats, nGood))
  for(ii in 1:nGood) Dvett[ , , ii] <- tcrossprod(res[ , ii], res[ , ii])
  
  # Matrix of expected SIGMA, first derivative of sigma wrt the parameters
  # and second derivatives. (Example: SIGMA12 = (D^2 SIGMA)/(D theta1 D theta2))     
  SIGMA_HAT <- matrix(0, nStats, nStats)
  firstSIG <- array(NA, c(nStats, nStats, nPar))
  
  # Now I regress each vector D_vett[,,i] on the parameters (using the same LINEAR regression)
  # The intercepts are the expected elements of the covariance matrix. I move along the upper triangular
  # part of the matrix. 
  # N.B. Since this is a LINEAR regression this code is not the same as Simon's notes. 
  
  for (ff in 1:nStats) for(ii in ff:nStats){ 
    tmp <- t(qr.coef(qrx2, as.matrix(Dvett[ff, ii, ])))
    tmp[is.na(tmp)] <- 0
    SIGMA_HAT[ff, ii] <- SIGMA_HAT[ii, ff] <- tmp[1]
    firstSIG[ff, ii, ] <- firstSIG[ii, ff, ] <- tmp[2:(nPar+1)]
  }
  
  # Identify statistics with very low variance (diagonal elements of SIGMA_HAT)
  # and removing them and their coefficients
  lowVar <- which( diag(SIGMA_HAT) < tolVar )
  if( length(lowVar) ) {
    if( addRegr ) stop( paste("The variance of one of the parameters is  <", tolVar, ""))
    nStats <- nStats - length(lowVar)
    if(nStats == 0) stop( paste("All the chosen statistics have variance <", tolVar) )
    beta <- beta[-lowVar, ]
    obserStats <- obserStats[-lowVar]
    firstSIG <- firstSIG[-lowVar, -lowVar, ] 
    SIGMA_HAT <- SIGMA_HAT[-lowVar, -lowVar] 
    warning( paste("There are", length(lowVar), "statistics with variance <", tolVar, "they were removed.") )
  }
  
  # Get the (scaled) QR decomposition of SIGMA_HAT, that later I'll use to 
  # solve several linear systems (so I avoid inverting it).
  D <- diag(diag(SIGMA_HAT)^-0.5, nStats)

  sigQR <- qr( D%*%SIGMA_HAT%*%D, tol = 0 )
  
  #############################################################################
  ###### Calculate gradient and Hessian of the log-likelhood wrt the parameters
  #############################################################################
  
  ######### Calculating the gradient of the synthetic likelihood
  # DmuDth1 = the first derivatives of mu wrt the parameters 
  # DmuDth2 = the second derivatives of mu wrt the parameters
  # currStat = the expected value of the statistic at the current position
  # obsRes = statistics(observed_path) - curr_stat
  # sigByRes = SIGMA_HAT^-1 %*% (y - currStat)
  DmuDth1 <- beta[ , 2:(nPar+1), drop = FALSE]
  DmuDth2 <- beta[ , (nPar+2):ncol(X1), drop = FALSE]
  currStat <- beta[ , 1]
  obsRes <- drop( obserStats - currStat )
  sigByRes <- drop(D %*% qr.solve(sigQR, D %*% obsRes, tol = 0))
  
  firstDeriv <- numeric(nPar)
  for(kk in 1:nPar)
  {
    firstDeriv[kk] <- 0.5 * (2 * crossprod(DmuDth1[ , kk], sigByRes)    +     
                               
                               crossprod(obsRes, D %*% qr.solve(sigQR, D%*%firstSIG[ , , kk]%*%sigByRes, tol = 0)) -
                               
                               eSaddle:::.Trace(D %*% qr.solve(sigQR, D %*% firstSIG[ , , kk], tol = 0)) )
  }
  
  
  ######### Calculating the Hessian of the synthetic likelihood
  # Different from Simon's document because the second regression is linear, and hence
  # all the second derivatives of the covariance matrix wrt the paramets have been put to zero.
  # I fill the hessian moving on the upper triangle and I use the switch to select the correct
  # matrix of the first derivatives of the covariance wrt the parameters. 
  
  secondDeriv <- matrix(0, nPar, nPar)
  
  if(nPar == 1){
    
    indexes <- matrix(1, 1, 1)
    
  }else{
    
    # Create matrix of indexes to manage the second derivarives (stored in DmuDth2)
    indexes <- diag(seq(1:nPar))
    entries <- seq(nPar + 1, nPar + factorial(nPar)/(factorial(2)*factorial(nPar-2)))
    zz <- 1
    for(jj in 1:nPar){
      indexes[jj, -(1:jj)] <- entries[zz:(zz + nPar - jj - 1)]
      zz <- zz + nPar - jj
    }
  }
  
  
  for(kk in 1:nPar) 
    for(ff in kk:nPar)
    {
      zz <- indexes[kk, ff]
      
      DsigDthK <- firstSIG[ , , kk]
      DsigDthL <- firstSIG[ , , ff]
      
      secondDeriv[kk, ff] <- secondDeriv[ff, kk] <-
        0.5 * (  
          2 * crossprod(DmuDth2[ , zz], sigByRes) - 
            2 * crossprod(DmuDth1[ , kk], D %*% qr.solve(sigQR, D%*%(DsigDthL%*%sigByRes), tol = 0)) -
            2 * crossprod(DmuDth1[ , kk], D %*% qr.solve(sigQR, D%*%DmuDth1[ , ff], tol = 0)) - 
            2 * crossprod(DmuDth1[ , ff], D %*% qr.solve(sigQR, D%*%(DsigDthK%*%sigByRes), tol = 0)) -  
            2 * crossprod(sigByRes, DsigDthL%*%(D %*% qr.solve(sigQR, D%*%(DsigDthK%*%sigByRes), tol = 0)) ) +
            eSaddle:::.Trace(D %*% qr.solve(sigQR, D%*%DsigDthL, tol = 0) %*% (D %*% qr.solve(sigQR, D%*%DsigDthK, tol = 0)) )
        )  
    }
  
  llk <- - 0.5 * crossprod(obsRes, sigByRes) - 0.5 * log( abs( prod( diag(qr.R(sigQR)) / diag(D^2) ) ) )
  
  list("gradient" = firstDeriv, "hessian" = secondDeriv, "llk" = llk)
  
})