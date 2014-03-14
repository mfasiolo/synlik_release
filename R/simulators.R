############################################################################################################################################
##################################                    S I M U L A T O R S                       ############################################
############################################################################################################################################
# All simulators and wrappers around C/C++ simulators are here

########
# Blowfly full model (Simon's simulators) _INTERNAL_
########
## lu,lu1 are noise matrices of the same dimension.
## they contain ilogistic(uniform) random deviates,
## fixed from rep to rep to avoid monte carlo error. 
## ncol(lu) is number of reps, nrow(lu) is length of sim 
## including burn.in period, theta is parameter vector
## first transform the ilogistic(uniform) deviates in lu and lu1 to
## gamma deviates with mean 1 and given variance, and do it quickly!
.blowfly <- function (theta, lu, lu1, burn.in, pop.siz) 
{
  ## lu,lu1 are noise matrices of the same dimension.
  ## they contain logit(uniform) random deviates,
  ## fixed from rep to rep to avoid monte carlo error. 
  ## ncol(lu) is number of reps, nrow(lu) is length of sim 
  ## including burn.in period, theta is parameter vector
  
  ## first transform the logit(uniform) deviates in lu and lu1 to
  ## gamma deviates with mean 1 and given variance, and do it quickly!
  var <- theta[4];
  var1 <- theta[6]
  
  luk <- seq(min(lu),max(lu),length=200) ## even spacing on logit(prob) scale
  uk <- exp(luk);uk <- uk/(1+uk)          ## on prob scale
  gk <- qgamma(uk,shape=1/var,scale=var) ## reference quantiles
  gk[!is.finite(gk)] <- 0                ## can underflow badly at lower end
  e <- approx(x=luk,y=gk,xout=lu)$y      ## fast approximate version of qgamma(u,...) 
  
  luk <- seq(min(lu1),max(lu1),length=200)
  uk <- exp(luk);uk <- uk/(1+uk)
  gk <- qgamma(uk,shape=1/var1,scale=var1)
  gk[!is.finite(gk)] <- 0 
  e1 <- approx(x=luk,y=gk,xout=lu1)$y 
  
  n.t <- nrow(lu)-burn.in
  n.reps <- ncol(lu)
  n <- matrix(0,n.t,n.reps)
  oo <- .C("blowC", n=as.double(n), as.double(theta), as.double(e), as.double(e1), as.integer(burn.in), as.integer(n.t),
           as.integer(n.reps), PACKAGE="synlik")
  t(matrix(oo$n, n.t, n.reps))
}

################
#' Simulates from the Ricker map
#' 
#' @param nObs Length of each simulated time serie.
#' @param nsim Number of simulations from the model.
#' @param A vector of length 4 or a matrix of size nsim by 4.
#'        If param is a vector each of the nsim simulations will use the same parameters, if it's a matrix the i-th simulation
#'        will use the i-th row of param as parameters.                    
#' @param extraArgs Only needed for compatibility with the "synlik" package.                 
#' @param nBurn Number of initial steps to be discarded 
#'        before saving the following nObs steps.                     
#' @param randInit If TRUE the initial values of each of the nsim run will be simulated using runif(nsim, 0, 1),
#'                 if FALSE each of the nsim simulations will start from initVal.                    
#' @param initVal If randInit == FALSE each of the nsim simulations will start from this value.  
#' @param ... Only needed for compatibility with the "synlik" package.
#' @return A matrix of size nsim by nObs where each row is a trajectory simulated from the model.             
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>   
#' @export
#' 
rickerSimul <- function(param, nsim, extraArgs, randInit = TRUE, initVal = 1.0, ...)
{
  if( !all( c("nObs", "nBurn") %in% names(extraArgs) ) ) stop("extraArgs should contain nBurn and nObs")
  nBurn <- extraArgs$nBurn
  nObs <- extraArgs$nObs
  
  if( !is.matrix(param) ) param <- matrix(param, 1, length(param))
  
  .Call( "simpleModelsWrap", model = "ricker", days = nObs, nSimul = nsim, param = param, nBurn = nBurn, randInit = randInit, initVal = initVal, PACKAGE = "synlik" )
}


########
#' Simulates from the blowfly model
#' 
#' @param nObs Length of each simulated time serie.
#' @param nsim Number of simulations from the model.
#' @param A vector of length 2 or a matrix of size nsim by 2.
#'        If param is a vector each of the nsim simulations will use the same parameters, if it's a matrix the i-th simulation
#'        will use the i-th row of param as parameters. 
#' @param extraArgs Only needed for compatibility with the "synlik" package.                   
#' @param nBurn Number of initial steps to be discarded before saving the following nObs steps.               
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon Wood.
#' @export
blowSimul <- function(param, nsim, extraArgs, ...)
{
  if( !is.vector(param) ) param <- as.vector(param)
  if( length(param) != 6) stop("Wrong number of parameters")
  
  if( !all( c("nObs", "nBurn", "steps") %in% names(extraArgs) ) ) stop("extraArgs should contain nBurn, nObs and steps")
  nBurn <- extraArgs$nBurn
  nObs <- extraArgs$nObs
  steps <- extraArgs$steps
  
  totStep <- nBurn + nObs * steps
  
  noise1 <- matrix(ilogistic(runif(nsim * totStep)), totStep, nsim)
  noise2 <- matrix(ilogistic(runif(nsim * totStep)), totStep, nsim)
  
  simul <- .blowfly(theta = exp(param), 
                    lu = noise1, 
                    lu1 = noise2, 
                    burn.in = nBurn, pop.siz = nsim)[ , cumsum(c(1,rep(steps, nObs-1)))]
  
  return( simul )
}














