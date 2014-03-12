############################################################################################################################################
##################################                    S I M U L A T O R S                       ############################################
############################################################################################################################################
# All simulators and wrappers around C/C++ simulators are here


############################################################################################################################################
############################
### Simon's simulators
############################

########
# Blowfly full model
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

########
# Bupalus
########

bup.para <- function(theta,e,u=0,burn.in) {
  ## e is noise matrix. ncol(e) is number of reps, nrow(e) is length of sim 
  ## including burn.in period, theta is parameter vector
  ## u is a sequence of uniform deviates, used for generating observed
  ## from expected while avoiding Monte Carlo error.
  
  n.t <- nrow(e)-burn.in
  n.reps <- ncol(e)
  n <- matrix(0,n.t,n.reps)
  oo <- .C("bup_par",n=as.double(n),as.double(theta),as.double(e),as.integer(burn.in),as.integer(n.t),
           as.integer(n.reps),PACKAGE="synlik")
  
  matrix(oo$n^.25+qnorm(u,0)*theta[6],n.t,n.reps) ## transformed response with obs error
}



############################################################################################################################################
###############################
### Wrappers for Rcpp simple model simulators
###############################


# Maynard Smith
#' Simulates from the Maynard-Smith map
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
maynardSimul <- function(param, nsim, extraArgs, randInit = TRUE, initVal = 1.0, ...)
{
  if( !all( c("nObs", "nBurn") %in% names(extraArgs) ) ) stop("extraArgs should contain nBurn and nObs")
  nBurn <- extraArgs$nBurn
  nObs <- extraArgs$nObs

  if(!is.matrix(param)) param <- matrix(param, 1, length(param)) 
  
  .Call( "simpleModelsWrap", model = "maynard", days = nObs, nSimul = nsim, param = param, nBurn = nBurn, randInit = randInit, initVal = initVal, PACKAGE = "synlik" )
}

# Ricker
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

# Varley
#' Simulates from the Varley map
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
varleySimul <- function(param, nsim, extraArgs, randInit = TRUE, initVal = 1.0, ...)
{
  if( !all( c("nObs", "nBurn") %in% names(extraArgs) ) ) stop("extraArgs should contain nBurn and nObs")
  nBurn <- extraArgs$nBurn
  nObs <- extraArgs$nObs
  
  if( !is.matrix(param) ) param <- matrix(param, 1, length(param))
  
  .Call( "simpleModelsWrap", model = "varley", days = nObs, nSimul = nsim, param = param, nBurn = nBurn, randInit = randInit, initVal = initVal, PACKAGE = "synlik" )
}

# Hassell
#' Simulates from the Hassell map
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
hassellSimul <- function(param, nsim, extraArgs, randInit = TRUE, initVal = 1.0, ...)
{
  if( !all( c("nObs", "nBurn") %in% names(extraArgs) ) ) stop("extraArgs should contain nBurn and nObs")
  nBurn <- extraArgs$nBurn
  nObs <- extraArgs$nObs 
  
  if( !is.matrix(param) ) param <- matrix(param, 1, length(param))
  
  .Call( "simpleModelsWrap", model = "hassell", days = nObs, nSimul = nsim, param = param, nBurn = nBurn, randInit = randInit, initVal = initVal, PACKAGE = "synlik" )
}

# Generalized Ricker
#' Simulates from the Generalized Ricker map
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
genRickerSimul <- function(param, nsim, extraArgs, randInit = TRUE, initVal = 1.0, ...)
{
  if( !all( c("nObs", "nBurn") %in% names(extraArgs) ) ) stop("extraArgs should contain nBurn and nObs")
  nBurn <- extraArgs$nBurn
  nObs <- extraArgs$nObs 
  
  if( !is.matrix(param) ) param <- matrix(param, 1, length(param))
  
  .Call( "simpleModelsWrap", model = "genRicker", days = nObs, nSimul = nsim, param = param, nBurn = nBurn, randInit = randInit, initVal = initVal, PACKAGE = "synlik" )
}

# Pennycuik
#' Simulates from the Pennycuik map
#' 
#' @param nObs Length of each simulated time serie.
#' @param nsim Number of simulations from the model.
#' @param A vector of length 4 or a matrix of size nsim by 4.
#'        If param is a vector each of the nsim simulations will use the same parameters, if it's a matrix the i-th simulation
#'        will use the i-th row of param as parameters.                    
#' @param extraArgs List containg [["nBurn"]] which is the number of initial steps to be discarded 
#'        before saving the following nObs steps.                     
#' @param randInit If TRUE the initial values of each of the nsim run will be simulated using runif(nsim, 0, 1),
#'                 if FALSE each of the nsim simulations will start from initVal.                    
#' @param initVal If randInit == FALSE each of the nsim simulations will start from this value.  
#' @param ... Only needed for compatibility with the "synlik" package.
#' @return A matrix of size nsim by nObs where each row is a trajectory simulated from the model.             
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>   
#' @export
#' 
pennySimul <- function(param, nsim, extraArgs, randInit = TRUE, initVal = 1.0, ...)
{
  if( !all( c("nObs", "nBurn") %in% names(extraArgs) ) ) stop("extraArgs should contain nBurn and nObs")
  nBurn <- extraArgs$nBurn
  nObs <- extraArgs$nObs 
  
  if( !is.matrix(param) ) param <- matrix(param, 1, length(param))
  
  .Call( "simpleModelsWrap", model = "penny", days = nObs, nSimul = nsim, param = param, nBurn = nBurn, randInit = randInit, initVal = initVal, PACKAGE = "synlik" )
}






############################################################################################################################################

###############
#### Simulator for stable distribution
###############
#' Simulates from the alpha stable family of distribution
#' 
#' @param nObs Length of each simulated time serie.
#' @param nsim Number of simulations from the model.
#' @param A vector of length 4 or a matrix of size nsim by 4.
#'        If param is a vector each of the nsim simulations will use the same parameters, if it's a matrix the i-th simulation
#'        will use the i-th row of param as parameters.
#' @param extraArgs Only needed for compatibility with the "synlik" package.
#' @param ... Only needed for compatibility with the "synlik" package.
#' @return A matrix of size nsim by nObs where each row is a sample from the distribution.                                  
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>   
#' @export
#' 
stableSimul <- function(param, nsim, extraArgs, ...)
{
  if( !is.loaded("stabledist") ) { library(stabledist) } 
  
  if( !( c("nObs") %in% names(extraArgs) ) ) stop("extraArgs should contain nObs")
  nObs <- extraArgs$nObs 
  
  output <- matrix(NA, nsim, nObs)
  
  if(is.vector(param))
  {
    param[c(1, 3)] <- exp(param[c(1, 3)])
    if(abs(param[1] - 1) < 0.01) stop("alpha == 1 is not allowed")
    stopifnot( length(param) == 4 )
    for(ii in 1:nsim) output[ii, ] <- rstable(nObs, alpha = param[1], beta = param[2], 
                                              gamma = param[3], delta = param[4])
    return(output)
  } 
  
  if(is.matrix(param))
  {
    param[ , c(1, 3)] <- exp(param[ , c(1, 3)])
    stopifnot( nrow(param) == nsim, ncol(param) == 4  )
    for(ii in 1:nsim) output[ii, ] <- rstable(nObs, alpha = param[ii, 1], beta = param[ii, 2], 
                                              gamma = param[ii, 3], delta = param[ii, 4])
    return(output)
  }
  
  stop("param has to be either a matrix of a vector")
}





#######
### Simulator for blowfly model
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





#######
### Simulator for Lokta-volterra model
########
#' Simulates from the Lokta-volterra model
#' 
#' @param nObs Length of each simulated time serie.
#' @param nsim Number of simulations from the model.
#' @param A vector of length 2 or a matrix of size nsim by 2.
#'        If param is a vector each of the nsim simulations will use the same parameters, if it's a matrix the i-th simulation
#'        will use the i-th row of param as parameters. 
#' @param extraArgs Only needed for compatibility with the "synlik" package.                   
#' @param nBurn Number of initial steps to be discarded before saving the following nObs steps.  
#' @param totTime Total period of time to which the simulations correspond (ex: 15 years).
#' @param nSteps Number of finite differences steps between 2 observations.
#' @param randInit random initialization for prey and predator?
#' @param x0 Initial state for the prey.
#' @param y0 Initial state for the predator.
#' @param onlyPrey If TRUE only the prey simulations are returned.                  
#' @param ... Only needed for compatibility with the "synlik" package. 
#' 
#' @return If onlyPrey == FALSE it returns a list where [["prey"]] and [["predator"]] are
#'         two matrices of size nsim by nObs, if onlyPrey == TRUE it return a matrix containing
#'         only the simulations for the prey.                    
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>   
#' @export
#' 
loktaSimul <- function(param, nsim, extraArgs, nBurn, 
                       totTime = 15, nSteps = 5,  
                       randInit = TRUE, x0 = 1, y0 = 0.5, onlyPrey = TRUE, ...)
{
  if(!is.matrix(param)) param <- matrix(param, 1, length(param))
  
  if( !( c("nObs") %in% names(extraArgs) ) ) stop("extraArgs should contain nObs")
  nObs <- extraArgs$nObs 
  
  tmp <- .Call("loktaCpp", nObs_ = nObs, nsim_ = nsim, param_ = param, 
               nBurn_ = nBurn, totTime_ = totTime, nSteps_ = nSteps, x0_ = x0, 
               y0_ = y0, randInit_ = randInit, PACKAGE = "synlik")
  
  if(onlyPrey) return(tmp$prey) else return(tmp)
  
}

# TEST
# dati <- loktaSimul(nObs = 100, nsim = 1000, param = matrix(log(c(1, 1, 0.1, 0.3)), 1, 4), 
#                    nBurn = 500, totTime = 15, nSteps = 5, x0 = 1, y0 = 0.5, randInit = TRUE, onlyPrey = FALSE)
# 
# ii <- 1
# plot(dati$prey[ii, ], type = 'l', ylim = c(0, 5))
# lines(dati$predator[ii, ], col = 2)
# ii = ii + 1

########
## Trait model
########
# Simulates from the trait model of "A stochastic dispersal-limited trait-based model of community dynamics" (2010)
# NB simulates statistics directly, so you don't need to specify object@@summaries
traitSimul <- function(param, nsim, extraArgs, ...)
{ 
  if( !is.loaded("EasyABC") ) library("EasyABC") 
  
  if( !is.matrix(param) ) param <- matrix(param, 1, length(param))
  if(ncol(param) != 4) stop("Wrong number of parameters")
  stopifnot(all(param[3] < 125 || param[3] > -25))
  
  seeds <- sample(1:10000, nsim, replace = FALSE)
  output <- matrix(NA, nsim, 4)
  
  nt <- 1
  J <- 500
  
  for(ii in 1:nsim) 
    output[ii, ] <- trait_model(input=c(seeds[ii], J, exp(param[ii, 1:2]), nt, param[ii, 3], exp(param[ii, 4])))
  
  return(output)
}




############################################################################################################################################

###############
#### Simulators for voles model (Turchin and Ellner (2000): LIVING ON THE EDGE OF CHAOS: POPULATION DYNAMICS OF FENNOSCANDIAN VOLES)
###############
##### Basic simulator 
#' Simulates from the voles models of Turchin and Ellner (2000).
#' 
#' @param nObs Length of each simulated time series in months.
#' @param nsim Number of simulations from the model.
#' @param param If model == "full" it must be a vector of length 9 or a matrix of size nsim by 9, if model == "standard" 
#'        it must be a vector of length 7 or a matrix of size nsim by 7.
#'        If param is a vector each of the nsim simulations will use the same parameters, if it's a matrix the i-th simulation
#'        will use the i-th row of param as parameters.
#' @param extraArgs Only needed for compatibility with the "synlik" package. 
#' @param model If model == "full" the function will simulate from the full model described in eq. 2 of Turchin and Ellner (2000),
#'              model == "standard" the function will simulate from the full model described in eq. 6 of Turchin and Ellner (2000).                   
#' @param nBurn Number of initial steps to be discarded before saving the following nObs steps.  
#' @param nSteps Number of finite differences steps between 2 observations (an interval which corresponds to a month).                
#' @param ... Only needed for compatibility with the "synlik" package. 
#' 
#' @return A list where [[""voles"]] and [["weasels"]] are two matrices of size nsim by nObs where each row is a trajectory
#'         simulated from the model.
#' @references Turchin and Ellner (2000), LIVING ON THE EDGE OF CHAOS: POPULATION DYNAMICS OF FENNOSCANDIAN VOLES.                    
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>   
#' @export volesSimulator
#' 
volesSimulator <- function(param, nsim, nObs, extraArgs = NULL, model = "full", nBurn = 12 * 10, nSteps = 10, T0 = 0, ...)
{
  if( !is.matrix(param) ) param <- matrix(param, 1, length(param))
  
  if( !is.null(extraArgs) ){ 
    if( !is.null(extraArgs$model) ) model <- extraArgs$model
    if( !is.null(extraArgs$nBurn) ) nBurn <- extraArgs$nBurn
    if( !is.null(extraArgs$T0) ) T0 <- extraArgs$T0
  }
  
  if( nBurn %% 12 != 0 ) warning("nBurn % 12 != 0 but nBurn must be a multiple of 12, 
                                  otherwise the first observed month doesn't correspond to January and so on!")
  
#   tmp <- switch(model,
#                 full = .Call("volesFullCpp", nMon_ = nObs, nSimul_ = nsim, nBurn_ = nBurn, param_ = param, nSteps_ = nSteps, T0_ = T0,
#                              randInit_ = TRUE, startVole_ = 0.0,  startWeas_ = 0.0, addObsNoise_ = TRUE, PACKAGE = "synlik"),
#                 standard = .Call("volesStdCpp", nMon_ = nObs, nSimul_ = nsim, nBurn_ = nBurn, param_ = param, nSteps_ = nSteps, T0_ = T0,
#                                  randInit_ = TRUE, startVole_ = 0.0,  startWeas_ = 0.0, addObsNoise_ = TRUE, PACKAGE = "synlik"),  
#                 stop("Wrong model name")
#   )
  
  if(!exists("volesFullCpp")) source("~/Desktop/All/Dropbox/Work/Reports/1st_Article_Review/Code/Classification/pomp_voles_model.R")
  tmp <- volesFullCpp(nMon = nObs, nSimul = nsim, nBurn = nBurn, param = param, nSteps = nSteps, T0 = T0,
                      randInit = TRUE, startVole = 0.0,  startWeas = 0.0, addObsNoise = TRUE)
  
  return(tmp)
}


### Wrappers used by synlik 
volesWrap <- function(param, nsim, extraArgs, nObs, ...)
{
  if( !is.loaded("synlik") ) library("synlik")
  nObs <- extraArgs$nObs
  stopifnot( !is.null(nObs) )
  
  simul <- volesSimulator(param = param, nsim = nsim, extraArgs = extraArgs, nObs = nObs, ...)[["voles"]]
  
  # Save only the months with indexes in extraArgs$monthsToSave
  if( !is.null(extraArgs$monthsToSave) ){ 
    return( simul[ , extraArgs$monthsToSave, drop = FALSE] )
  } else {
    return( simul )
  }

}


############################################################################################################################################

#######
#### Simulator of the Malaria Model of Bhadhra et at. 2011.
#######
#' Simulator for the S2EI2 malaria model of Bhadhra et at. 2011.
#'
#' @param procparam: matrix of size (nSimul by 10), the i-th row is the set of process parameters 
#'                    used for the i-th simulation. The columns represent, in order, the _LOG_ of parameters
#'                    mu_EI1, mu_I1I2, mu_I2S2, mu_S2S1, mu_I1S1, q, covarCoeff, sigma, c and tau.  
#' @param obsparam: matrix of size (nSimul by 2), the i-th row is the set of parameters of the observational 
#'                   process used for the i-th simulation. The columns represent, in order, the _LOG_ of parameters
#'                   rho and psiSquare.                  
#' @param splineparam: matrix of size (nSimul by splineSettings[["nBasis"]]), the i-th row is the set of spline coefficients 
#'                      used for the i-th simulation.
#' @param initStates:  matrix of size (nSimul by 7), the i-th row is the set of initial values 
#'                     used for the i-th simulation. The columns represent, in order, the states
#'                     S1, S2, E, I1, I2, lambda1 and lambda2.
#'                     NB: S1, S2, E, I1 and I2 have to be integers that sum up to population_[0]
#' @param population: vector of length (nObs) representing the total population at each month. In this model I'm assuming
#'              that the population stays constant between months (which is roughly true).
#' @param rainfall: vector of length (nObs * nSteps) representing the rainfall at each step.
#' @param nSteps: finite difference step in which one month is divided (mininum nStep = 1, i.e. you directly step from
#'          one month to the next).
#' @param nSimul: number of simulated paths
#' @param obsInterval: the length of the interval between 2 observation (which will be subdivided in nSteps). Typically
#'                     it is going to be 1/12 so each observation represents a month.
#' @param nObs: number of observations (months) of each simulated trajectory.
#' @param T0: initial time of the simulations. It is used to set up the spline basis.
#' @param splineSettings: named list where [["nBasis"]], [["degree"]] and [["period"]] are the (self-explanatory)
#'                        settings used to set up the spline basis for the seasonal component of the model.
#'
#' @return A named list where S1, S2, E, I1, I2, lambda1 and lambda2 are matrices of size (nSimul by nObs) each row
#'         representing a simulated trajectory of a particular state.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>                         
#' @export malSirSimul

malSirSimul <- function(procparam, obsparam, splineparam, initStates, population,
                        rainfall, nSimul, nObs, 
                        nSteps = 30,  obsInterval = 1/12,  T0 = 0.0, 
                        splineSettings = list("nBasis" = as.integer(6),  "degree" = as.integer(2), "period" = 1.0))
{  
  stopifnot(
    nSteps >= 1,
    nObs   >= 1,
    
    length(population) == nObs,
    length(rainfall) == (nObs * nSteps + 1),
    
    nrow(procparam) == nSimul,
    nrow(obsparam) == nSimul,
    nrow(splineparam) == nSimul,
    nrow(initStates) == nSimul,
    
    ncol(initStates) == 7,
    ncol(splineparam) == splineSettings["nBasis"]
  )
  
  # Rounding to integer S1, S2, E, I1 and I2 is important, otherwise the euler multinomial ampling gives NaNs
  initStates[ , 1:5] <- round( exp(initStates[ , 1:5, drop = FALSE]) )
  
  .Call("malSirSimulCpp", procparam_ = procparam, obsparam_ = obsparam, splineparam_ = splineparam, initStates_ = initStates, population_ = population,
                          rainfall_ = rainfall, nSteps_ = nSteps, nSimul_ = nSimul, obsInterval_ = obsInterval, nObs_ = nObs, T0_ = T0, splineSettings_ = splineSettings,
        PACKAGE = "synlik")
  
}


########
# Wrapper for synlik
########

malSirCreator <- function(population, rainfall, nSteps, obsInterval, T0, splineSettings)
{
  
  tmpFun <- function(param, nsim, extraArgs, ...)
  {
    if( !is.loaded("synlik") ) library("synlik")
    
    if( !is.matrix(param) ) param <- matrix(param, nsim, length(param), byrow = TRUE)
    
    if( is.null(extraArgs$nObs) ) stop("is.null(extraArgs$nObs)")
    nObs <- extraArgs$nObs

    procparam <- param[ , 1:10, drop = FALSE]
    obsparam <- param[ , 11:12, drop = FALSE]
    splineparam <- param[ , 13:18, drop = FALSE]
    
    # Dealing with initial values of states S1, S2, E, I1, I2 and
    # 1. Apply logistic transformation (because each number is a fraction)
    # 2. Normalized them so the sum to 1
    # 3. Multiply by the total population
    # 4. Round them to the nearest integer
    tmp <- logistic(param[ , 19:23, drop = FALSE])
    tmp <- round( tmp/rowSums(tmp) * population[1] )
    
    # Add lambda1 and lambda2 to initial states
    initStates <- cbind(tmp, exp(param[ , 24:25, drop = FALSE]))
    
    stopifnot(
      nSteps >= 1,
      nObs   >= 1,
      
      length(population) == nObs,
      length(rainfall) == nObs * nSteps,
      
      nrow(procparam) == nsim,
      nrow(obsparam) == nsim,
      nrow(splineparam) == nsim,
      nrow(initStates) == nsim,
      
      ncol(initStates) == 7,
      ncol(splineparam) == splineSettings["nBasis"]
    )
    
    tmp <- .Call("malSirSimulCpp",
                  procparam_ = procparam, 
                  obsparam_ = obsparam, 
                  splineparam_ = splineparam, 
                  initStates_ = initStates, 
                  population_ = population,
                  rainfall_ = rainfall, 
                  nSteps_ = nSteps,
                  nSimul_ = nsim,   
                  obsInterval_ = obsInterval, 
                  nObs_ = nObs,
                  T0_ = T0, 
                  splineSettings_ = splineSettings,
                  PACKAGE = "synlik")
    
    return( tmp[["Obs_Cases"]] )
    
  }
  
  return(tmpFun)
}



###############
#### Simulator for stable distribution
###############
#' Simulates from the exponential distribution
#' 
#' @param nObs Length of each simulated time serie.
#' @param nsim Number of simulations from the model.
#' @param param A vector of length 1 or a matrix of size nsim by 1.
#'        If param is a vector each of the nsim simulations will use the same parameters, if it's a matrix the i-th simulation
#'        will use the i-th row of param as parameters.
#' @param extraArgs Only needed for compatibility with the "synlik" package.
#' @return A matrix of size nsim by nObs where each row is a sample from the distribution.                                  
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>   
#' @export
#' 
expSimul <- function(param, nsim, extraArgs, nObs, ...)
{
  
  output <- matrix(NA, nsim, nObs)
  
  if(is.vector(param))
  {
    output <- matrix(rexp(nObs*nsim, exp(param)), nsim, nObs)
    return(output)
  } 
  
  if(is.matrix(param))
  {
    param <- exp(param)
    stopifnot( nrow(param) == nsim, ncol(param) == 1  )
    for(ii in 1:nsim) output[ii, ] <- rexp(nObs, exp(param[ii, ]))
    return(output)
  }
  
  stop("param has to be either a matrix of a vector") 
}













