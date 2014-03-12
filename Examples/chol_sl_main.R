##########################
#### Cholera with SL
##########################

library(pompExamples)
library(synlik)
library(parallel)

source("~/Desktop/All/Dropbox/Work/Programming/New_SL/synlik/Examples/chol_sl_fun.R")

#######
# Create model
#######

pompExample(dacca)
cholPomp <- dacca

chol_sl <- new("synlik",
              simulator = cholSimul,
              summaries = cholStats,
              param = cholPomp@params,
              extraArgs = list("pompModel" = cholPomp)
)

#### Simulate from the object
chol_sl@data <- drop(cholPomp@data)
chol_sl@extraArgs$obsData <- chol_sl@data

plot(chol_sl)

# Simulating statistics
a <- simulate(chol_sl, nsim = 1000, stats = T, clean = TRUE)
plotMatrix(cor(a), Lst = list())

# Estimating likelihood
synlikEval(chol_sl, 
           param  = chol_sl@param,
           nsim    = 1000)

constr <- list("indexes" = which(names(cholPomp@params) %in% 
                                   c("gamma", "eps", "deltaI", "sd.beta", "tau", "S.0", "I.0", "R1.0", "R2.0", "R3.0")), 
               "upper" = rep(1e3, 10), "lower" = rep(1e-8, 10))

# Statistics selection
chol_sl@simulator <- cholSimulMulti
a <- shrinkCoef(chol_sl, nsim = 1000, mu = chol_sl@param, sigma = 0.2*proposal, constr = constr, type = "lasso", clean = TRUE)
apply(a$regrCoef, 2, function(input) sum(input, na.rm = TRUE))
chol_sl@simulator <- cholSimul
 
#### Plotting a slice of the logLikelihood
N <- 14
synSlice(object = chol_sl, 
         ranges = list("gamma" =       seq(    10,    30,  length.out = N), 
                       "eps"   =       seq(    10,   140,  length.out = N), 
                       "deltaI" =      seq(  0.04,   0.2,   length.out = N), 
                       "beta.trend" =  seq( -0.02,  0.01, length.out = N),
#                        "log.beta1" =   seq(   0.2,   1.5, length.out = N),
#                        "log.beta2" =   seq(     5,     7, length.out = N),
#                        "log.beta3" =   seq(    -5,  -2.5, length.out = N),
#                        "log.beta4" =   seq(     2,     6, length.out = N),
#                        "log.beta5" =   seq(     1,     6, length.out = N),
#                        "log.beta6" =   seq(     2,     7, length.out = N),
#                        "log.omega1" =  seq(    -5,     0, length.out = N),
#                        "log.omega2" =  seq(    -7,  -0.5, length.out = N),
#                        "log.omega3" =  seq(    -6,    -1, length.out = N),
#                        "log.omega4" =  seq(    -6,    -3, length.out = N),
#                        "log.omega5" =  seq(   -20,    -4, length.out = N),
#                        "log.omega6" =  seq(   -12,     0, length.out = N),
#                        "sd.beta"    =  seq(     2,   4.5, length.out = N),
                       "tau"      =    seq(   0.1,   0.4, length.out = N),
#                        "S.0"        =  seq(   0.1,  0.99, length.out = N),
#                        "I.0"        =  seq(  0.01,   0.5, length.out = N),
#                        "R1.0"        = seq(0.0002, 0.001, length.out = N),
#                        "R2.0"        = seq(0.0002, 0.002, length.out = N),
                       "R3.0"        = seq(0.0002, 0.002, length.out = N)),  
         nsim = 1000, multicore = T, saddle = F)

####################
###### Maximum Synthetic likelihood
####################

model <- "reserv"
if(model == "twoPath")
{
  initPar <- c("gamma" = 8, "eps" = 0.7, "rho" = 7.1, "delta" = 0.02, "deltaI" = 9.238,
               "clin" = 0.0051, "alpha" = 1, "beta.trend" = -6.8, "log.beta1" = 6.5,
               "log.beta2" = 11.4, "log.beta3" = 2.1, "log.beta4" = 9, "log.beta5" = 8.6,
               "log.beta6" = 9.5, "log.omega1" = -4.5, "log.omega2" = -10, "log.omega3" = -10,
               "log.omega4" = -10, "log.omega5" = -10, "log.omega6" = -10, "sd.beta" = 639.6,
               "tau" = 0.23, "S.0" = 0.62, "I.0" = 0.378, "Rs.0" = 0.00001, "R1.0" = 0.000843,
               "R2.0" = 0.000972, "R3.0" = 0.000000116, "nbasis" = 6L, "nrstage" = 3)
  cholPomp@covar[ , 1] <- 1
} else{
  if(model == "reserv")
  {
    data(chol_sl_results)
    tmp <- Reduce("rbind", lapply(chol_sl_results, function(input) input@mcmcChain[-(1:5000), ]))
    initPar <- dacca@params
    proposal <- 0.1 * cov(tmp) / length(chol_sl_results)
    
    #proposal[c("rho", "clin", "Rs.0",                           # Just two-stage model 
    #           "delta", "alpha", "nbasis", "nrstage")] <- 0     # Fixed parameters
    
    cholPrior <- function(input, ...){
      names(input) <- names(initPar)
      toCheck <- c("gamma", "eps", "deltaI", "sd.beta", "tau", "S.0", "I.0", "R1.0", "R2.0", "R3.0")
      sum(dunif(input[toCheck], 1e-8, 1e3, log = T)) +
        sum(input[c("log.beta1", "log.beta2", "log.beta3", "log.beta4", "log.beta5", "log.beta6",
                    "log.omega1", "log.omega2", "log.omega3", "log.omega4", "log.omega5", "log.omega6")])
    }
  }
}

# Put the parameters in a transformed scale
tim <- proc.time()
chol_sl    <- synMcmc(chol_sl,
                     initPar = chol_sl@param, 
                     nIter = 1000, 
                     nsim = 1000,
                     propCov = 0.1 * proposal, 
                     burnIn = 0,
                     priorFun = cholPrior,
                     recompute = FALSE,
                     multicore = TRUE,
                     control= list("ncores" = 6))
tim <- proc.time() - tim

plot(chol_sl)


############
##### Parallel runs
############

nCores <- 6
tim <- proc.time()
results    <-   mclapply(1:nCores, 
                         function(input){
                           
                           chol_sl@data <- simulate(chol_sl)
                           chol_sl@extraArgs$obsData <- chol_sl@data
                                                      
                           synMcmc(chol_sl,
                                   initPar = initPar, 
                                   nIter = 15000, 
                                   nsim = 1000,
                                   propCov = proposal, 
                                   burnIn = 0,
                                   priorFun = cholPrior,
                                   recompute = FALSE)
                         },
                         mc.cores = nCores)
tim <- proc.time() - tim

# chol_sl_results <- results
# save(file = "~/Desktop/chol_sl_results.RData", chol_sl_results)


#tmpCov <- matrix(0, 30, 30)
for(pathInd in 1:nCores)
{
  tmpMat <- results[[pathInd]]@mcmcChain
  par(mfrow = c(6, 5))
  for(ii in 1:ncol(tmpMat))
  { 
    tmp <- tmpMat[ , ii]
    plot(tmp, type = 'l', main = names(results[[pathInd]]@param)[ii]) 
    abline(h = chol_sl@param[ii], col = 2)
  }
  
  #tmpCov <- tmpCov + cov(tmpMat)
  Sys.sleep(1)
}

allChain <- Reduce("rbind", lapply(results, function(input) input@mcmcChain[-(1:5000), ]))
allLik <- Reduce("rbind", lapply(results, function(input) input@LogLikChain[-(1:5000) ]))

par(mfrow = c(6, 5))
for(ii in 1:ncol(tmpMat))
{ 
  plot(allChain[ , ii], type = 'l', main = names(results[[pathInd]]@param)[ii]) 
  abline(h = chol_sl@param[ii], col = 2)
}

par(mfrow = c(1, 1))
par <- colMeans(allChain)
names(par) <- names(chol_sl@param)
plot(drop(simulate(chol_sl, param = par)), type = 'l')
lines(drop(dacca@data), col = 2)

plot(drop(simulate(chol_sl, param = par))^0.2, drop(dacca@data)^0.2)

N <- 14
synSlice(object = chol_sl, 
         ranges = list("eps"  =   seq(    10,   140,  length.out = N)),
         param = par,
         nsim = 1000, multicore = T, saddle = F)


seqEps <- seq(10, 140, length.out = 10)
tmp <- seqEps
for(ii in 1:10){
  param <- par
  param[2] <- seqEps[ii]
  tmp[ii] <- pfilter(cholPomp, params = param, Np = 1000)@loglik
}

plot(seqEps, tmp, type = 'l')






####################
###### Maximum Synthetic likelihood
####################

data(chol_prop_pmcmc)

constr <- list("indexes" = which(names(cholPomp@params) %in% 
                                       c("gamma", "eps", "deltaI", "sd.beta", "tau", "S.0", "I.0", "R1.0", "R2.0", "R3.0")), 
                "upper" = rep(1e3, 10), "lower" = rep(1e-8, 10))

nP <- 30
nIter <- 2000
## New procedure
tim <- proc.time()
res <- sml(object = chol_sl, 
           init = cholPomp@params, 
           initCov = diag(diag(proposal)),
           nP = nP,
           nsim = 500, 
           nIter = nIter, 
           alpha = 0.999,
           recycle = FALSE,
           multicore = TRUE,
           ncores = 7,
           constr = constr,
           saddle = FALSE
)
tim <- proc.time() - tim

# Plotting results
par(mfrow = c(6, 5))
for(ii in 1:ncol(res$estim)) 
{
  plot(rep(1:nIter, each = nP), res$simPar[ , ii], ylab = names(cholPomp@params)[ii], main = names(cholPomp@params)[ii],
       ylim = c(min(c(res$estim[ , ii], res$simPar[ , ii])), max(c(res$estim[ , ii], res$simPar[ , ii]))))
  lines(res$estim[ , ii], type = 'l', col = 2, lwd = 3)
}

par(mfrow = c(1, 1))
plot(simulate(chol_sl, param = res$estim[nIter, ]), type = 'l')

