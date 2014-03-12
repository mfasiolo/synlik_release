##########################
#### Budmoth with Pomp
##########################

library(pompExamples)
library(synlik)
library(parallel)

source("~/Desktop/All/Dropbox/Work/Programming/New_SL/synlik/Examples/bud_pomp_functions.R")
source("~/Desktop/All/Dropbox/Work/Programming/New_SL/synlik/Examples/bud_synlik_functions.R")

#######
# Create model
#######

pompModel <- budModCreator("tri")
truepar <- pompModel@params 

bud_sl <- new("synlik",
                 simulator = budSimul,
                 summaries = budStats,
                 param = pompModel@params,
                 extraArgs = list("pompModel" = pompModel)
)

#### Simulate from the object
bud_sl@data <- simulate(bud_sl)
bud_sl@extraArgs$obsData <- bud_sl@data

par(mfrow = c(2, 2))
for(ii in 1:3)
  plot(bud_sl@data[[ii]], type = 'l', ylab = names(bud_sl@data)[ii], main = names(bud_sl@data)[ii])

a <- simulate(bud_sl, nsim = 10, stats = F)
plotMatrix(cor(a), Lst = list())
abs(cor(a)) > 0.9

checkNorm(bud_sl)

synlikEval(bud_sl, 
           param  = bud_sl@param,
           nsim    = 1e3)

##################
########## Slicing the synthetic likelihood
##################

# Slicing the Synthetic Likelihood
nPoints <- 10
ranges <- list( list("alpha",  seq(0.3, 0.7, length.out = nPoints) ),       # 1
                list("sig.alpha", seq(0.01, 0.3, length.out = nPoints) ),   # 2 
                list("gam", seq(20, 100, length.out = nPoints) ),           # 3  
                list("lambda",  seq(20, 30, length.out = nPoints) ),      # 4
                list("sig.lambda", seq(0.05, 0.5, length.out = nPoints) ),       # 5
                list("g",   seq(0.04, 0.16, length.out = nPoints) ),       # 6
                list("delta",  seq(5, 15, length.out = nPoints) ), # 7 
                list("a",    seq(1, 2.5, length.out = nPoints) ), # 7
                list("sig.a",  seq(0.01, 0.5, length.out = nPoints) ),       # 8
                list("w",  seq(0.05, 0.3, length.out = nPoints) ),        # 9 
                list("beta0",  seq(-1, 1, length.out = nPoints) ),       # 10
                list("beta1",  seq(25, 45, length.out = nPoints) ),       # 11
                list("u",  seq(0.7, 0.99, length.out = nPoints) ), # 12
                list("sigQobs", seq(0.01, 0.1, length.out = nPoints) ),               # 13
                list("sigNobs", seq(0.1, 1, length.out = nPoints) ),               # 14
                list("sigSobs", seq(0.01, 0.3, length.out = nPoints) ),                 # 15
                list("Q.0", seq(0.5, 0.99, length.out = nPoints) ),               # 16
                list("N.0", seq(0.001, 0.05, length.out = nPoints) ),                # 17
                list("S.0", seq(0.05, 0.5, length.out = nPoints) )                # 18
)


indexes <- c(1:19)
tmp <- mclapply(ranges[indexes], 
                function(val){
                  try(synSlice(object = bud_sl, parName = val[[1]], range = val[[2]], 
                               nsim = 1e3))
                },
                mc.cores = 4)

par(mfrow = c(5, 4))
kk <- 1
for(ii in indexes)
{
  plot(ranges[[ii]][[2]], tmp[[kk]], type = 'l', 
       main = ranges[[ii]][[1]], xlab = ranges[[ii]][[1]], ylab = "log-likelihood")
  abline(v = bud_sl@param[[ii]], lty = 2)
  kk <- kk + 1
}


##########
#### MCMC estimation
##########
# Put the parameters in a transformed scale
initPar <- pompModel@params 
initPar[pompModel@userdata$logitvar] <- ilogistic(initPar[pompModel@userdata$logitvar])
initPar[pompModel@userdata$logvar] <- log(initPar[pompModel@userdata$logvar])

tim <- proc.time()
bud_sl    <- synMcmc(bud_sl,
                     initPar = initPar, 
                     nIter = 10, 
                     nsim = 500,
                     propCov = 0.1 * diag(rep(0.001, 19)^2), 
                     burnIn = 0,
                     priorFun = function(input, ...){
                       names(input) <- names(pompModel@params)
                       sum(input[pompModel@userdata$logvar]) + # inverting log-transformation
                       sum(-input[pompModel@userdata$logitvar] - # inverting inverse logistic transformation 
                             2*log(1 + exp(-input[pompModel@userdata$logitvar])))
                     },
                     recompute = FALSE,
                     trans = TRUE,
                     multicore = FALSE,
                     control= list("ncores" = 2))
tim <- proc.time() - tim

plot(bud_sl,
     trans = list(
     parName = c(pompModel@userdata$logvar, pompModel@userdata$logitvar), 
     transform = c(rep("exp", length(pompModel@userdata$logvar)), rep("logistic", length(pompModel@userdata$logitvar)))
     ))



########## Parallel runs

# Put the parameters in a transformed scale
initPar <- pompModel@params 
initPar[pompModel@userdata$logitvar] <- ilogistic(initPar[pompModel@userdata$logitvar])
initPar[pompModel@userdata$logvar] <- log(initPar[pompModel@userdata$logvar])

data(prop_bud_synlik)

nCores <- 3
tim <- proc.time()
results    <-   mclapply(1:nCores, 
                     function(input){
                       
                       bud_sl@data <- simulate(bud_sl)
                       bud_sl@extraArgs$obsData <- bud_sl@data
                       
                       synMcmc(bud_sl,
                               initPar = initPar, 
                               nIter = 45000, 
                               nsim = 1000,
                               propCov = 0.01 * prop_bud_synlik, 
                               burnIn = 0,
                               priorFun = function(input, ...){
                                 names(input) <- names(pompModel@params)
                                   sum(input[pompModel@userdata$logvar]) + # inverting log-transformation
                                   sum(-input[pompModel@userdata$logitvar] - # inverting inverse logistic transformation 
                                         2*log(1 + exp(-input[pompModel@userdata$logitvar]))) +
                                   dunif(exp(input[c("sig.alpha", "sig.lambda", "sig.a", "sigQobs", "sigNobs", "sigSobs")]), 0.01, 100, log = TRUE)
                               },
                               recompute = FALSE,
                               trans = TRUE)
                       },
                     mc.cores = nCores)
tim <- proc.time() - tim

tmpCov <- matrix(0, 19, 19)
for(pathInd in 1:nCores)
{
  tmpMat <- results[[pathInd]]@mcmcChain
  par(mfrow = c(5, 4))
  for(ii in 1:ncol(tmpMat))
  { 
    tmp <- tmpMat[ , ii]
    if(names(initPar[ii]) %in% pompModel@userdata$logvar) tmp <- exp(tmp)
    if(names(initPar[ii]) %in% pompModel@userdata$logitvar) tmp <- logistic(tmp)
    plot(tmp, type = 'l', main = names(results[[pathInd]]@param)[ii]) 
    abline(h = pompModel@params[ii], col = 2)
  }
   
  tmpCov <- tmpCov + cov(tmpMat)
  #Sys.sleep(1)
}
### Checking convergence of the chains
# prop_bud_synlik <- tmpCov / nCores
# save(file = "~/Desktop/prop_bud_synlik.RData", prop_bud_synlik)

upper <- sapply(ranges, function(input) tail(input[[2]], 1))
lower <- sapply(ranges, function(input) head(input[[2]], 1))
names(lower) <- names(upper) <- names(bud_sl@param)
newStats <- shrinkStat(bud_sl, 
                       nsim = 100, 
                       upper = upper, 
                       lower = lower,
                       method = "lasso")




#############################
# New stoc optim
#############################

data(prop_bud_synlik)

#### Simulate from the object
bud_sl@data <- simulate(bud_sl)
bud_sl@extraArgs$obsData <- bud_sl@data

# Put the parameters in a transformed scale

constr <- list("indexes" = c(which(names(pompModel@params) %in% pompModel@userdata$logitvar),
                             which(names(pompModel@params) %in% pompModel@userdata$logvar)),
               "upper" = c(rep(1, 4), rep(100, 14)) , "lower" = rep(0, 18))

initPar <- pompModel@params 
initPar <- c(0.2, 0.2, 35, 15, 0.4, 0.3, 6, 3, 0.3, 0.2, 0, 25, 0.8, 0.1, 0.4, 0.2, 0.6, 0.1, 0.4)

nIter = 300
nP = 50

## New procedure
res <- sml(object = bud_sl, 
           init = initPar, 
           initCov = 0.2 * diag(pompModel@params),
           nP = nP,
           nsim = 500, 
           nIter = nIter, 
           alpha = 0.99,
           recycle = F,
           multicore = T,
           ncores = 5,
           saddle = F,
           constr = constr#,
           #temper = c(seq(0.1, 1, length.out = nIter/5), rep(1, nIter * 0.8))
)

par(mfrow = c(5, 4))
for(ii in 1:length(pompModel@params)) {
  tmp <- res$simPar[ , ii]
  poi <- res$estim[ , ii]

  plot(rep(1:nIter, each = nP), tmp, ylab = names(bud_sl@param)[ii], main = names(bud_sl@param)[ii],
       ylim = c(min(c(poi, tmp)), max(c(poi, tmp))))
  lines(poi, type = 'l', col = 2, lwd = 3)
  abline(h = pompModel@params[ii], lty = 2, lwd = 2, col = 3)
  points(1, initPar[ii], col = 2, pch = 2, lwd = 3)
  #points(rep(1:nIter, each = nP), tmp)
}

sum(res$simLik == 0)