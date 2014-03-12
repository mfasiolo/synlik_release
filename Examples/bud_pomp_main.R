##########################
#### Budmoth with Pomp
##########################

library(synlik)
library(pompExamples)

source("~/Desktop/All/Dropbox/Work/Programming/New_SL/synlik/Examples/bud_pomp_functions.R")
source("~/Desktop/All/Dropbox/Work/Programming/New_SL/synlik/Examples/bud_synlik_functions.R")

########
#### Create model and simulate data
########

pompModel <- budModCreator("tri")
pompModel <- simulate(pompModel)
plot(pompModel)

truepar <- pompModel@params 

# Trying particle filter
pompModel <- pfilter(pompModel, Np = 1000)
pompModel@loglik

###########
## MCMC
###########

# Put the parameters in a transformed scale
initPar <- pompModel@params 
initPar[pompModel@userdata$logitvar] <- ilogistic(initPar[pompModel@userdata$logitvar])
initPar[pompModel@userdata$logvar] <- log(initPar[pompModel@userdata$logvar])

tmp <- budPMCMC(object = pompModel, 
                initPar = initPar, 
                nIter = 1000, 
                nsim = 500,
                propCov = diag(rep(0.0005, length(initPar))^2), 
                burnIn = 0,
                priorFun = function(input, ...){
                  names(input) <- names(pompModel@params)
                  sum(input[pompModel@userdata$logvar]) + # inverting log-transformation
                    sum(-input[pompModel@userdata$logitvar] - # inverting inverse logistic transformation 
                          2*log(1 + exp(-input[pompModel@userdata$logitvar])))
                },
                recompute = FALSE,
                trans = TRUE)



par(mfrow = c(5, 4))
for(ii in 1:ncol(tmp$mcmcSample)){
  parNam <- names(pompModel@params)[ii]
  res <- tmp$mcmcSample[ , ii]
  
  if(parNam %in% pompModel@userdata$logitvar ) res <- logistic(res)
  if(parNam %in% pompModel@userdata$logvar ) res <- exp(res)
  
  plot(res, type = 'l', main = parNam)
  abline(h = truepar[ii], col = 2)
}

########## Parallel runs

# Put the parameters in a transformed scale
initPar <- pompModel@params 
initPar[pompModel@userdata$logitvar] <- ilogistic(initPar[pompModel@userdata$logitvar])
initPar[pompModel@userdata$logvar] <- log(initPar[pompModel@userdata$logvar])

data(prop_bud_pmcmc)

nCores <- 4
tim <- proc.time()
results    <-   mclapply(1:nCores, 
                         function(input){
                           
                           pompModel <- simulate(pompModel)
                           
                           budPMCMC(pompModel,
                                   initPar = initPar, 
                                   nIter = 10000, 
                                   nsim = 500,
                                   propCov = 0.1 * prop_bud_pmcmc, #diag(rep(0.0005, length(initPar))^2), 
                                   burnIn = 0,
                                   priorFun = function(input, ...){
                                     names(input) <- names(pompModel@params)
                                     sum(input[pompModel@userdata$logvar]) + # inverting log-transformation
                                       sum(-input[pompModel@userdata$logitvar] - # inverting inverse logistic transformation 
                                             2*log(1 + exp(-input[pompModel@userdata$logitvar]))) +
                                       dunif(exp(input[c("sig.alpha", "sig.lambda", "sig.a", "sigQobs", "sigNobs", "sigSobs")]), 0.001, 100, log = TRUE)
                                   },
                                   recompute = FALSE,
                                   trans = TRUE)
                         },
                         mc.cores = nCores)
tim <- proc.time() - tim

tmpCov <- matrix(0, 19, 19)
for(pathInd in 1:nCores)
{
  tmpMat <- results[[pathInd]]$mcmcSample
  par(mfrow = c(5, 4))
  for(ii in 1:ncol(tmpMat))
  { 
    tmp <- tmpMat[ , ii]
    if(names(initPar[ii]) %in% pompModel@userdata$logvar) tmp <- exp(tmp)
    if(names(initPar[ii]) %in% pompModel@userdata$logitvar) tmp <- logistic(tmp)
    plot(tmp, type = 'l', main = names(pompModel@params)[ii]) 
    abline(h = pompModel@params[ii], col = 2)
  }
  
  tmpCov <- tmpCov + cov(tmpMat)
  #Sys.sleep(1)
}

# prop_bud_pmcmc <- tmpCov / nCores
# save(file = "~/Desktop/prop_bud_pmcmc.RData", prop_bud_pmcmc)




