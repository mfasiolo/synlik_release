
# Delta = death rate FIXED
# DeltaI = death rate of infected (m in paper)
# Alpha = fixed to 1

## Two stages only
# rho = duration of shorth term immunity
# clin = share that go from S to I (put it to 1 for other)
# Rs.0 = Y

## Seasonal
# log.omega2 to 6 

library(pompExamples)

pompExample(dacca)

cholPomp <- dacca

par(mfrow = c(1, 1))
plot(drop(dacca@data), type = 'l')
lines(drop(simulate(cholPomp)@data), col = 2)

tmp <- pfilter(cholPomp, Np = 1000)
tmp@loglik


#########################
######## MCMC estimation
#########################

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
    data(chol_prop_pmcmc)
    initPar <- dacca@params
    #proposal <- 0.1 * chol_prop_pmcmc  # 0.02 * initPar
    
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


tim <- proc.time()
tmp <- malariaPMCMC(cholPomp, 
                    initPar = initPar,
                    nIter = 10, 
                    burnIn = 0,
                    priorFun = cholPrior, 
                    propCov = 0.1 * proposal,
                    nsim = 1e3)
tim <- proc.time() - tim


# Plotting results
par(mfrow = c(6, 5))
for(ii in 1:ncol(tmp$mcmcSample)) 
{
  plot(tmp$mcmcSample[ , ii], type = 'l', ylab = names(cholPomp@params)[ii])
}


###############
# Multiple Runs
###############
library(parallel)

tim <- proc.time()
nCores <- 4
results <- mclapply(1:nCores,
                function(input){
                  
                  par <- cholPomp@params
                  cholPomp <- simulate(cholPomp, params = par)
                  
                  malariaPMCMC(cholPomp, 
                    initPar = initPar,
                    nIter = 10000, 
                    burnIn = 0,
                    priorFun = cholPrior, 
                    propCov = 0.1 * proposal,
                    nsim = 1e3)
                },
                mc.cores = nCores
)
tim <- proc.time() - tim

# chol_pmcmc_results <- results
# save(file = "~/Desktop/chol_pmcmc_results.RData", chol_pmcmc_results)

# Plotting results
tmpCov <- matrix(0, 30, 30)
for(pathInd in 1:nCores)
{
  tmpMat <- results[[pathInd]]$mcmcSample
  par(mfrow = c(6, 5))
  for(ii in 1:ncol(tmpMat))
  { 
    tmp <- tmpMat[ , ii]
    #if(names(initPar[ii]) %in% cholPomp@userdata$logvar) tmp <- exp(tmp)
    #if(names(initPar[ii]) %in% cholPomp@userdata$logitvar) tmp <- logistic(tmp)
    plot(tmp, type = 'l', main = names(cholPomp@params)[ii]) 
    abline(h = cholPomp@params[ii], col = 2)
  }
  
  tmpCov <- tmpCov + cov(tmpMat)
  #Sys.sleep(1)
}
# proposal <- tmpCov / nCores
# save(file = "~/Desktop/chol_prop_pmcmc.RData", proposal)

par(mfrow = c(6, 5))
 for(ii in 1:30)
 {
   tmp <- sapply(results, function(input) input$mcmcSample[-(1:1e4) , ii])
   hist(tmp, main = names(cholPomp@params)[ii])
   abline(v = cholPomp@params[ii], col = 2, lwd = 2)
 }


#################################
########## Trying MIF
#################################

tmp <- mif(object = cholPomp, 
           Nmif = 10, 
           start = initPar, 
           rw.sd = pmax(proposal, 1e-8),
           var.factor = 1,
           ic.lag = 40,
           cooling.fraction = 0.95,
           Np = 1000,
           trans = TRUE)
           
plot(tmp)

