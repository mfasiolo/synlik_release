###################################
####### synlik package example use
###################################

library(synlik)

####################################################################################
##################### Example use with stable map
####################################################################################

# Stable distributions (see ?rstable) have intractable likelihoods, hence we use synthetic likelihood.
# We are interested in estimating 4 parameters: alpha = 1.5, beta = 0.1, gamma = 1, delta = 2

#### Create "synlik" object
# We need: the number of observations, 
#          a function that simulates data ( a wrapper for rstable() ), 
#          a function that transforms the data into parameters ( 15 empirical quantiles in this case )
#          a named vector of parameters

stable_sl <- new("synlik", 
                 simulator = stableSimul,
                 summaries = stableStats,
                 param = c(alpha = log(1.5), beta = 0.1, gamma = log(1), delta = 2),
                 extraArgs = list("nObs" = 1000)
                 
)

#### Simulate from the object
stable_sl@data <- simulate(stable_sl)

#### Plotting the object
plot(stable_sl)

#### Checking multivariate normality
checkNorm(stable_sl)

#### Evaluate the likelihood
synlikEval(stable_sl, 
          param  = c(alpha = log(1.5), beta = 0.1, gamma = log(1), delta = 2),
          nsim    = 1e3)

#### Plotting a slice of the logLikelihood on multiple cores
synSlice(object = stable_sl, parName = "alpha", range = log(seq(1.2, 1.7, by = 0.05)), 
        param = c(alpha = log(1.5), beta = 0.1, gamma = log(1), delta = 2), 
        nsim = 1e3, saddle = TRUE, multicore = TRUE)

#### Plotting a slice of the logLikelihood wrt 2 parameters on multiple cores (might take 1 minute)
synSlice(object = stable_sl, parName = c("alpha", "gamma"),
        range = list("alpha" = log(seq(1.2, 1.7, by = 0.05)), "gamma" = seq(0.9, 1.1, by = 0.02)), 
        param = c(alpha = log(1.5), beta = 0.1, gamma = log(1), delta = 2), 
        nsim = 1e3, multicore = TRUE)

####################################
########### MCMC
####################################

#### MCMC estimation on multiple cores
stable_sl <- synMcmc(stable_sl, 
                    initPar = c(alpha = log(1.7), beta = -0.1, gamma = log(1.3), delta = 1.5),
                    nIter = 50, 
                    burnIn = 10,
                    priorFun = function(input) { input[1] > log(1) && input[1] < log(2) && abs(input[2]) < 1 }, 
                    propCov = diag(c(0.1, 0.1, 0.1, 0.1))^2, 
                    nsim = 1e3, 
                    multicore = T)

# Continue from were we left the mcmc estimation, with the same settings.
stable_sl <- continue(stable_sl)

#### Plotting MCMC chains
plot(stable_sl, trans = list("parIndex" = c(1, 3), "transform" = rep("exp", 2)))


###########################################
############## Synthetic Maximum likelihood
###########################################

stable_sl@data <- simulate(stable_sl)

# We try to maximize the synthetic likelihood (simulation on multiple cores) (might take 1 minute or 2)
#set.seed(432)
stable_sl <- synMaxlik(object  = stable_sl, 
                        nIter = 40,
                        nsim = 1000, 
                        initCov = diag(c(0.1, 0.1, 0.1, 0.1))^2, 
                        initPar = c(alpha = log(1.7), beta = -0.1, gamma = log(1.3), delta = 1.5), 
                        addRegr = TRUE,
                        constr = list("indexes" = c(1, 2), "upper" = c(log(2), 1), "lower" = c(log(1), -1) ),
                        multicore = T)

# We do some more iterations
stable_sl <- continue.synMaxlik(stable_sl)

plot(stable_sl, trans = list("parIndex" = c(1, 3), "transform" = rep("exp", 2)))

coef(stable_sl)

## Checking coverage

a <- checkCoverage.synMaxlik(object  = stable_sl, 
                             nRep = 100,
                             nIter = 50,
                             nsim = 1000, 
                             initCov = diag(c(0.1, 0.1, 0.1, 0.1))^2, 
                             initPar = c(alpha = log(1.7), beta = -0.1, gamma = log(1.3), delta = 1.5), 
                             addRegr = TRUE,
                             constr = list("indexes" = c(1, 2), "upper" = c(log(2), 1), "lower" = c(log(1), -1) ),
                             multicore = FALSE)

