#######################################################################
##############    PMCMC for malaria model
####################################################################### 

library(synlik)
library(pomp)
library(parallel)

source("mal_pomp_fun.R")

############
# Set up the pomp model
############

T0 <- - 3 / 12
nObs <- 240
obsInterval <- 1 / 12
nSteps <- 30

data(malaria_data)

tbasis <- seq(T0, T0 + nObs * obsInterval, by = obsInterval / nSteps)
basis <- periodic.bspline.basis(tbasis, nbasis=6, degree=2, period=1)


# Create the pomp model, rprocess uses a closure, so that a lot of data (population, rain, etc) are fixed.
malariaPomp <- pomp(malariaPomp, 
                    time = 1:240,
                    rprocess = malSirPomCreator(population = malaria_data$population, 
                                                rainfall = setupRain(malaria_data$rain, malaria_data$year, nSteps, obsInterval, lag = 5, v = 200),
                                                nSteps   = nSteps, 
                                                obsInterval = obsInterval, 
                                                splineBasis = periodic.bspline.basis(tbasis, nbasis=6, degree=2, period=1)), 
                    params = c(log(c(mu_EI1 = 8.902, 
                                     mu_I1I2 = 5.511, 
                                     mu_I2S2 = 0.035, 
                                     mu_S2S1 = 0.334, 
                                     mu_I1S1 = 6.563,
                                     q = 9.424 * 1e-4,
                                     covarCoeff = 0.512,
                                     sigma = 0.243,
                                     c = 0.01,
                                     tau = 0.03)), 
                               log(c( rho = 0.015, psiSquare = 0.395^2)), 
                               c("spl1" = 1.201, "spl2" = 2.088, "spl3" = 3.866, "spl4" = 2.808, "spl5" = 5.996, "spl6" = 5.333),
                               ilogistic(c("S1.0" = 0.138, "S2.0" = 0.775, "E.0"  = 0.004, "I1.0" = 0.002, "I2.0" = 0.08)),
                               log(c("LAM1.0" = 0.0171, "LAM2.0" = 0.0061)), "H_Cases.0" = 0)
                    )


noRainPar <- c(log(c(mu_EI1 = 7.408, 
                     mu_I1I2 = 11.544, 
                     mu_I2S2 = 0.004, 
                     mu_S2S1 = 0.23, 
                     mu_I1S1 = 2.32,
                     q = 4.763 * 1e-4,
                     covarCoeff = 1e-6,
                     sigma = 0.309,
                     c = 0.04,
                     tau = 0.022)), 
               log(c( rho = 0.03, psiSquare = 0.39^2)), 
               c("spl1" = -2.469, "spl2" = 2.001, "spl3" = 4.227, "spl4" = 2.786, "spl5" = 6.534, "spl6" = 7.08),
               ilogistic(c("S1.0" = 0.164, "S2.0" = 0.765, "E.0"  = 0.002, "I1.0" = 0.02, "I2.0" = 0.067)),
               log(c("LAM1.0" = 0.0133, "LAM2.0" = 0.0045)), "H_Cases.0" = 0)

#########
# Try the simulator
#########
malariaPomp@params <- noRainPar
tmp <- simulate(malariaPomp, nsim = 100)

tmp <- t(sapply(tmp, function(input) input@data))

par(mfrow = c(1, 1))
{
  plot((1:nObs) * obsInterval,  colMeans(tmp), type = 'l', xlab = "Year", ylab = "Obs_cases", 
       main = "Obs_cases", xaxt = 'n', ylim = c(0, 0.6e4) )
  polygon(x = c((1:nObs) * obsInterval, rev((1:nObs) * obsInterval)), 
          y = c(apply(tmp, 2, quantile, 0.1), rev(apply(tmp, 2, quantile, 0.9))), col = 'grey',
          border = FALSE)
  lines((1:nObs) * obsInterval,  colMeans(tmp))
  axis(side = 1, at=0:20, labels=0:20 + 87)
  lines((1:nObs) * obsInterval, malaria_data$cases, col = 2)
}


##########
# Put some data into the object
##########

malariaPomp <- simulate(malariaPomp)
plot(malariaPomp)


###########
# Trying the particle filter
###########

malariaPomp <- pfilter(malariaPomp, max.fail = nObs, Np=1000)
malariaPomp$loglik

####################################
########### MCMC
####################################
malariaPomp@data <- matrix(malaria_data$cases, 1, 240)
rownames(malariaPomp@data) <- c("Y")

initPar <- malariaPomp@params
initPar["psiSquare"] = 8 * 0.395^2

#### MCMC estimation
tim <- proc.time()
tmp <- malariaPMCMC(malariaPomp, 
                            initPar = initPar,
                            nIter = 10, 
                            burnIn = 0,
                            priorFun = function(input){
                              sum(input[c(1:12, 24:25)]) + # inverting log-transformation
                                sum(-input[19:23] - # inverting inverse logistic transformation 
                                      2*log(1 + exp(-input[19:23]))) +
                                dunif(exp(input[1]), 1, 30, log = T)   +
                                dunif(exp(input[2]), 0.5, 10, log = T) +
                                dunif(exp(input[3]), 0, 0.5, log = T)  +
                                dunif(exp(input[4]), 0, 2, log = T)    +
                                dunif(exp(input[5]), 0, 20, log = T)   +
                                dunif(exp(input[6]), 0, 0.2, log = T)  +
                                dunif(exp(input[7]), 0.1, 1, log = T)  +
                                dunif(exp(input[8]), 0, 2, log = T)    +
                                dunif(exp(input[9]), 0, 0.2, log = T)  +
                                dunif(exp(input[11]), 0.001, 0.1, log = T)  +
                                dunif(exp(input[24]), 0, 0.1, log = T) +
                                dunif(exp(input[25]), 0, 0.1, log = T) 
                            }, 
                            recompute = FALSE,
                            propCov = diag(c(rep(0.01, 12), 
                                             rep(0.01, 6 ), 
                                             rep(0.01, 7), 1e-3)) ^ 2,
                            nsim = 1e3)
tim <- proc.time() - tim

# Plotting results
par(mfrow = c(5, 5))
for(ii in 1:(ncol(tmp$mcmcSample)-1)) 
{
  plot(tmp$mcmcSample[ , ii], type = 'l', ylab = names(malariaPomp@params)[ii])
}



############################
###### Debugging the simulator
############################

# 
# T0 <- - 3 / 12 # Different from synlik in order to be in phase
# nObs <- 240
# obsInterval <- 1 / 12
# nSteps <- 30
# 
# data(malaria_data)
# 
# tbasis <- seq(T0, T0 + nObs * obsInterval, by = obsInterval / nSteps)
# basis <- periodic.bspline.basis(tbasis, nbasis=6, degree=2, period=1)
# 
# stepFun <- malSirPomCreator(population = malaria_data$population, 
#                             rainfall = setupRain(malaria_data$rain, malaria_data$year, nSteps, obsInterval, lag = 5, v = 200),
#                             nSteps   = nSteps, 
#                             obsInterval = obsInterval, 
#                             splineBasis = periodic.bspline.basis(tbasis, nbasis=6, degree=2, period=1))
# 
# nsimul <- 2
# mah <- stepFun(xstart = t(matrix(rep(c(ilogistic(c(S1 = 0.138, S2 = 0.775, E  = 0.004, I1 = 0.002, I2 = 0.08)),
#                                        log(c(LAM1 = 0.0133, LAM2 = 0.0045)), 0), nsimul), nsimul, 8, byrow = TRUE)), 
#                times = 0:240, 
#                params = matrix(rep(c(log(c(mu_EI1 = 8.902, 
#                                            mu_I1I2 = 5.511, 
#                                            mu_I2S2 = 0.035, 
#                                            mu_S2S1 = 0.334, 
#                                            mu_I1S1 = 6.563,
#                                            q = 9.424 * 1e-4,
#                                            covarCoeff = 0.512,
#                                            sigma = 0.243,
#                                            c = 0.01,
#                                            tau = 0.03)), 
#                                      log(c( rho = 0.015, psiSquare = 0.395^2)), 
#                                      c("spl1" = 1.201, "spl2" = 2.088, "spl3" = 3.866, "spl4" = 2.808, "spl5" = 5.996, "spl6" = 5.333),
#                                      ilogistic(c(S1 = 0.138, S2 = 0.775, E  = 0.004, I1 = 0.002, I2 = 0.08)),
#                                      log(c(LAM1 = 0.0133, LAM2 = 0.0045)), 0), nsimul),
#                                26, nsimul))
# 
# par(mfrow = c(1, 1))
# plot(apply(mah, 3, function(input) input[8, 2]), type = 'l')




