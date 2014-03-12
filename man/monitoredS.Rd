\docType{data}
\name{monitoredS}
\alias{monitoredS}
\title{monitoredS: Results of a synthetic likelihood MCMC}
\usage{
  MCMC(object, nbSimul = 600, FunctionSample = OmniSample,
    upFreq = 100, saveFreq = 20, sdProp = 0.4,
    adaptOK = FALSE, checkAdapt = 20, lowAcceptRate = 0.15,
    highAcceptRate = 0.4, useAutoStop = TRUE,
    checkAutoStop = 600,
    monitorFile = "thetasamples_all.txt", ...)

  MCMC
}
\arguments{
  \item{myData}{All data necessary to the Model as a list.
  It should contain: \itemize{ \item{\code{initValues}} a
  vector with the starting values for each parameter in
  order \item{\code{parmNames}} the names of the parameters
  sampled, in order \item{\code{monNames}} the names of the
  monitored variables, as returned by \code{Model} }
  Additional data can an should be provided to Model() as
  other named items in the list.}

  \item{Model}{The function returning a likelihood given a
  set of parameter and myData. #' The only requirement for
  the output of this function is to be compatible with
  \code{FunctionSample}, which use this output, normally a
  likelihood, to guide the mcmc sampling. If using the
  default \code{OmniSample} please follow the
  specifications in the help of this function.  Finally, if
  the output is a list and a slot \code{monitor} is defined
  it gives the values of variables to save and check the
  convergence on by default the parameter set is used.}

  \item{FunctionSample}{A function wich returns the next
  sampled parameter set as a vector and take as entries, in
  order: \itemize{ \item{Model} as defined above
  \item{myData} as defined above \item{oldTheta} the
  current set of parameters \item{nameParam} the name of
  the parameter to be sampled \item{sdProp} the standard
  deviation for this parameter \item{repeatSample} should
  the likelihood of the current theta be recomputed (0/1 or
  FALSE/TRUE) } This function is the heart of the MCMC
  procedure and is the one that can be specially adapted to
  handle any singularities of your model/likelihood.
  OmniSample is our one size fit all answer to the
  currently identified singularities of the synthetic
  likelihood.  Note that \code{myData} should contain and
  additional slot \code{myData$sampling} if the default
  \code{OmniSample} is used.}

  \item{nbSimul}{this is the number of iterations to be
  performed in non-autostop mode \code{useAutoStop =
  FALSE}. In auto-stop mode it should be left at the
  default value.}

  \item{upFreq}{screen update frequency: number of
  iterations between two display of the state of the chain.
  0 is no updates.}

  \item{saveFreq}{how often to save monitored variables to
  a file (in number of iterations)?}

  \item{sdProp}{the initial standard deviation for all
  parameters proposal TODO(Corentin): should be understood
  in rate of the initial value for each parameter and
  default should be set to 0.1 or 10% of the initial value
  for each parameter}

  \item{adaptOK}{should the sampling variance be
  automatically set?  This is very recommended, specially
  with multiple parameters as this step can be painful
  manually and bad proposal variance can lead to poor
  mixing. The sampler really began sampling only once all
  parameters are in the recommended acceptance range}

  \item{checkAdapt}{How often (in iterations) should the
  acceptance be checked and sdProp adapted}

  \item{lowAcceptRate}{minimum acceptance rate (see
  AdaptSdProp)}

  \item{highAcceptRate}{maximum acceptance rate (see
  AdaptSdProp)}

  \item{useAutoStop}{should the chain automatically stop
  when convergence and enough iterations have been
  sampled?}

  \item{checkAutoStop}{Number of iterations before first
  checking if enough iterations (if
  \code{useAutoStop==TRUE}).  Unless you use an insane
  number of monitoring variables, the default should be
  just fine.}

  \item{monitorFile}{The base name of the files recording
  monitored variables and acceptance rate}

  \item{...}{Additional parameters to pass to
  FunctionSample, to customize the way the acceptance is
  handle. In the default case look at OmniSample to set
  heating, repeated sampling etc.}
}
\value{
  data.frame with monitored values, one per column and one
  line per iteration. In addition it can have the
  attributes.

  If there is auto-adjustement of sdprop \itemize{
  \item{lastAdapt}{ last iteration of the adjustement of
  \code{sdprop}} } If there is auto-stopping \itemize{
  \item{lastBurnIn}{ last iteration of the burn-In}
  \item{Kthin}{ adequate thining} } As estimated by the
  Raftery-Lewis criterium (see \code{RafteryLewis})
}
\description{
  Results of the not-run example of MCMC(myData,ModelS)
  saved here as it takes a couple of minutes to run. Also
  needed for the tests.
}
\details{
  This function should perform standard sampling for any
  model when \code{repeatSample=FALSE}. In the synthetic
  likelihood framework it is meant to work in conjonction
  with \code{OmniSample}.

  For the auto-adjustment of the sampling variances check
  AdaptSdProp

  For the auto-stopping behavior check CombinedDiag TODO:
  add CombinedDiag here
}
\examples{
### We will fit the mean and standard deviation
### for a gaussian distribution on a sample of 10 draws
### First with an analytical likelihood
### Second with the synthetic likelihood

# make data
set.seed(1000)
nbDrawsInData <- 10
sample<-rnorm(nbDrawsInData,mean=10,sd=2)

# format the options for the MCMC chain
myData <- list(
               initValues=c(10,2),
               parmNames=c("mean","sd"),
               sampling=c("norm","lnorm"),
               monNames=c("mean","sd","lik"),
               data=sample
               )
# First model: analytical likelihood
ModelA<-function(theta,myData,...){
  LP<-sum(dnorm(myData$data,mean=theta[1],sd=theta[2],log=TRUE))
  return(list(LP=LP,monitor=c(theta,LP)))
}
# Second model: synthetic likelihood
ModelS<-function(theta,myData,...){
  dataSimul<-BasicNSimulator(50,theta,10)
  statsSimul<-BasicNStatsMoments(dataSimul,NULL)
  statsData<-BasicNStatsMoments(myData$data,NULL)
  LP<-synLik(sY=t(statsSimul),sy=t(statsData))
  return(list(LP=LP,monitor=c(theta,LP)))
}

### run the two models
# run the MCMC on analytical model
system.time(monitoredA<-MCMC(myData,ModelA)) # should run in ~5 seconds

# run the MCMC on synthetic likelihood model
# system.time(monitoredS<-MCMC(myData,ModelS)) # should run in ~4 minutes
data(monitoredS)

### normalize marginal distributions to allow comparison
# normalize the marginal distributions
outMeanA<-reweight.marginal.distributions(monitoredA[,1],exp(monitoredA[,3]))
outSdA<-reweight.marginal.distributions(monitoredA[,2],exp(monitoredA[,3]))
outMeanS<-reweight.marginal.distributions(monitoredS[,1],exp(monitoredS[,3]))
outSdS<-reweight.marginal.distributions(monitoredS[,2],exp(monitoredS[,3]))

### plot to compare
par(mfrow=c(1,2))
# mean
plot(outMeanA$param,outMeanA$lik,main="",xlab="Mean",
   ylab="Density",cex=0.2)
lines(outMeanS$param,outMeanS$lik,type="p",col="blue",cex=0.2)
abline(v=mean(sample))

# sd
plot(outSdA$param,outSdA$lik,main="",xlab="Standard deviation",
   ylab="Density",cex=0.2)
lines(outSdS$param,outSdS$lik,type="p",col="blue",cex=0.2)
abline(v=sd(sample))
legend("topright",pch=c(19,19,-1),lty=c(0,0,1),
   c("analytic","synthetic","sample"),col=c("black","blue","black"))
}
\seealso{
  OmniSample,AdaptSdProp
}
\keyword{data}
\keyword{datasets}

