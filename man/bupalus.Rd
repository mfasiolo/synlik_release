\name{bupalus}
\alias{bupalus}
\alias{bup.para}
\alias{bupalus.ll}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Fit model to data on bupalus}
\description{
A stochastic host parasitoid model for Bupalus. See example section for fitting code.
} %- end description

\usage{
bup.para(theta,e,u,burn.in)
bupalus.ll(theta,y,e,u,burn.in,trans=NULL,stats=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{theta}{parameter vector containing log.r, sig.e and phi.}
\item{y}{population data series.}
\item{e}{matrix of process noise/random effect deviates. Number of columns gives 
         number of replicate series. Number of rows gives total length of each run, of which
         first \code{burn.in} is discarded.}
\item{u}{vector of uniform random deviates, used to generate data from expected data. Unused
         at present}
\item{burn.in}{Number of steps to discard from start of iteration.}
\item{trans}{an object returned by \code{\link{get.trans}} for transforming statistics to normality. \code{NULL} for no transformation.}
\item{stats}{Return statistics, not the log synthetic likelihood.}
}

\details{ The \code{tentsmuir} data frame comes from file \code{tentspar.txt}. The transformed data in \code{tpop} are 
population measurements to the power 0.18. Almost identical data are in \code{tens2.dat} where the transormed 
population data reported is Box-cox 0.25.


}
\value{matrix of replicate data series
}
\references{


\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  
ricker, blowfly, stork
}

\examples{
  library(sl)
  data(tentsmuir) ## load bupalus at Tentsmuir data
  theta <- c(2.5,20,.6,.6,.2,.1)

  pname <- names(theta) <- c("r0","K","a","b","sig.e","sigma")
  burn.in <- 50
  n.t <- 38
  n.reps <- 500  
  
  ## deviates to use in simulation....
  e <- matrix(rnorm((burn.in+n.t)*n.reps),burn.in+n.t,n.reps)
  u <- runif(n.t*n.reps)

  ## Example runs...
  Y <- bup.para(theta,e,u,burn.in)

  par(mfrow=c(2,1),mar=c(4,4,1,1))
  with(tentsmuir,plot(year,pop^.25,type="l")) ## data
  plot(1:38,Y[,1],type="l")                   ## model


  ## Initial MVN checking...
  bpop <- tentsmuir$pop^.25
  sY <- bupalus.ll(theta,bpop,e,u,burn.in,stats=TRUE)
  sy <- attr(sY,"observed")
  par(mfrow=c(2,2))
  MVN.check(sY,sy)


  ## section the log likelihood
  lo <- c(1.5,5,.2,.2,.01,.01)
  hi <- c(4,50,1,1,.4,.4)

  k <- 6 ## which parameter to look at
  n.s <- 100
  thetap <- theta
  ll <- xx <- seq(lo[k],hi[k],length=n.s)
  for (i in 1:n.s) {
    thetap[k] <- xx[i]
    ll[i] <- bupalus.ll(thetap,bpop,e,u,burn.in)
  }
  plot(xx,ll,type="l")

  #######################
  ## The MCMC chain...
  #######################

  n.mc <- 600 ## for package checking
  \dontrun{
  n.mc <- 50000 ## for serious analysis
  }

  th <- matrix(0,6,n.mc)
  llr <- rep(NA,n.mc)
  th[,1] <- log(theta)
  prop.sd <- c(.1,.1,.1,.1,.1,.1)/2
  ll <- c(-100,log(2),-100,-100,-100,log(.05))
  ul <- c(100,log(100),100,100,100,log(.5))

  butrans <- NULL
  ## load transforms to normality for bupalus. Comment out
  ## for no transform....
  data(butrans)

  llr[1] <- bupalus.ll(exp(th[,1]),bpop,e,u,burn.in,trans=butrans)
  reject <- 0
 
  uni <- runif(n.mc)
  for (i in 2:n.mc) { ## MCMC loop
    th[,i] <- th[,i-1] + rnorm(6)*prop.sd
   
    llr[i] <- bupalus.ll(exp(th[,i]),bpop,e,u,burn.in,trans=butrans)
  
    alpha <- min(1,exp(llr[i]-llr[i-1]))
    if (sum(th[,i]<ll|th[,i]>ul)) alpha <- 0
   
    if (uni[i]>alpha) { ## reject
      th[,i] <- th[,i-1]
      llr[i] <- llr[i-1]
      reject <- reject + 1
    }
    
    if (i\%\%200==0) { ## plot progress 
       par(mfrow=c(2,4))
       for (j in 1:6) {
         xx <- exp(th[j,1:i]) 
         if (i<5000) ylim=range(xx) else ylim=range(xx[(i/2):i])
         plot(1:i,xx,type="l",ylim=ylim,ylab=pname[j])
       }
       xx <- llr[1:i]
       if (i<5000) ylim=range(xx) else ylim=range(xx[(i/2):i])
       plot(1:i,xx,type="l",ylim=ylim,ylab="log l_s",main=round(1-reject/i,digits=2))

    }
  } ## end of MCMC loop

  1-reject/n.mc ## acceptance rate


  ## histograms of parameter dist from chain
  par(mfrow=c(2,3))
  start <- round(n.mc/4)
  ind <- start:n.mc
  for (j in 1:6) hist(exp(th[j,ind]),xlab=pname[j])


  ## some example reps....
  plot(1:38,bpop,type="l",ylim=range(bpop))
  
  for (j in 0:4) {
    k <- n.mc-2000*j
    Y <- bup.para(exp(th[,k]),e,u,burn.in)
    plot(1:38,Y[,j+30],type="l",ylim=range(bpop))
  }

}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





