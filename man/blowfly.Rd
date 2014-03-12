\name{blowfly}
\alias{blowfly}
\alias{ng.bf}
\alias{blowfly.ll}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Fit model to Nicholson's Blowfly data}
\description{
A stochastic Nisbet and Gurney type model for Nicholson's blowfly data. See example for fitting code.
} %- end description

\usage{
ng.bf(theta,lu,lu1,burn.in)
blowfly.ll(theta,y,lu,lu1,burn.in,trans=NULL,stats=FALSE,step=2)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{theta}{parameter vector containing delta, P, N0, sig.p, tau, sig.d}
\item{y}{observed population data.}
\item{lu}{matrix of process noise/random effect deviates. Number of columns gives 
         number of replicate series. Number of rows gives total length of each run, of which
         first \code{burn.in} is discarded. The deviates are logits of uniforms. These deviates
         are used to generate the gamma deviates actually required in the model.}
\item{lu1}{as \code{lu}, just with a second set of deviates.}
\item{burn.in}{Number of steps to discard from start of iteration.}
\item{trans}{an object returned by \code{\link{get.trans}} for transforming statistics to normality. \code{NULL} for no transformation.}
\item{stats}{Return statistics, not the log synthetic likelihood.}
\item{step}{time step in data.}
}

\details{ The data in the \code{blowfly} dataframe are population I from fig 8 of Nicholson (1957),
and are tabulated in Brillinger et al. (1980). They are the data used in Kendall et al. (2001). See \code{\link{bf}} 
for other blowfly series.

}
\value{matrix of replicate data series
}
\references{
Brillinger, D. R., J. Guckenheimer, P. Guttorp, and G. Oster. 1980. Empirical modelling of population time series data: the case of age and density dependent vital rates. Lectures on Mathematics in the Life Sciences13:65-90. 

Nicholson, A. J. 1957. The self-adjustment of populations to change. Cold Spring Harbor Symposia on Quantitative Biology22:153-173. 


\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  
bupalus, bf, ricker, stork
}

\examples{
  library(sl)
 
  logistic <- function(x,k=.5) k/(1+exp(-x))
  logit <- function(x,k=.5) log((x/k)/(1-x/k))

  set.seed(3) ## for strict repeatablity -- remove for reps
  
  ## load data: uncomment data set required, below
  data(bf1);bf <- bf1[seq(1,275,by=2),] ## not clear that daily in bf1 is raw 
  # data(blowfly);bf <- blowfly
  # data(bf2);bf <- bf2
  # data(bf3);bf <- bf3
 
  step <- bf$day[2] - bf$day[1] ## step length in data
  n.y <- nrow(bf)
  n.t <- n.y*step
  
  theta <- c(.16,6.5,400,.1,14,.1) 
  pname <- names(theta) <- c("delta","P","N0","sig.p","tau","sig.d")
  
  ## the logit uniform deviates to use in sims...
  burn.in <- 100;n.reps <- 500  
  lu <- matrix(logit(runif((burn.in+n.t)*n.reps),k=1),burn.in+n.t,n.reps)
  lu1 <- matrix(logit(runif((burn.in+n.t)*n.reps),k=1),burn.in+n.t,n.reps)
  
  ## an initial run of the blowfly model...
  Y <- ng.bf(theta,lu,lu1,burn.in)[1:n.y*step,]
  
  ## plot data and some model replicates...
  par(mfrow=c(2,1),mar=c(5,5,1,1))
  with(bf,plot(day,pop,type="l"))
  plot(1:n.y*step,Y[,1],type="l",xlab="day")
  for (i in 2:5) lines(1:n.y*step,Y[,i],col=i)

  ## initial MVN check of statistics....
  par(mfrow=c(2,2))
  sY <- blowfly.ll(theta,bf$pop,lu,lu1,burn.in,stats=TRUE,step=step)
  sy <- attr(sY,"observed")
  MVN.check(sY,sy)  


  ## section the log likelihood....
  lo <- c(.1,3,200,.1,10,0)
  hi <- c(.2,9,600,.4,16,.5)

  k <- 4 ## which parameter to vary
  n.s <- 100
  thetap <- theta
  ll <- xx <- seq(lo[k],hi[k],length=n.s)
  for (i in 1:n.s) {
    thetap[k] <- xx[i]
    ll[i] <- blowfly.ll(thetap,bf$pop,lu,lu1,burn.in,step=step)
  }
  plot(xx,ll,type="l")

  ## get transformation info
  data(bftrans) ## not worthwhile
  bftrans <- NULL 
 
  n.mc <- 600 ## for checking only 
  \dontrun{
  n.mc <- 50000 ## for serious analysis
  }
  th <- matrix(0,6,n.mc)  ## storage for parameters
  rss <- llr <- rep(NA,n.mc) ## storage for l_s and diagnostic

  theta0 <- c(.1,4,200,.1,10,.1) ## ok start (long way from truth)
  # theta0 <- c(.5,3,200,.1,10,.1) ## start in middle of G&N Quasi cycle range  

  th[,1] <- c(log(theta0[1:4]),theta0[5],log(theta0[6])) ## transform to log scale
  prop.sd <- c(.1,.1,.1,.1,.1)*.2  ## poposal standard deviation
  pind <- c(1:4,6) ## index params with normal proposal
  thetap <- c(exp(th[1:4,1]),th[5,1],exp(th[6,1])) ## params on original scale
  
  ## initial l_s
  ll <- blowfly.ll(thetap,bf$pop,lu,lu1,burn.in,trans=bftrans,step=step)
  llr[1] <- ll
  rss[1] <- attr(ll,"rss")
  reject <- 0
  uni <- runif(n.mc)
  ll <- c(-100,-100,-100,-3,8,-10)
  ul <- c(100,100,100,2,20,10)

  for (i in 2:n.mc) { ## main MCMC loop
    ## make proposal....
    th[pind,i] <- th[pind,i-1] + rnorm(5)*prop.sd
    th[5,i] <- th[5,i-1] + sample(c(-1,0,0,1),1)
    thetap <- c(exp(th[1:4,i]),th[5,i],exp(th[6,i]))
   
    ## l_s for proposal...
    ll <- blowfly.ll(thetap,bf$pop,lu,lu1,burn.in,trans=bftrans,step=step)
   
    ## get acceptance probability...
    llr[i] <- ll
    rss[i] <- attr(ll,"rss")  
    alpha <- min(1,exp(llr[i]-llr[i-1]))
    if (sum(th[,i]<ll|th[,i]>ul)) alpha <- 0
    
    ## accept/reject...
    if (uni[i]>alpha) { ## reject
      th[,i] <- th[,i-1]
      llr[i] <- llr[i-1]
      reject <- reject + 1
    }
    ## Plot results every so often....
    if (i\%\%200==0) {
      ind <- 1:i
      if (i>1000) ind2 <- round(i/2):i else ind2 <- ind ## for plot range 
      par(mfrow=c(3,3),mar=c(4,4,1,1))
      for (k in 1:6) plot(ind,th[k,ind],type="l",ylim=range(th[k,ind2]),ylab=pname[k])
      plot(ind,llr[ind],type="l",ylim=range(llr[ind2]))
      plot(ind,rss[ind],type="l",ylim=range(rss[ind2]),main=round(1-reject/i,digits=2))
      thetap <- c(exp(th[1:4,i]),th[5,i],exp(th[6,i]))
      Y <- ng.bf(thetap,lu,lu1,burn.in)[1:n.y*step,]
      lim <- range(bf$pop)
      plot(1:n.y,Y[,1],type="l",ylim=lim)
      lines(1:n.y,bf$pop,col=2)
      lines(1:n.y,Y[,2],type="l",ylim=lim,lty=2)
    }
  }
  1-reject/n.mc ## acceptance rate

  ## MVN checking...
  thetap <- c(exp(th[1:4,n.mc]),th[5,n.mc],exp(th[6,n.mc]))
  sY <- blowfly.ll(thetap,bf$pop,lu,lu1,burn.in,trans=bftrans,stats=TRUE,step=step)
  sy <- attr(sY,"observed")
  par(mfrow=c(2,2))
  MVN.check(sY,sy)

  ## plot chains again...
  ind <- 1:n.mc
  par(mfrow=c(2,4),mar=c(4,4,1,1))
  for (k in 1:6) plot(ind,th[k,ind],type="l")
  plot(ind,llr[ind],type="l")

  ## plot examples...
  p <- 1;k <- n.mc
  lim <- range(bf$pop)
  par(mfrow=c(3,2))
  plot(1:n.y,bf$pop,ylim=lim,type="l")
  Y <- ng.bf(theta0,lu,lu1,burn.in)[1:n.y*step,]
  plot(1:n.y,Y[,1],ylim=lim,type="l")
  thetap <- c(exp(th[1:4,k]),th[5,k],exp(th[6,k]))
  Y <- ng.bf(thetap,lu,lu1,burn.in)[1:n.y*step,]
  for (i in 1:4) { plot(1:n.y,Y[,p],type="l",ylim=lim);p <- p+1}


}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





