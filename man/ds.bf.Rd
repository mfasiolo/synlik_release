\name{ds.bf}
\alias{ds.bf}
\alias{dsbf.ll}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Fit model to Nicholson's Blowfly data}
\description{
A Nisbet and Gurney type model for Nicholson's blowfly data, where the only stochasticity is demographic. 
See example for fitting code.
} %- end description

\usage{
ds.bf(theta,burn.in=100,n.t=100,n.rep=500)
dsbf.ll(theta,y,burn.in,n.rep=500,trans=NULL,stats=FALSE,step=2)
}
%- maybe also `usage' for other objects documented here.
\arguments{

\item{theta}{parameter vector containing delta, P, N0 and tau.}
\item{y}{observed population data.}
\item{burn.in}{Number of steps to discard from start of iteration.}
\item{n.t}{number of days to simulate, post burn in.}
\item{n.rep}{number of replicates}
\item{trans}{an object returned by \code{\link{get.trans}} for transforming statistics to normality. \code{NULL} for no transformation.}
\item{stats}{Return statistics, not the log synthetic likelihood.}
\item{step}{time step in data.}
}

\details{ In this model births are a Poisson proccess and death is a binomial proccess. This is the only stochasticity. 


The data in the \code{blowfly} dataframe are population I from fig 8 of Nicholson (1957),
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
  
  theta <- c(.16,6.5,400,14) 
  pname <- names(theta) <- c("delta","P","N0","tau")
  
  burn.in <- 200;n.reps <- 500  
 
  ## an initial run of the blowfly model...
  Y <- ds.bf(theta,burn.in,n.y*step,n.reps)[1:n.y*step,]
  
  ## plot data and some model replicates...
  par(mfrow=c(2,1),mar=c(5,5,1,1))
  with(bf,plot(day,pop,type="l"))
  plot(1:n.y*step,Y[,1],type="l",xlab="day")
  for (i in 2:5) lines(1:n.y*step,Y[,i],col=i)

  ## initial MVN check of statistics....
  par(mfrow=c(2,2))
  sY <- dsbf.ll(theta,bf$pop,burn.in,n.rep=n.reps,stats=TRUE,step=step)
  sy <- attr(sY,"observed")
  MVN.check(sY,sy)  


  ## section the log likelihood....
  lo <- c(.1,3,200,10)
  hi <- c(.2,9,600,16)

  k <- 3 ## which parameter to vary
  n.s <- 100
  thetap <- theta
  ll <- xx <- seq(lo[k],hi[k],length=n.s)
  for (i in 1:n.s) {
    thetap[k] <- xx[i]
    ll[i] <-  dsbf.ll(thetap,bf$pop,burn.in,n.rep=n.reps,step=step) 
  }
  plot(xx,ll,type="l")

  ## get transformation info
  data(bftrans) ## not worthwhile
  bftrans <- NULL 
 
  n.mc <- 600 ## for checking only 
  \dontrun{
  n.mc <- 50000 ## for serious analysis
  }
  th <- matrix(0,4,n.mc)  ## storage for parameters
  rss <- llr <- rep(NA,n.mc) ## storage for l_s and diagnostic

  theta0 <- c(.1,4,200,10) ## ok start (long way from truth)
  # theta0 <- c(.5,3,200,10) ## start in middle of G&N Quasi cycle range  

  th[,1] <- c(log(theta0[1:3]),theta0[4]) ## transform to log scale
  prop.sd <- c(.1,.1,.1)*.2  ## poposal standard deviation
  pind <- c(1:3) ## index params with normal proposal
  thetap <- c(exp(th[1:3,1]),th[4,1]) ## params on original scale
  
  ## initial l_s
  llt <-  dsbf.ll(thetap,bf$pop,burn.in,n.rep=n.reps,trans=bftrans,step=step) 
  llr[1] <- attr(llt,"true")
  reject <- 0
  uni <- runif(n.mc)
  ll <- c(-100,-100,-100,8)
  ul <- c(100,100,100,20)

  for (i in 2:n.mc) { ## main MCMC loop
    ## make proposal....
    th[pind,i] <- th[pind,i-1] + rnorm(3)*prop.sd
    th[4,i] <- th[4,i-1] + sample(c(-1,0,0,1),1)
    thetap <- c(exp(th[1:3,i]),th[4,i])
   
    ## l_s for proposal...
    llp <- dsbf.ll(thetap,bf$pop,burn.in,n.rep=n.reps,trans=bftrans,step=step) 
   
    ## get acceptance probability...
    
    alpha <- min(1,exp(llp-llt))
    if (sum(th[,i]<ll|th[,i]>ul)) alpha <- 0
    
    ## accept/reject...
    if (uni[i]>alpha) { ## reject
      th[,i] <- th[,i-1]
      llr[i] <- llr[i-1]
      reject <- reject + 1
    } else { ## accept
      llr[i] <- attr(llp,"true")
      llt <- llp
    }
    ## Plot results every so often....
    if (i\%\%200==0) {
      ind <- 1:i
      if (i>1000) ind2 <- round(i/2):i else ind2 <- ind ## for plot range 
      par(mfrow=c(2,3),mar=c(4,4,1,1))
      for (k in 1:4) plot(ind,th[k,ind],type="l",ylim=range(th[k,ind2]),ylab=pname[k])
      plot(ind,llr[ind],type="l",ylim=range(llr[ind2]))
      ##plot(ind,rss[ind],type="l",ylim=range(rss[ind2]),main=round(1-reject/i,digits=2))
      thetap <- c(exp(th[1:3,i]),th[4,i])
      Y <- ds.bf(thetap,burn.in,n.y*step,n.reps)[1:n.y*step,]
      lim <- range(bf$pop)
      plot(1:n.y,Y[,1],type="l",ylim=lim)
      lines(1:n.y,bf$pop,col=2)
      lines(1:n.y,Y[,2],type="l",ylim=lim,lty=2)
    }
  }
  1-reject/n.mc ## acceptance rate

  ## MVN checking...
  thetap <- c(exp(th[1:3,n.mc]),th[4,n.mc])
  sY <- dsbf.ll(thetap,bf$pop,burn.in,n.rep=n.reps,trans=bftrans,stats=TRUE,step=step)
  sy <- attr(sY,"observed")
  par(mfrow=c(2,2))
  MVN.check(sY,sy)

  ## plot chains again...
  ind <- 1:n.mc
  par(mfrow=c(2,3),mar=c(4,4,1,1))
  for (k in 1:4) plot(ind,th[k,ind],type="l")
  plot(ind,llr[ind],type="l")

  ## plot examples...
  p <- 1;k <- n.mc
  lim <- range(bf$pop)
  par(mfrow=c(3,2))
  plot(1:n.y,bf$pop,ylim=lim,type="l")
  Y <- ds.bf(theta0,burn.in,n.y*step,n.reps)[1:n.y*step,]
  plot(1:n.y,Y[,1],ylim=lim,type="l")
  thetap <- c(exp(th[1:3,k]),th[4,k])
  Y <-  ds.bf(thetap,burn.in,n.y*step,n.reps)[1:n.y*step,]
  for (i in 1:4) { plot(1:n.y,Y[,p],type="l",ylim=lim);p <- p+1}


}


\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..





