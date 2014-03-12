\name{stork}
\alias{stork}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stochastic Ricker Example}
\description{
Full worked example with completely stochastic ricker.
} %- end description


\details{ Model is N[t+1] ~ Poi(r N[t] exp(- a N[t])). Example shows complete run.

}

\references{


\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}




\examples{
 
stork <- function(theta,burn.in=50,n.t=50,reps=500) {
## Ricker with demographic stochasticity
  r <- theta[1]
  a <- theta[2]
  N = rep(10,reps)
  for (i in 1:burn.in) {
    EB <- N * r * exp(-a * N) 
    B <- rpois(EB,lambda=EB)
    N <-  B 
  } 
  n <- matrix(NA,reps,n.t)
  for (i in 1:n.t) {
    EB <- N * r * exp(-a * N) 
    B <- rpois(EB,lambda=EB)    ## births
    N <-  B 
    n[,i] <- N 
  } 
  t(n)
}

stork.ll <- function(ltheta,y,stats=FALSE,trans=NULL) {
## log likelihood for stochastic Ricker above
  theta <- exp(ltheta)
  n.t <- length(y)
  Y <- stork(theta,burn.in=20,n.t=n.t,reps=600)
  Y[Y==0] <- NA ## replicates are in columns

  ## get run lengths...

  run.length <- colSums(!is.na(Y))

  ## junk reps that never existed...
  
  Y <- Y[,run.length>9]
  if (ncol(Y)>500) Y <- Y[,1:500]  
  
  ## gather the statistics...
  ## nlar and sl.acf have buggy NA handling...
  if (!is.matrix(y)) y <- matrix(y,length(y),1)

  s0 <- nlar(y,lag=c(1,1),power=c(1,2))
  S0 <- nlar(Y,lag=c(1,1),power=c(1,2))
  
  s1 <- sl.acf(y)
  S1 <- sl.acf(Y)
  S1[is.na(S1)] <- 0

  sy <- c(as.numeric(s0),as.numeric(s1),
          mean(y),mean(y)-median(y))
  
  cmeanY <- colMeans(Y,na.rm=TRUE) 
  cmedianY <- apply(Y,2,median,na.rm=TRUE)
  sY <- rbind(S0,S1,cmeanY,cmeanY-cmedianY)
  
  cs <- colSums(sY)
  sY <- sY[,is.finite(cs)] ## ditch any col with non finite stats

  if (!is.null(trans)) {
    sy <- trans.stat(sy,trans)
    sY <- trans.stat(sY,trans)
  }

  if (stats) {
    attr(sY,"observed") <- sy
    return(sY)
  }

  er <- robust.vcov(sY)

  rss <- sum((er$E\%*\%(sy-er$mY))^2)
  ll <- -rss/2 - er$half.ldet.V 
  ll
}


library(sl)
set.seed(1)

## simulate data....
theta <- c(3.8,.18) 
names(theta) <- c("r","a")
Y <- stork(theta,burn.in=20,n.t=100)

sum(rowSums(Y==0)>0) ## number that went extinct

##plot reps
i <- 1
plot(Y[,i],type="l");i<-i+1 ## repeat this line to see more

y <- Y[,1] ## simulated data vector

## show that a quadratic AR carries information...

y1 <- y[2:100]
y0 <- y[1:99]
b <- lm(y1~y0+I(y0^2))
plot(y0,y1);points(y0,fitted(b),col=2,pch=19)

## check Normality....

sY <- stork.ll(log(theta),y,stats=TRUE)
sy <- attr(sY,"observed")
par(mfrow=c(2,2))
MVN.check(sY,sy)


########################################################
## The MCMC chain exploring the synthetic log likelihood
########################################################

  srktrans <- NULL ## default value (stats used raw)

## load transformations of statistics obtained after a pilot run with 
## untransformed statistics. Comment out this line to use 
## statistics untransformed....
  
 # data(srktrans) 

  set.seed(3) ## comment out for replicates!!
  n.mc <- 600 ## just used for checking
  \dontrun{
  n.mc <- 50000 ## length of MCMC run
  }
  theta <- c(4,0.1) ## gets dynamics somewhere in right ball park
  n.para <- length(theta)
  th <- matrix(0,n.para,n.mc) ## storage for chain results
  llr <- rep(NA,n.mc)         ## storage for log synthetic likelihood
  th[,1] <- log(theta)        ## initial state
  prop.sd <- c(.1,.1)    ## proposal standard deviations

 
  llr[1] <- stork.ll(th[,1],y,trans=srktrans)
  reject <- 0
  uni <- runif(n.mc)
  pname <- c("r","a")
  for (i in 2:n.mc) { ## MCMC loop
    th[,i] <- th[,i-1] + rnorm(n.para)*prop.sd
  
    llr[i] <- stork.ll(th[,i],y,trans=srktrans)
  
    alpha <- min(1,exp(llr[i]-llr[i-1]))

    if (uni[i]>alpha) { ## reject
      th[,i] <- th[,i-1]
      llr[i] <- llr[i-1]
      reject <- reject + 1
    }

    if (i\%\%200==0) { ## plot the chains so far (and report acceptance rate)
       par(mfrow=c(2,2))
       start <- 1
       for (j in 1:n.para) { 
         xx <- exp(th[j,1:i])
         if (i < 5000) ylim=range(xx) else ylim=range(xx[(i/2):i])
         plot(1:i,xx,type="l",ylab=pname[j],ylim=ylim)
       }
       xx <- llr[1:i]
       if (i < 5000) ylim=range(xx) else ylim=range(xx[(i/2):i])
       plot(1:i,xx,type="l",ylab="log l_s",ylim=ylim,main=round(1-reject/i,digits=2))
    }
  } ## end of the MCMC chain

  1-reject/n.mc ## acceptance rate
  
  ## checking MVN etc...
  sY <- stork.ll(th[,n.mc],y,stats=TRUE)
  sy <- attr(sY,"observed")
  par(mfrow=c(2,2))
  MVN.check(sY,sy)


}


\keyword{models}%-- one or more ..





