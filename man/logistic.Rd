\name{logistic}
\alias{logistic}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Fit logistic and Ricker models to same data.}
\description{
Example implementing synthetic likelihood for logistic in R, and comparing fit to Ricker. 
} %- end description



\details{ See example code, noting that a robust version of the synthetic likelihood is used
}

\references{
\url{http://www.maths.bath.ac.uk/~sw283/}
}
\author{ Simon N. Wood \email{simon.wood@r-project.org}}


\seealso{  
not.sq, ricker, stork
}

\examples{

lg.map <- function(theta,e,u,burn.in=50) {
## A stochastic version of the logistic map, 
## model is n_t+1 = r n_t (1 - n_t), r~U(b-a,b), a<b<4,
## y_t ~ Poi(\phi n_t)

  a <- theta[1]
  b <- theta[2]
  phi <- theta[3]
  n.t <- nrow(e) - burn.in
  n.reps <- ncol(e)
  n <- matrix(0, n.t, n.reps)
  N <- rep(.2,n.reps)
  for (i in 1:burn.in) {
    r <- b - e[i,]*a
    N <- r*N*(1-N)
  }
  for (i in 1:n.t) {
    r <- b -  e[i,]*a
    N <- r*N*(1-N)
    n[i,] <- N
  }  
  y <- rpois(phi*n,phi*n)
  matrix(y, n.t, n.reps)
}


logistic.ll <- function(theta,y,e,u,burn.in,trans=NULL,stats=FALSE,ricker=FALSE) {
## function to obtain log synthetic likelihood for the stochastic logistic model
## as coded in `lg.map', or for the stochastic ricker

## simulate from model with current parameter values
  if (ricker) Y <- ricker(theta, e, u, burn.in) else
  Y <- lg.map(theta,e,u,burn.in)

## Now assemble the relevant statistics
  if (!is.matrix(y)) y <- matrix(y,length(y),1)
  acf.Y <- sl.acf(Y,max.lag=5)
  acf.y <- sl.acf(y,max.lag=5)

  b0.Y <- nlar(Y^.3,lag=c(1,1),power=c(1,2))
  b0.y <- nlar(y^.3,lag=c(1,1),power=c(1,2))

  b1.Y <- order.dist(Y,y,np=3,diff=1)
  b1.y <- order.dist(y,y,np=3,diff=1)   

## combine the statistics...

  sy <- c(as.numeric(acf.y),
          as.numeric(b0.y),
          as.numeric(b1.y),
          mean(y)#,sum(y==0)
         )
  sY <- rbind(acf.Y,
              b0.Y, 
              b1.Y,
              colMeans(Y)#,colSums(Y==0)
             )

  if (!is.null(trans)) {
    sy <- trans.stat(sy,trans)
    sY <- trans.stat(sY,trans)
  }

## get the log synthetic likelihood

  sY <- sY[,is.finite(colSums(sY))]
  if (stats) {  ## return statistics
   attr(sY,"observed") <- sy
   return(sY) 
  }

  er <- robust.vcov(sY)

  ## robustify the likelihood...
  rss <- sum((er$E\%*\%(sy-er$mY))^2)
  
  ll0 <- -rss/2 - er$half.ldet.V ## true l_s

  d0 <- qchisq(.99,nrow(sY))^.5

  rss <- not.sq(sqrt(rss),alpha=.1,d0=d0)
 
  ll <- -rss/2 - er$half.ldet.V ## robustified l_s
  attr(ll,"true") <- ll0 ## append the true l_s
  ll
  
} ## end of logistic.ll

library(sl)
set.seed(3)
n.t <- 50 ## length of data run
burn.in <- 50 ## burnin length for model
n.reps <- 500 ## slow mixing if this lowered 

## simulate random deviates to re-use for fitting...
e <- matrix(runif((burn.in+n.t)*n.reps),burn.in+n.t,n.reps)
u <- runif(n.t*n.reps) 

## simulate from stochastic logistic....

theta <- c(.3,3.3,100)
names(theta) <- c("a","b","phi")
Y <- lg.map(theta,e,u,burn.in=burn.in)
yl <- y <- Y[,1]

## ... or simulate from stochastic ricker...

theta <- c(3.2,0.3,10)
names(theta) <- c("r","sig.e","phi")
Y <- ricker(theta,e,u,burn.in)
yr <- y <- Y[,1]



###################
# MCMC chain....
###################

  set.seed(3) ## comment out for replicates!!
  ricker <- FALSE ## fit ricker model....
  if (ricker) {
    theta <- c(3,0.3,10)
    pname <- names(theta) <- c("r","sig.e","phi")
  } else {
    theta <- c(.2,3,80)
    pname <- names(theta) <- c("a","b","phi")
  }
  rktrans <- NULL
  n.mc <- 600 ## just for checking

\dontrun{ 
  n.mc <- 50000 ## length of MCMC run
}
  n.para <- length(theta)
  th <- matrix(0,n.para,n.mc) ## storage for chain results
  llr <- rep(NA,n.mc)         ## storage for log synthetic likelihood
  th[,1] <- log(theta)        ## initial state
  if (ricker) {
    prop.sd <- c(.02,.1,.05) 
  } else {
    prop.sd <- c(.02,.01,.01)*2    ## proposal standard deviations
  }
 
  llt <- logistic.ll(exp(th[,1]),y,e,u,trans=rktrans,burn.in=burn.in,ricker=ricker)
  
  llr[1] <- attr(llt,"true")
  reject <- 0
  uni <- runif(n.mc)
  
  for (i in 2:n.mc) { ## MCMC loop
    
    ok <- FALSE
    while (!ok) {
      th[,i] <- th[,i-1] + rnorm(n.para)*prop.sd
      pth <- exp(th[,i])
      if (ricker||(pth[1]<pth[2]&&pth[2]<4)) ok <- TRUE
    }  

    llp <- logistic.ll(pth,y,e,u,trans=rktrans,burn.in=burn.in,ricker=ricker)
  
    alpha <- min(1,exp(llp - llt))
    if (uni[i]>alpha) { ## reject
      th[,i] <- th[,i-1]
      llr[i] <- llr[i-1]
      reject <- reject + 1
    } else { ## accept
      llr[i] <- attr(llp,"true")
      llt <- llp
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

  rownames(th) <- pname
  thl <- rbind(th,llr)
  start <- round(n.mc*.4)
  ml <- chain2ll(thl,ll="llr",start=start)
  
  
  sY <- logistic.ll(pth,y,e,u,trans=rktrans,burn.in=burn.in,ricker=ricker,stats=TRUE)
  sy <- attr(sY,"observed")
  par(mfrow=c(2,2))
  MVN.check(sY,sy)

}


\keyword{models}




