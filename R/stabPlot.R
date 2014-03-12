
stabPlot <- function(ps = FALSE, ...)
{
  
  ## routine to reproduce N&G blowfly stability plot
  
  f <- function(d,P,mu,omega) {
    ## stability equations cf. p 296 Nisbet and Gurney
    f <- rep(0,2);J <- matrix(0,2,2)
    lPd1 <- log(P/d)-1;emu <- exp(mu)
    f[1] <- d*(1 + lPd1*emu*cos(omega)) - mu
    f[2] <- d*lPd1*emu*sin(omega) - omega
    J <- matrix(0,2,2) ## df_i/d_theta_j
    J[1,1] <- d*lPd1*emu * cos(omega) - 1 ## df_1/d mu
    J[1,2] <- - d * lPd1 * emu * sin(omega) ## df_1/d omega
    J[2,1] <- d * lPd1 * emu * sin(omega)
    J[2,2] <- d * lPd1 * emu * cos(omega) - 1
    list(f=f,J=J)
  }
  
  mu.omega <- function(d,P,mu=0,omega=0.1) {
    ## NOTE: assumes P>d --- otherwise not feasible.
    ## Solves equations for mu and omega, given d (delta tau)
    ## and P (P tau).
    
    x0 <- c(mu,omega)
    F <- f(d,P,x0[1],x0[2])
    fnorm0 <- sum(F$f^2)
    for (i in 1:100) {
      dx <- solve(F$J,-F$f)
      while (1) { 
        x <- x0 + dx
        F <- f(d,P,x[1],x[2])
        fnorm <- sum(F$f^2) 
        if (fnorm<1e-24) return(list(mu=x[1],omega=x[2]))
        if (fnorm < fnorm0) break 
        dx <- dx/2
      }
      if (sum(dx^2)<1e-16) break
      x0 <- x;fnorm0 <- fnorm 
    } ## end of Newton Raphson loop
    list(mu=x[1],omega=x[2])
  }
  
  ## search for the over/under-damped boundary....
  
  objP <- function(P,d) {
    mu.omega(d,P)$omega - 1e-6  
  }
  
  objd <- function(P,d) {
    mu.omega(d,P)$mu  
  }
  
  dT <- exp(seq(log(.01),log(13),length=200)) 
  
  PT <- NA * dT;lims <- c(1,2000)
  for (i in 37:200) {
    if (i>37) { ## bracket
      lims <- rep(PT[i-1]*.99,2*1.01)
      while (objP(lims[1],dT[i])>0) lims[1] <- lims[1]*.97
      while (objP(lims[2],dT[i])<0) lims[2] <- lims[2]*1.03
    }
    PT[i] <- uniroot(objP,lims,d=dT[i])$root
  }
  
  ## Now search for the stable/unstable boundary
  
  PT1 <- NA * dT
  lims <- c(10,2000)
  for (i in 86:200) {
    if (i>86) { ## bracket
      lims <- rep(PT[i-1]*.99,2*1.01)
      while (objd(lims[1],dT[i])<0) lims[1] <- lims[1]*.97
      while (objd(lims[2],dT[i])>0) lims[2] <- lims[2]*1.03
    }
    PT1[i] <- uniroot(objd,lims,d=dT[i],tol=1e-8)$root
  }
  
  start <- 30000;stop <- 50000
  
  #par(mfrow=c(1,1),mar=c(5,5,1,1))
  ylim <- c(1,1000);xlim=c(.01,20)
  #ylim <- c(20,1000);xlim=c(.2,10)
  #ylim <- c(5,500);xlim=c(.1,10)
  if (ps) {
    cex1 <- 1.7;cexl <- 1.8
    cexa <- 1.4;cex2 <- .3
  } else cex1 <- cex2 <- cexa <- cexl <- 1
  plot(dT,PT,type="l",log="xy",ylim=ylim,xlim=xlim,lty=2,
       xlab=expression(paste(delta,tau)),
       ylab=expression(paste(P,tau)),
       cex.axis = cexa,cex.lab=cexl, ...)
  lines(dT,PT1,lty=2) ## any solutions boundary
  abline(0,1,lty=2)
  points(1,40,cex=cex1)
  #points(5,30,pch=19,cex=cex1)
  
  text(4.4,1.1,"no equilibrium",cex=cex1)
  text(.03,3,"overdamped",cex=cex1)
  text(.4,10,"underdamped",cex=cex1)
  text(1,500,"cyclic",cex=cex1)
}