## routine's for transformation of statistics to better meet normality assumptions,
## and for checking the MVN approximation

getPieceTrans <- function(object,
                          param = object@param, 
                          nsim = 1e3,
                          ...) {
  ## object is an ns by n.reps matrix of statistics. This routine works through
  ## its rows finding piecewise linear transformations to normality, by 
  ## interactive use of `locator'.
  
  if( is(object, "synlik") )
  {
    S <- t( simulate(object, nsim, param, stats = TRUE, ...) )  
  } else{
    if( !is.numeric(object) ) stop("object should be either of class \"synlik\", a numeric matrix or vector")
    if( !is.matrix(object) ) object <- matrix(object, length(object), 1)
    
    if( nrow(object) < ncol(object) ) warning("The number of variables is higher then the number of simulations (nsim) used")
    S <- t( object )
  }
  
  op <- par(mfrow=c(1,1))
  if (!is.matrix(S)) S <- matrix(S,1,length(S))
  ns <- nrow(S)    ## the number of statistics
  n.rep <- ncol(S) ## the number of replicates
  z <- qnorm((1:n.rep-0.5)/n.rep) ## standard normal quantiles
  trans <- list()
  for (i in 1:ns) { ## loop through stats...
    plot(sort(S[i,]),z)
    tr <- locator( , type="l", col=2)
    if (length(tr$x)<2) { 
      warning(paste("no transform for variable", i));
      trans[[i]] <- list(x = NA, y = NA)  
    } else
      ## now extend end segments
      if (length(tr$x)==2) { ## single segment --- extend both ends
        xr <- tr$x[2]-tr$x[1]
        slope <- (tr$y[2]-tr$y[1])/xr
        tr$x[1] <- tr$x[1] - 1000*xr 
        tr$y[1] <- tr$y[1] - slope*1000*xr
        tr$x[2] <- tr$x[2] + 1000*xr
        tr$y[2] <- tr$y[2] + slope*1000*xr
        trans[[i]] <- tr      
      } else { ## extend end segments
        xr <- max(tr$x) - min(tr$x)
        slope <- (tr$y[2]-tr$y[1])/(tr$x[2]-tr$x[1])
        tr$x[1] <- tr$x[1] - 1000*xr 
        tr$y[1] <- tr$y[1] - slope*1000*xr
        nt <- length(tr$x)
        slope <- (tr$y[nt]-tr$y[nt-1])/(tr$x[nt]-tr$x[nt-1])
        tr$x[nt] <- tr$x[nt] + 1000*xr
        tr$y[nt] <- tr$y[nt] + slope*1000*xr
        trans[[i]] <- tr     
      }
  }
 # trans[[ns+1]] <- NA
  par(op)
  
  rm(S)
  
  transFun <- function(S, ...) {
    ## apply a piecewise linear `trans' object to the rows 
    ## of a statistics matrix, S
    if (!is.matrix(S)) S <- matrix(S, length(S), 1)
    
    for (i in 1:ncol(S)) {
      if (!is.na(trans[[i]]$x[1]))
        S[ , i] <- approx(trans[[i]]$x, trans[[i]]$y, S[ , i], rule=2)$y
    }
    if (ncol(S)==1) S <- as.numeric(S)
    S
  }
  
  return(transFun)
}