not.sq <- function(x,alpha= .5,d0 = 5) {
  ## not.sq(x) = x^2             if x <= d0
  ##           = k*x^alpha + c      otherwise
  ## k and c are chosen to give continuity to first derivative
  ind <- x < d0
  x[ind] <- x[ind]^2
  k <- 2*d0^(2-alpha)/alpha
  cc <- d0^2 - k * d0^alpha
  x[!ind] <- k * x[!ind]^ alpha + cc  
  x
}
