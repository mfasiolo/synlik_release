# rb <- function(x,z) {
#   100*(z-x^2)^2 + (1-x)^2 
# }
# 
# ## b)
# n <- 100
# xm <- seq(-1.5,1.5,length=n)
# zm <- seq(-.5,1.5,length=n)
# f <- outer(xm,zm,rb)
# contour(xm,zm,matrix(f,n,n))
# 
# ## c)
# contour(xm,zm,matrix(log10(f),n,n),levels=(1:10/2))
# 
# ## d) 
# rb.grad <- function(x,z) {
# ## gradient of rosenbrocks function at a single point
#   g <- rep(NA,2)
#   g[2] <- 200*(z-x^2)
#   g[1] <- 400*(x^3-z*x) + 2*(x-1)
#   g
# }
# 
# ## e)
# x0 <- .5; z0 <- 1; eps <- 1e-7
# f0 <- rb(x0,z0)
# g <- g0 <- rb.grad(x0,z0) ## exact gradiant
# ## put FD gradiant in g0
# g0[1] <- (rb(x0+eps,z0)-f0)/eps
# g0[2] <- (rb(x0,z0+eps)-f0)/eps
# g;g0 ## compare
# 
# ## f) 
# rb.hess <- function(x,z) {
# ## Hessian of Rosenbrock's function, at a single point
#   H <- matrix(NA,2,2)
#   H[2,2] <- 200
#   H[1,1] <- 1200*x^2 - 400*z + 2
#   H[2,1] <- H[1,2] <- -400*x
#   H
# }
# 
# ## g)
# H <- H0 <- rb.hess(x0,z0)
# H0[,1] <- (rb.grad(x0+eps,z0)-g)/eps
# H0[,2] <- (rb.grad(x0,z0+eps)-g)/eps
# H;H0
# 
# ## stochastic approximation...
# contour(xm,zm,matrix(log10(f),n,n),levels=(1:10/2))
# n <- 100;x0 <- -1; z0 <- 0;
# nIter <- 100
# H1 <- diag(c(100,100))
# eh <- eigen(H1)
# d <- eh$values;
# yOld <- xOld <- zOld <- matrix(NA, nIter, n)
# 
# 
# for (j in 1:100) {
#   #X <- eh$vectors%*%(sqrt(d)*matrix(rnorm(n*2),2,n))*.002
#   #x <- X[1,]
#   #z <- X[2,]
#   if( j == 1){ cholDec <- chol( eh$vectors %*% t(eh$vectors) / eh$values ) 
#   } else { cholDec <- 0.9 * cholDec + 0.1 * chol( eh$vectors %*% t(eh$vectors) / eh$values ) }
# 
#   tmp <- t(cholDec) %*% matrix(rnorm(2 * n), 2, n)
#   x <- xOld[ , j] <- tmp[1, ]
#   z <- zOld[ , j] <- tmp[2, ]
#   
#   y <- yOld[ , j] <- rb(x+x0,z+z0) + rnorm(n)*2
#   points(x+x0,z+z0,pch=".")
#   dat <- data.frame(y=y,x=x,z=z)
#   ## linear regression for the grad and hessian
#   resamp <- TRUE
#   g1 <- rep(0,2)
#   if (resamp) {
#     evmm <- -1e40
#     for (i in 1:100) {
#       b <- lm(y ~ x + z + I(x*z) + I(x^2) + I(z^2),
#         data=dat[sample(1:n,n,replace=TRUE),])
#       g <- coef(b)[2:3]
#       gnorm <- sum((g-g1)^2)
#       if (i == 1) { 
#         gnorm.best <- gnorm
#         g1 <- g 
#       } else {
#         if (gnorm<gnorm.best) {
#           gnorm.best <- gnorm
#           g1 <- g 
#         }
#       }
#       H1 <- matrix(coef(b)[c(5,4,4,6)]*c(1,.5,.5,1),2,2)*2
#       evm <- min(eigen(H1)$values)
#       if (evm>evmm) {
#         evmm <- evm
#         Hbest <- H1
#       }
#     }
#   } 
#   if (j==1 || !resamp) {
#     b <- lm(y ~ x + z + I(x*z) + I(x^2) + I(z^2),
#             data=dat)
#     g1 <- coef(b)[2:3]
#   }
#   H1 <- matrix(coef(b)[c(5,4,4,6)]*c(1,.5,.5,1),2,2)*2
#   if (resamp) H1 <- Hbest
#   eh <- eigen(H1)
#   d <- eh$values;
#   mind <- abs(max(eh$values))*1e-2
#   d[d<mind] <- mind
#   
#   if(rbinom(1, 1, 0.2)){ g1 <- g1 * 100; }
#   
#   dp <- .getUpdate( g1, eh )
#   
#   #dp <- -solve(H1,g1)
#   lines(c(x0,x0+dp[1]),c(z0,z0+dp[2]),col=3)
#   x0 <- x0 + dp[1]
#   z0 <- z0+ dp[2]
# }
# 
# .getUpdate <- function(grad, currHess, oldHess, quant = 0.8)
# {
#   delta <- -Hess$vectors%*%((t(Hess$vectors)%*%grad)/Hess$values)
#   quad <- drop( -crossprod(delta, grad) )
#   limit <- qchisq(quant, length(grad))
#   
#   if(quad > limit) delta <- delta / sqrt(quad / limit)
#   
#   delta
# }
# 
# .getUpdate(grad = c(1, 2) * 10, Hess = eigen(-diag(c(3, 3))))
# 
# 
# g1;g <- rb.grad(x0,z0);g1
# H1;H<-rb.hess(x0,z0);H
# eigen(H)
# 
# contour(xm,zm,matrix(log10(f),n,n),levels=(1:10/2))
# points(x0,z0)
# dp <- -solve(H,g)
# lines(c(x0,x0+dp[1]),c(z0,z0+dp[2]),col=2)
# dp <- -solve(H1,g1)
# lines(c(x0,x0+dp[1]),c(z0,z0+dp[2]),col=3)
# 
# eh <- eigen(H1);d <- eh$values; U <- eh$vectors
# d[d<1] <- 1
# dp <- -U%*%((t(U)%*%g1)/d)
# lines(c(x0,x0+dp[1]),c(z0,z0+dp[2]),col=4)
# 
# 
# 
# 
# 
