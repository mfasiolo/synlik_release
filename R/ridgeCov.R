# 
# library("Rcpp")
# 
# cppFunction('
#             
# List ridgeCov(NumericMatrix Y, NumericMatrix T, double lambda)
# {
#  int dim = Y.ncol();
#  int nObs = Y.nrow();
# 
#  NumericVector D(dim);
# 
# 
#  double pred;
#  int theEnd;
# 
#  for(int iD = 0; iD < dim; iD++)
#  {
# 
#   for(int iOb = 0; iOb < nObs; iOb++)
#   {
#   pred = 0;
# 
#   for(int ii = 0; ii < iD; ii++) pred += Y(iOb, ii) * T(ii, iD);
# 
#   D[iD] += pow( Y(iOb, iD) - pred , 2.0);
#   }
# 
#   D[iD] /= nObs;
# 
#  }
# 
# 
#  return List::create(_["T"] = T, _["D"] = D);
# 
# }
#             
#             
# ')
# 
# 
# library(MASS)
# nDim <- 3
# nObs <- 500
# Y <- mvrnorm(n = nObs, rep(0, 3), matrix(c(10, 3, 0, 3, 2, 0, 0, 0, 1), 3, 3))# matrix(rnorm(nDim * nObs), nObs, nDim) # 
# Y <- t( t(Y) - colMeans(Y))
# 
# phi <- oldphi <- matrix(0, nDim, nDim)
# for(ii in 1:100)
# {
#   tmp <- ridgeCov(Y, T = phi, lambda = 0)
#   phi <- second(Y = Y, D = tmp$D, lambda = 0)
#   
#   if( sum(abs(phi - oldphi)) < 1e-15 ) break
#   oldphi <- phi
# }
# 
# C <- t(phi) %*% diag(tmp$D) %*% phi
# C - crossprod(Y, Y) / nObs
# 
# D <- diag(diag(C)^-0.5, nDim)
# mean(D%*%C%*%D)
# 
# 
# 
# 
# C
# cov(Y)
# 
# 
# 
# 
# second <- function(Y, D, lambda)
# {
#   nDim = ncol(Y)
#   nObs = nrow(Y)
#   
#   output <- matrix(0, nDim, nDim)
#   
#   for(iD in 2:nDim)
#   {
#     H <- matrix(0, iD-1, iD-1)
#     G <- numeric(iD - 1)
#     
#     tmpY <- Y[ , 1:(iD-1), drop = FALSE]
#     
#     # Calculate H and G
#     for(iO in 1:nObs) {
#       H <- H + tcrossprod( tmpY[iO, ] )
#       G <- G + Y[iO, iD] * tmpY[iO, ]
#     }
#     H <- H / D[iD] + lambda * diag(iD - 1)
#     
#     G <- G / D[iD]
#     
#     output[1:(iD - 1), iD] <- solve( H, G ) 
#   }
#   
#   diag(output) <- 1
#   
#   return( output )
# }
#   
# 
# 
