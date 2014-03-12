# Get a correlation matrix out of a covariance matrix
extractCorr <- function(mat)
{
  diag(sqrt( diag(mat) )^-1, nrow(mat)) %*% mat %*% diag(sqrt( diag(mat) )^-1, nrow(mat))  
}
