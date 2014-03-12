#####
#' Empirical cumulant generating function
#' @description Calculates the empirical cumulant generating function (CGF) and its derivatives
#'               given a sample of n d-dimentional vectors
#'
#' @param lambda Point at which the CGF is evaluated (d-dimensional vector).
#' @param X (n by d) matrix containing the data.
#' @param grad If grad == 0 only the value of the CGF is returned, 
#'             if grad == 1 also its first derivative wrt lambda 
#'             and if grad == 2 also the second derivarive wrt lambda.
#' @param mix Mixture of empirical and normal CGF to use (if 1 only empirical CGF is used).
#'            Default value is 0.9.
#' @return A list with entries:
#'         \itemize{
#'         \item{ \code{K} }{The value of the empirical CGF at lambda;}
#'         \item{ \code{dK} }{The value of the gradient empirical CGF wrt lambda at lambda;}
#'         \item{ \code{d2K} }{The value of the hessian of the empirical CGF wrt lambda at lambda;}
#'         }
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon Wood.
#' @examples 
#' X <- matrix(rnorm(2 * 1e3), 1e3, 2)
#' ecgf(lambda = c(0, 0), X = X) 
#' @export
#'
ecgf <- function(lambda, X, grad = 0, mix = 0.9) {
  ## X[i,j] is ith rep of jth variable. Evaluate observed KGF 
  ## and its derivs w.r.t. lambda, without overflow...
  
  .ecgf(lambda = lambda, 
        X = X, 
        kum1 = colMeans(X), 
        kum2 = .robCov(t(X), alpha2 = 4, beta2 = 1.25)$COV, 
        grad = grad, 
        mix = mix, 
        deriv = FALSE,  
        addList = list(NaN) )
  
}