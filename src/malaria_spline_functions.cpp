#include "synlik.h"

using namespace Rcpp;

// Description later
inline void my_bspline_internal(       double & output, 
                                 const double & x, 
                                 const int & ii, 
                                 const int & p, 
                                 const NumericVector knots, 
                                 const int & nknots);
                                 
// Description later                                
inline void my_periodic_bspline_basis_eval (      double x,   
                                            const double & period, 
                                            const int & degree, 
                                            const int & nBasis, 
                                                   NumericVector & output);

/*
  * Function that CREATES the spline basis
  *
  * ARGS LIST
  * covariate: vector of the covariate (for ex. time)
  * nBasis:    number of basis functions
  * degree:    degree of the bsplines
  * period:    period of the splines 
  *
  * OUTPUT: matrix (length(covariate) by nBasis) of Bsplines.  
*/
  inline NumericMatrix my_periodic_bspline_basis_creator (const NumericVector & covariate, 
                                                          const int & nBasis, 
                                                          const int & degree, 
                                                          const double & period) 
{
    if (degree < 0)        stop("periodic_bspline_basis error: must have degree >= 0");
    if (nBasis <= 0)       stop("periodic_bspline_basis error: must have nbasis > 0");
    if (nBasis < degree)  stop("periodic_bspline_basis error: must have nbasis >= degree");
    if (period <= 0.0)     stop("periodic_bspline_basis error: must have period > 0");
    
    int totLength = covariate.size();
    
    NumericMatrix output(totLength, nBasis); 
    
    // Vector giving some working memory for my_periodic_bspline_basis_eval
    NumericVector workingVett(nBasis); 
    
    /* Fill up the splines (output) one row at the time */
      for (int ii = 0; ii < totLength; ii++) {
        my_periodic_bspline_basis_eval(covariate[ii], period, degree, nBasis, workingVett);
        output(ii, _) = workingVett;
      }
    
    return output;
  }


/*
  * Function that EVALUATES the spline basis
  *
  * ARGS LIST
  * x:         point at which the spline is evaluated (for ex. time = 3.5), x will be modified so no reference!
  * nBasis:    number of basis functions
  * degree:    degree of the bsplines
  * period:    period of the splines 
  *
  * OUTPUT: vector of length(nBasis) elements where the i-th element is an evaluation
  *         of the i-th basis at point x.  
*/
  inline void my_periodic_bspline_basis_eval (      double x,   
                                                    const double & period, 
                                                    const int & degree, 
                                                    const int & nBasis, 
                                                    NumericVector & output)
{
    int nKnots = nBasis + 2*degree + 1;
    int shift = (degree-1)/2;
    NumericVector knots(nKnots); 
    NumericVector tmpWork(nKnots);
    double dx;
    
    if (period <= 0.0)    stop("my_periodic_bspline_basis_eval error: must have period > 0");
    if (nBasis <= 0)      stop("my_periodic_bspline_basis_eval error: must have nBasis > 0");
    if (degree < 0)       stop("my_periodic_bspline_basis_eval error: must have degree >= 0");
    if (nBasis < degree)  stop("my_periodic_bspline_basis_eval error: must have nBasis >= degree");
    
    dx = period/((double) nBasis);
    
    int kk;
    for ( kk = -degree; kk <= nBasis+degree; kk++) {
      knots[degree+kk] = kk*dx;
    }
    
    x = fmod(x,period);
    if (x < 0.0) x += period;
    
    for (kk = 0; kk < nKnots; kk++) {
      my_bspline_internal(tmpWork[kk], x, kk, degree, knots, nKnots);
    }
    
    for (kk = 0; kk < degree; kk++) tmpWork[kk] += tmpWork[nBasis+kk];
    
    int jj;
    for (kk = 0; kk < nBasis; kk++) {
      jj = (shift+kk) % nBasis;
      output[kk] = tmpWork[jj];
    }
  }

/*
  * Function to be called by my_periodic_bspline_basis_eval spline basis
  *
  * ARGS LIST
  * output:    output of whatever this function does
  * x:         I don't know
  * .....      I really don't know what this function does...
  *
  * OUTPUT: the double in "output".  
*/
  
  inline void my_bspline_internal(      double & output, 
                                        const double & x, 
                                        const int & ii, 
                                        const int & p, 
                                        const NumericVector knots, 
                                        const int & nknots)
{
    double a, b;
    double y1, y2;
    
    int ii2, p2;
    
    if (p == 0) {
      output = (double) ((knots[ii] <= x) && (x < knots[ii+1]));
    } else {
      ii2 = ii+1;
      p2 = p-1;
      my_bspline_internal(y1,x,ii,p2,knots,nknots);
      my_bspline_internal(y2,x,ii2,p2,knots,nknots);
      
      a = (x-knots[ii]) / (knots[ii+p]-knots[ii]);
      b = (knots[ii2+p]-x) / (knots[ii2+p]-knots[ii2]);
      output = a * y1 + b * y2;
    }
  }