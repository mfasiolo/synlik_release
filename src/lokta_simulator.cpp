
#include "synlik.h"

SEXP loktaCpp(SEXP nObs_, SEXP nsim_, SEXP params_, SEXP nBurn_, 
           SEXP totTime_, SEXP nSteps_, SEXP x0_, SEXP y0_, SEXP randInit_)
{  
  using namespace Rcpp;
  
  try{
    int nObs = as<int>(nObs_);
    int nsim = as<int>(nsim_);
    NumericMatrix params = as<NumericMatrix>(params_);
    int nBurn = as<int>(nBurn_);
    int totTime = as<int>(totTime_);
    int nSteps = as<int>(nSteps_);
    double x0 = as<double>(x0_);
    double y0 = as<double>(y0_);
    bool randInit = as<bool>(randInit_);
    
    RNGScope scope;
    
    int nParams = params.ncol(); 
    bool multiParams = false;
    
    if(nParams != 4) stop("Wrong number of parameters");
    if(params.nrow() > 1) { multiParams = true; }
    if(multiParams == true && params.nrow() != nsim) 
      stop("Number of parameters vectors is different from the number of simulations");
    
    double a = exp(params(0, 0));
    double b = exp(params(0, 1));
    double sigmaProc = exp(params(0, 2)) / nSteps;
    double sigmaObs = exp(params(0, 3));
    
    int totSteps = nSteps * (nObs + nBurn);
    int burnSteps = nSteps * nBurn;
    int obsSteps = nSteps * nObs;
    
    double h = static_cast<double>(totTime) / (nObs * nSteps);
    double dth = h / 2.0;
    
    NumericMatrix resultX(nsim, nObs);
    NumericMatrix resultY(nsim, nObs);
    
    NumericVector initX(nsim);
    NumericVector initY(nsim);
    
    if(randInit)
    {
      initX = runif(nsim, 0.75, 1.5);
      initY = initX - runif(nsim, -0.5, 0.0);
    } else {
      initX = initX + x0;
      initY = initY + y0;
    }
    
    NumericVector procNoise( rnorm( (totSteps+1)*nsim ) );
    NumericVector obsNoise( rnorm( (nObs+1) * nsim ) );
    
    NumericVector::iterator procNoiseIter = procNoise.begin();
    NumericVector::iterator obsNoiseIter = obsNoise.begin();
    NumericVector::iterator initXiter = initX.begin();
    NumericVector::iterator initYiter = initY.begin();
    
    int obsCounter = 0;
    
    double x, y;   
    double f1, f21, f2, f22, f23, f3, f4, f24;
    
    for(int iSimul = 0; iSimul < nsim; iSimul++, initXiter++, initYiter++)
    {
      if(multiParams)
      {
        a = exp(params(iSimul, 0));
        b = exp(params(iSimul, 1));
        sigmaProc = exp(params(iSimul, 2));
        sigmaObs = exp(params(iSimul, 3)) / nSteps;
      }
      
      x = *initXiter;
      y = *initYiter;
      
      for(int ii = 1; ii <= burnSteps; ii++, procNoiseIter++)
      {
        f1 = a*x - x*y;  
        f21 = b*x*y - y;
        
        f2 = a * ( x + f1 * dth ) - ( x + f1 * dth) * ( y + f21 * dth );
        f22 = b * ( x + f1 * dth )*( y + f21 * dth ) - ( y + f21 * dth );
        
        f3 = a * ( x + f2 * dth ) - ( x + f2 * dth ) * ( y + f22 * dth );
        f23 = b * ( x + f2 * dth )*( y + f22 * dth ) - ( y + f22 * dth );
        
        f4 = a * ( x + f3 * h ) - ( x + f3 * h ) * ( y + f23 * h );
        f24 = b * ( x + f3 * h )*( y + f23 * h ) - ( y + f23 * h );
        
        x = (x + (1.0/6.0) * ( f1 + 2*f2 + 2*f3 + f4 ) * h) * exp( *procNoiseIter * sigmaProc ) ;
        y = y + (1.0/6.0) * ( f21 + 2*f22 + 2*f23 + f24 ) * h;
        
        if(x < 0.0 || y < 0.0) 
        {
          stop("Going below zero");   
        }
      }
      
      obsCounter = 0;
      
      for(int ii = 0; ii < obsSteps; ii++, procNoiseIter++)
      {
        
        f1 = a*x - x*y;  
        f21 = b*x*y - y;
        
        f2 = a * ( x + f1 * dth ) - ( x + f1 * dth) * ( y + f21 * dth );
        f22 = b * ( x + f1 * dth )*( y + f21 * dth ) - ( y + f21 * dth );
        
        f3 = a * ( x + f2 * dth ) - ( x + f2 * dth ) * ( y + f22 * dth );
        f23 = b * ( x + f2 * dth )*( y + f22 * dth ) - ( y + f22 * dth );
        
        f4 = a * ( x + f3 * h ) - ( x + f3 * h ) * ( y + f23 * h );
        f24 = b * ( x + f3 * h )*( y + f23 * h ) - ( y + f23 * h );
        
        x = (x + (1.0/6.0) * ( f1 + 2.0*f2 + 2.0*f3 + f4 ) * h) * exp( *procNoiseIter * sigmaProc ) ;
        y = y + (1.0/6.0) * ( f21 + 2.0*f22 + 2.0*f23 + f24 ) * h;
        
        if(x < 0.0 || y < 0.0) 
        {
          stop("Going below zero");   
        }
        
        if( (ii % nSteps) == 0)
        {
          resultX(iSimul, obsCounter) = x * exp( *obsNoiseIter * sigmaObs );
          resultY(iSimul, obsCounter) = y;
          obsCounter++;
          obsNoiseIter++;
        }
        
      }
    }
    
    return List::create( _["prey"] = resultX, _["predator"] = resultY);
    
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
}
