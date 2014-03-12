#include "synlik.h"

using namespace Rcpp;

/*
  * Euler-multinomial sampling as in pomp package
  *
  * ARGS LIST:
  * size:      number of items extracted from the n = lenght(rates) boxes
  * rates:     rates of exit from each box (rates of decay of exponential function)
  * dt:        time during which the exits can occur
  */
inline void myReulermulti(double size, 
                          const  NumericVector & rates, 
                          const  double & dt,
                                 NumericVector & outSample) 
{
    
  int nBoxes = rates.size();
  int ii;
        
  // Comment from *here*  [a]
    if ( (size < 0) || (dt < 0.0) ) // || ( round(size + 0.5 - 1e-1) != size ) )  //SAFETY CHECK
    {
     Rcout << "size = " << size << ", dt = " << dt << 
     ", round(size + 0.5) = " << round(size + 0.5) << ", size = " << size << std::endl;
     stop("reulermultinom: either size < 0.0, dt < 0.0 or size in not int");
    } 
  // [a] to here for speed
  
  double totRate = std::accumulate(rates.begin(), rates.end(), 0.0);
  
  if (totRate > 0.0) 
  {
   size = R::rbinom(size, 1.0 - exp(-totRate*dt) ); // total number of events happening
  
  // Comment from *here*  [b]
  /* Often we size = NaN probably just because prob = 1.0 - exp(-totRate*dt) is almost 0
  so if that's the only reason, I don't issue a warning */
    if ( !R_FINITE(size) && ( 1.0 - exp(-totRate*dt) > 1e-4) )                      // SAFETY CHECK
{
  Rcout << "Warning euler-multi: size N = " << size << ", prob = " 
  << 1.0 - exp(-totRate*dt) << std::endl;
  size = 0.0;
}
// [b] to here for speed
  
  nBoxes -= 1;   // -1 because we can skip the last loop
  for (ii = 0; ii < nBoxes; ii++) {
    if (rates[ii] > totRate) totRate = rates[ii];
    outSample[ii] = ( (size > 0.0) && (totRate > 0.0) ) ? R::rbinom(size, rates[ii]/totRate) : 0.0;
    
    // Comment from *here*  [b]
    if ( !( R_FINITE(size)&&R_FINITE(totRate)&&                                //SAFETY CHECK
            R_FINITE(rates[ii])&&R_FINITE(outSample[ii]) ) )  
    Rcout << "Warning euler-multi: we are getting NaNs " << std::endl; 
    // [c] to here for speed
    
    size -= outSample[ii];
    totRate -= rates[ii];
  }
  outSample[nBoxes] = size;
} else {
  outSample = outSample * 0.0;
}
  }



/* 
 * Adds binomial noise to the observations
 */
  
  inline void binomObsNoise(const NumericMatrix & vectRho, 
                            const NumericMatrix & hidStates,
                                  NumericMatrix & obsStates)
{ 
    double rho; 
    
    int nSimul = hidStates.nrow();
    int nObs   = hidStates.ncol();
        
    for(int simulIndex = 0; simulIndex < nSimul; simulIndex++)
    {
      rho = exp(vectRho(simulIndex, 0));
      
      for(int obsIndex = 0; obsIndex < nObs; obsIndex++)
      {
        obsStates(simulIndex, obsIndex) = R::rbinom(hidStates(simulIndex, obsIndex), rho);
      }
    }
}
  
  
/* 
 * Adds Negative-Binomial noise to the observations
 */
  
inline void NegBinomObsNoise(const NumericMatrix & params, 
                             const NumericMatrix & hidStates,
                                   NumericMatrix & obsStates)
{ 
    double rho; 
    double size;
    
    int nSimul = hidStates.nrow();
    int nObs   = hidStates.ncol();
        
    for(int simulIndex = 0; simulIndex < nSimul; simulIndex++)
    {
      rho = exp(params(simulIndex, 0));
      size = 1 / exp(params(simulIndex, 1)); // size = 1 / (psiSquare^2)
      
      for(int obsIndex = 0; obsIndex < nObs; obsIndex++)
      {
        obsStates(simulIndex, obsIndex) = rnbinom_mu(size, rho * hidStates(simulIndex, obsIndex));
      }
    }
    
}


