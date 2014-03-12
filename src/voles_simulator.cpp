
#include "synlik.h"

/* 
 * Simulator for the full Voles model (Turchin and Ellner (2000): 
 * LIVING ON THE EDGE OF CHAOS: POPULATION DYNAMICS OF FENNOSCANDIAN VOLES)
 *
 * Using Runga-Kutta integrator of order 2.
 *
 * ARGS:
 *
 * nMon = (int) number of months which is the length of each trajectory .   
 * nSimul = (int) number of simulated trajectories.
 * nBurn = (int) number of burn-in semesters
 * params = (NumericMatrix) matrix either 1 x 9 or nSimul by 9 of parameters
 * nSteps = (int) the number of steps that divides one month, mininum is 1 (you step direcly from one observation to the next).
 * T0 = (double) this is the starting day. If T0 == 0 than the first observed month will correspond to January. 
 * randInit = (bool) if true the initial values of voles and weasels are simulated
 * startVole, startWeas = (NumericVectors) vectors of length nSimul constaining the initial values of log(voles) and log(weasels).
                          Needed only if randInit == true.
 * addObsNoise = (bool) if true the observational noise is added, otherwise no.
 *
 * NOTE - we have noise only on the voles process, not on the weasels (which will be affected indirectly).
 */
SEXP volesFullCpp(SEXP nMon_, 
                  SEXP nSimul_, 
                  SEXP nBurn_, 
                  SEXP params_, 
                  SEXP nSteps_, 
                  SEXP T0_, 
                  SEXP randInit_, 
                  SEXP startVole_, 
                  SEXP startWeas_, 
                  SEXP addObsNoise_)
{
  using namespace Rcpp;
  
  try{
    
    RNGScope scope;
    
    const int nMon = as<int>(nMon_);
    const int nSimul = as<int>(nSimul_);
    const int nBurn = as<int>(nBurn_);
    const int nSteps = as<int>(nSteps_);
    const NumericMatrix params = as<NumericMatrix>(params_); 
    const double T0 = as<double>(T0_); 
    const bool randInit = as<bool>(randInit_); 
    const NumericVector startVole = as<NumericVector>(startVole_);
    const NumericVector startWeas = as<NumericVector>(startWeas_);
    const bool addObsNoise = as<bool>(addObsNoise_);
    
    const int totMon = nMon + nBurn;
    const double dt = (1.0 / 12.0) / nSteps;
    double currVole, currWeas, oldVole, oldWeas, expOldVole, expOldWeas;
    
    int nParams = params.ncol(); 
    bool multiParams = false;
    
    if(nParams != 9) stop("Wrong number of parameters");
    if(params.nrow() > 1) { multiParams = true; }
    if(multiParams == true && params.nrow() != nSimul) 
      stop("Number of parameters vectors is different from the number of simulations");
    
    double r = exp(params(0, 0));
    double e = exp(params(0, 1));
    double g = exp(params(0, 2));
    double h = exp(params(0, 3));
    double a = exp(params(0, 4));
    double d = exp(params(0, 5));
    double s = exp(params(0, 6));
    double sigmaProc = exp(params(0, 7));
    double phi = exp(params(0, 8));
    
    double Hsquared = pow(h, 2.0);
    
    NumericVector procNoise( rnorm(totMon*nSteps*nSimul, 0, sigmaProc) );     
    
    NumericMatrix Vole(nSimul, nMon);
    NumericMatrix Weas(nSimul, nMon);
    
    NumericVector initVole(nSimul);
    NumericVector initWeas(nSimul);
    
    if(randInit)
    {
      initVole = log( runif(nSimul) / 2.0 + 0.5 );                            
      initWeas = log( runif(nSimul) / 4.0 + 0.01 );                           
    }
    else {
      if( startVole.size() != nSimul || startWeas.size() != nSimul) stop("If randInit == false startVole and startWeas should be vector of length nSimul");
      initVole = clone(startVole);
      initWeas = clone(startWeas);
    }
    
    // Calculating the seasonal component in order not to calculate it for every simulated path.
    // sineValues in 12 * nSteps long so it misses the final day, but that is equal to the beginning.
    // There is an offset of 1/12 so that the seasonal component peaks when currTime == 8 / 12 (August) and 
    // has mininum when currTime == 2 / 12 (February)
  std::vector<double> sineValues(nSteps * 12, 0.0);          
  std::vector<double>::iterator theBegin = sineValues.begin();
  std::vector<double>::iterator theEnd = sineValues.end();
  std::vector<double>::iterator sineIter = theBegin;
  
  for(double currTime = T0; sineIter != theEnd; sineIter++, currTime += dt)
  {
    *sineIter = ( 1.0 - e * sin( 2 * PI * (currTime + 1.0 / 12.0) ) );  
  }
  
  // Iterators over process and observational noise, used for efficiency.
  NumericVector::iterator procNoiseIter = procNoise.begin();
  
  int dayIndex = 0;                   // Index of which observation we have to store
  int observationClock = nSteps;      // We store results every observationClock simulations (corresponding to 6 months)
  
  // These are used to stop the two inner loops (computed here for efficiency)
  int totalBurn =  nBurn * nSteps;    
  int totalSimul = nMon * nSteps;      
  
  for(int ii = 0; ii < nSimul; ii++)
  {
    oldVole = currVole = initVole[ii];
    oldWeas = currWeas = initWeas[ii];
    expOldVole = exp(initVole[ii]);
    expOldWeas = exp(initWeas[ii]);
    
    if(multiParams)
    {
     r = exp(params(ii, 0));
     e = exp(params(ii, 1));
     g = exp(params(ii, 2));
     h = exp(params(ii, 3));
     a = exp(params(ii, 4));
     d = exp(params(ii, 5));
     s = exp(params(ii, 6));
     sigmaProc = exp(params(ii, 7));
     phi = exp(params(ii, 8));
     Hsquared = pow(h, 2.0);
    }
    
    sineIter = theBegin;
    
    // Burn-in simulations start
    for(int jj = 0; jj < totalBurn; jj++, sineIter++, procNoiseIter++)
    {
      if(sineIter == theEnd){ sineIter = theBegin; }
      
      // First half step forward of Runga Kutta: calculating the intermediate point. No process noise needed here.
      expOldVole = exp(oldVole + (1.0 / expOldVole) * ( r * *sineIter * expOldVole - r * expOldVole * expOldVole - 
                   (g * expOldVole * expOldVole )/( expOldVole * expOldVole + Hsquared) -
                   (a * expOldWeas * expOldVole)/(expOldVole + d) ) * dt / 2.0);
                   
      expOldWeas = exp(oldWeas + (1.0 / expOldWeas) * 
                   (  s * *sineIter * expOldWeas  - (s * expOldWeas * expOldWeas)/expOldVole  ) * dt / 2.0);
      
      // Final step step forward using the derivative at the intermediate point. 
      // After calculating the step, process noise is added to the currVole.
      currVole = oldVole + ((1.0 / expOldVole) * ( r * *sineIter * expOldVole - r * expOldVole * expOldVole - 
                   (g * expOldVole * expOldVole )/( expOldVole * expOldVole + Hsquared) -
                   (a * expOldWeas * expOldVole)/(expOldVole + d) ) + *procNoiseIter) * dt; 
                   
      currWeas = oldWeas + (1.0 / expOldWeas) * (  s * *sineIter * expOldWeas  - 
                 (s * expOldWeas * expOldWeas)/expOldVole) * dt;    
                 
      oldVole = currVole;
      oldWeas = currWeas;
      expOldVole = exp(currVole);
      expOldWeas = exp(currWeas);
    }
    
    
    // Real simulations start
    dayIndex = 0;
    for(int jj = 1; jj <= totalSimul; jj++, sineIter++, procNoiseIter++) 
    {
      if(sineIter == theEnd){ sineIter = theBegin; }
      
      // First half step forward of Runga Kutta: calculating the intermediate point.
      expOldVole = exp(oldVole + (1.0 / expOldVole) * ( r * *sineIter * expOldVole - r * expOldVole * expOldVole - 
                   (g * expOldVole * expOldVole )/( expOldVole * expOldVole + Hsquared) -
                   (a * expOldWeas * expOldVole)/(expOldVole + d) ) * dt / 2.0);
      expOldWeas = exp(oldWeas + (1.0 / expOldWeas) * 
                   (  s * *sineIter * expOldWeas  - (s * expOldWeas * expOldWeas)/expOldVole  ) * dt / 2.0);
            
      // Final step step forward using the derivative at the intermediate point.      
      currVole = oldVole + ((1.0 / expOldVole) * ( r * *sineIter * expOldVole - r * expOldVole * expOldVole - 
                   (g * expOldVole * expOldVole )/( expOldVole * expOldVole + Hsquared) -
                   (a * expOldWeas * expOldVole)/(expOldVole + d) ) + *procNoiseIter) * dt; 
                   
      currWeas = oldWeas + (1.0 / expOldWeas) * (  s * *sineIter * expOldWeas  - 
                 (s * expOldWeas * expOldWeas)/expOldVole) * dt;  
                 
      oldVole = currVole;
      oldWeas = currWeas;
      expOldVole = exp(currVole);
      expOldWeas = exp(currWeas);
      
      // Store an observation every 6 months.
       if( ( jj % observationClock ) == 0) 
         {   
         Vole(ii, dayIndex) = addObsNoise? R::rpois(phi * exp(currVole)) : ( exp(currVole) );
         Weas(ii, dayIndex) = exp(currWeas);
         dayIndex++;
         }
        
    }
    
      
  }
  
  return List::create( _["voles"] = Vole, _["weasels"] = Weas);
  
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
}



/* 
 * Simulator for the standard Voles model (Turchin and Ellner (2000): page 3108)
 *
 * Using Runga-Kutta integrator of order 2.
 *
 * ARGS: same as for the full model apart from params = (NumericMatrix) matrix either 1 x 7 or nSimul by 7 of parameters
 *
 */
SEXP volesStdCpp(SEXP nMon_, SEXP nSimul_, SEXP nBurn_, SEXP params_, SEXP nSteps_, SEXP T0_, 
                  SEXP randInit_, SEXP startVole_, SEXP startWeas_, SEXP addObsNoise_)
{
  using namespace Rcpp;
  
  try{
    
    RNGScope scope;
    
    int nMon = as<int>(nMon_);
    int nSimul = as<int>(nSimul_);
    int nBurn = as<int>(nBurn_);
    int nSteps = as<int>(nSteps_);
    NumericMatrix params = as<NumericMatrix>(params_); 
    double T0 = as<double>(T0_); 
    bool randInit = as<bool>(randInit_); 
    NumericVector startVole = as<NumericVector>(startVole_);
    NumericVector startWeas = as<NumericVector>(startWeas_);
    bool addObsNoise = as<bool>(addObsNoise_);
    
    const int totMon = nMon + nBurn;
    const double dt = (1.0 / 12.0) / nSteps;
    double currVole, currWeas, oldVole, oldWeas, expOldVole, expOldWeas;
    
    int nParams = params.ncol(); 
    bool multiParams = false;
    
    if(nParams != 7) stop("Wrong number of parameters");
    if(params.nrow() > 1) { multiParams = true; }
    if(multiParams == true && params.nrow() != nSimul) 
      stop("Number of parameters vectors is different from the number of simulations");
    
    double r = exp(params(0, 0));
    double e = exp(params(0, 1));
    double beta = exp(params(0, 2));
    double d = exp(params(0, 3));
    double alpha = exp(params(0, 4));
    double sigmaProc = exp(params(0, 5));
    double phi = exp(params(0, 6));
    
    NumericVector procNoise( rnorm(totMon*nSteps*nSimul, 0, sigmaProc) );     
    
    NumericMatrix Vole(nSimul, nMon);
    NumericMatrix Weas(nSimul, nMon);
    
    NumericVector initVole(nSimul);
    NumericVector initWeas(nSimul);
    
    if(randInit)
    {
      initVole = log( runif(nSimul) / 2.0 + 0.5 );                            
      initWeas = log( runif(nSimul) / 4.0 + 0.01 );                           
    }
    else {
      if( startVole.size() != nSimul || startWeas.size() != nSimul) stop("If randInit == false startVole and startWeas should be vector of length nSimul");
      initVole = clone(startVole);
      initWeas = clone(startWeas);
    }
    
    // Calculating the seasonal component in order not to calculate it for every simulated path.
    // sineValues in 12 * nSteps long so it misses the final day, but that is equal to the beginning.
    // There is an offset of 1/12 so that the seasonal component peaks when currTime == 8 / 12 (August) and 
    // has mininum when currTime == 2 / 12 (February)
    std::vector<double> sineValues(nSteps * 12, 0.0);          
    std::vector<double>::iterator theBegin = sineValues.begin();
    std::vector<double>::iterator theEnd = sineValues.end();
    std::vector<double>::iterator sineIter = theBegin;
    
    for(double currTime = T0; sineIter != theEnd; sineIter++, currTime += dt)
    {
      *sineIter = ( 1.0 - e * sin( 2 * PI * (currTime + 1.0 / 12.0) ));
    }
    
    // Iterators over process and observational noise, used for efficiency.
    NumericVector::iterator procNoiseIter = procNoise.begin();
    
    int dayIndex = 0;                   // Index of which observation we have to store
    int observationClock = nSteps;      // We store results every observationClock simulations (corresponding to 6 months)
    
    // These are used to stop the two inner loops (computed here for efficiency)
    int totalBurn =  nBurn * nSteps;    
    int totalSimul = nMon * nSteps;       
    
    for(int ii = 0; ii < nSimul; ii++)
    {
      oldVole = currVole = initVole[ii];
      oldWeas = currWeas = initWeas[ii];
      expOldVole = exp(initVole[ii]);
      expOldWeas = exp(initWeas[ii]);
      
      if(multiParams)
      {
        r = exp(params(ii, 0));
        e = exp(params(ii, 1));
        beta = exp(params(ii, 2));
        d = exp(params(ii, 3));
        alpha = exp(params(ii, 4));
        sigmaProc = exp(params(ii, 5));
        phi = exp(params(ii, 6));
      }
      
      sineIter = theBegin;
      
      // Burn-in simulations start
      for(int jj = 0; jj < totalBurn; jj++, sineIter++, procNoiseIter++)
      {
        if(sineIter == theEnd){ sineIter = theBegin; }
        
        // First half step forward of Runga Kutta: calculating the intermediate point.
        expOldVole = exp(oldVole + (1.0 / expOldVole) * ( r * *sineIter * expOldVole - r * expOldVole * expOldVole - 
                                                            (expOldVole * expOldWeas)/(expOldVole + d) ) * dt / 2.0);
        expOldWeas = exp(oldWeas + (1.0 / expOldWeas) * (  (alpha * expOldWeas * expOldVole) / (expOldVole + d) - 
                                                             beta * (2.0 - *sineIter) * expOldWeas ) * dt / 2.0);
        
        // Final step step forward using the derivative at the intermediate point.  
        currVole = oldVole + ((1.0 / expOldVole) * ( r * *sineIter * expOldVole - r * expOldVole * expOldVole - 
                                                       (expOldVole * expOldWeas)/(expOldVole + d) ) + *procNoiseIter) * dt;  
        currWeas = oldWeas + (1.0 / expOldWeas) * (  (alpha * expOldWeas * expOldVole) / (expOldVole + d) - 
                                                       beta * (2.0 - *sineIter) * expOldWeas ) * dt;
        
        oldVole = currVole;
        oldWeas = currWeas;
        expOldVole = exp(currVole);
        expOldWeas = exp(currWeas);
      }
      
      // Real simulations start
      dayIndex = 0;
      for(int jj = 1; jj <= totalSimul; jj++, sineIter++, procNoiseIter++) 
      {
        if(sineIter == theEnd){ sineIter = theBegin; }
        
        // First half step forward of Runga Kutta: calculating the intermediate point.
        expOldVole = exp(oldVole + (1.0 / expOldVole) * ( r * *sineIter * expOldVole - r * expOldVole * expOldVole - 
                                                            (expOldVole * expOldWeas)/(expOldVole + d) ) * dt / 2.0);
        expOldWeas = exp(oldWeas + (1.0 / expOldWeas) * (  (alpha * expOldWeas * expOldVole) / (expOldVole + d) - 
                                                             beta * (2.0 - *sineIter) * expOldWeas ) * dt / 2.0);
        
        // Final step step forward using the derivative at the intermediate point.  
        currVole = oldVole + ((1.0 / expOldVole) * ( r * *sineIter * expOldVole - r * expOldVole * expOldVole - 
                                                       (expOldVole * expOldWeas)/(expOldVole + d) ) + *procNoiseIter) * dt;  
        currWeas = oldWeas + (1.0 / expOldWeas) * (  (alpha * expOldWeas * expOldVole) / (expOldVole + d) - 
                                                       beta * (2.0 - *sineIter) * expOldWeas ) * dt;
        
        oldVole = currVole;
        oldWeas = currWeas;
        expOldVole = exp(currVole);
        expOldWeas = exp(currWeas);
        
        // Store an observation every 6 months.
       if( ( jj % observationClock ) == 0) 
         {   
         Vole(ii, dayIndex) = addObsNoise? R::rpois( phi * exp(currVole) ) : ( exp(currVole) );
         Weas(ii, dayIndex) = exp(currWeas);
         dayIndex++;
         }
        
      }
      
    }
    
    return List::create( _["voles"] = Vole, _["weasels"] = Weas);
    
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
}



