
#include "synlik.h"
#include "malaria_header.h"

using namespace Rcpp;

/*
  * stepSIR = simulated one step of a seasonal SIR model with gamma noise
*           which corresponds to the last model described in "Introduction to pomp"
*
  * ARGS LIST
* params:      nSimul by nParams each row is the set of parameters used for one simulation
* splineParams: nSimul by nBasis its a matrix
* splineBasis: (nObs * nSteps) by nBasis its a matrix of spline basis covariates.
* initStates:  nSimul by 3 each row is the set of initial values 
               (susceptibles, infected and recovered) used for one simulation
* nStep:       finite difference step in which one obsInterval is divided (mininum nStep = 1)
* nSimul:      number of simulated paths
* obsInterval: the interval between 2 observation (which will be subdivided in nSteps)
* T0:          initial time, useful if there is a periodic component in the system
*
  * NOTES:  
  * - susceptibles, infected and recovered in general don't sum up exactly to N but their sum
*   should oscillate around N
*/

NumericMatrix stepSIR(const NumericMatrix & params, 
                      const NumericMatrix & splineParams,
                      const NumericMatrix & splineBasis,
                      const NumericMatrix & initStates, 
                      const int & nSteps, 
                      const int & nSimul, 
                      const double & obsInterval,
                      const int & currIteration)
{
 RNGScope scope;

 if(nSteps < 1) stop("nSteps cannot be < 1");
 
 /* Each simulated path should have it's own set of parameters */
  if(params.nrow() != nSimul)       stop("params.nrow() != nSimul");
  if(splineParams.nrow() != nSimul) stop("initStates.nrow() != nSimul");
  if(initStates.nrow() != nSimul)   stop("initStates.nrow() != nSimul");

const int nParams = params.ncol();
const double dt = obsInterval / nSteps;     // nSteps = 1 means no intermediate points
double Sus, Inf, Rec;                       // Susceptibles, Infected and Recovered
double CumNewRec = 0.0;                     // Cumulative number of recoveries in obsInterval

if(nParams != 5) stop("Wrong number of parameters");

/* unpack the parameters (for legibility only) */
  double N;           // population size
double gamma;       // recovery rate
double mu;          // birth rate = death rate
double beta;        // birth rate = death rate
double iota;        // latent force of infection
double betaSd;      // stdev of the gamma noise

double betaVar;
double birthRate; 

// Vector of spline coeffients and vector of spline covariates
NumericVector splineCoeffs(splineParams.ncol()); 
NumericVector splineCovar(splineBasis.ncol()); 

/* Output Susceptibles, Infected, Recovered and Cumulative new recoveries*/
  NumericMatrix output(nSimul, 4);

// Rates of exit from each class
NumericVector ratesSus(2);
NumericVector ratesInf(2);
NumericVector ratesRec(1);

// Number of transitions from each class
double births;
NumericVector perturbSus(2);
NumericVector perturbInf(2);
NumericVector perturbRec(1);

// Gamma noise vector 
NumericVector gammaNoise(nSteps);
NumericVector::iterator gammaIter;

/* Simulate all the nSimul paths */
  for(int simulIndex = 0; simulIndex < nSimul; simulIndex++)
  {
    // Setting parameters for the simulIndex-th path (beta determined later!)
    N       =  exp(params(simulIndex, 0));         
    gamma   =  exp(params(simulIndex, 1));     
    mu      =  exp(params(simulIndex, 2)); 
    iota    =  exp(params(simulIndex, 3));      
    betaSd  =  exp(params(simulIndex, 4)); 
    
    betaVar = betaSd * betaSd;
    birthRate = mu * N * dt;
    
    /* Setting specific spline coefficients, taken from the simulIndex-th 
    * row of splineParams and stored in splineCoeffs                      */
      RowToVector(splineParams, splineCoeffs, simulIndex);
    RowToVector(splineBasis, splineCovar, currIteration);
    
    // Calculating specific beta as inner prodoct of spline coefficient and spline basis (covariate) 
    beta  = exp( std::inner_product(splineCoeffs.begin(), splineCoeffs.end(), splineCovar.begin(), 0.0) );
    
    // Setting initial states  
    Sus = initStates(simulIndex, 0);
    Inf = initStates(simulIndex, 1);
    Rec = initStates(simulIndex, 2);
    
    // Setting specific rates of exit
    ratesSus[1] = mu;
    ratesInf[0] = gamma; 
    ratesInf[1] = mu;
    ratesRec[0] = mu;
    
    // Generating specific gamma noise and resetting cumulative recoveries to 0
    CumNewRec = 0.0;
    gammaNoise = rgamma(nSteps, dt/betaVar, betaVar);
    gammaIter = gammaNoise.begin();
    
    for(int stepIndex = 0; stepIndex < nSteps; stepIndex++, gammaIter++) 
    {
      ratesSus[0] =  ( beta * Inf / N + iota ) * *gammaIter / dt;
      
      births = R::rpois(birthRate);                   // New births
      myReulermulti(Sus, ratesSus, dt, perturbSus);   // exits from S in perturbSus
      myReulermulti(Inf, ratesInf, dt, perturbInf);   // exits from I in perturbInf
      myReulermulti(Rec, ratesRec, dt, perturbRec);   // exits from R in perturbRec
      
      Sus += ( births - perturbSus[0] - perturbSus[1] );
      Inf += ( perturbSus[0] - perturbInf[0] - perturbInf[1] );
      Rec += ( perturbInf[0] - perturbRec[0] );
      CumNewRec += perturbInf[0];
    }
    
    output(simulIndex, 0) = Sus;
    output(simulIndex, 1) = Inf;
    output(simulIndex, 2) = Rec;
    output(simulIndex, 3) = CumNewRec;
  }

return output;
}



/*
  * The SIR simulator: it simulates a general SIR model
* 
  * ARGS LIST
* procParams: parameters of the hidden process
* obsParams:  parameters of the observation process
* initStates:  nSimul by 2 (or 1 by 2) each row is the set of initial values (susceptibles and infected) used for one simulation
* nStep:       finite difference step in which one obsInterval is divided (mininum nStep = 1)
* nSimul:      number of simulated paths
* obsInterval: the interval between 2 observation (which will be subdivided in nSteps)
* T0:          initial time, useful if there is a periodic component in the system
*/
  
SEXP sirSimulCpp(SEXP procParams_, 
                 SEXP obsParams_,
                 SEXP splineParams_,
                 SEXP initStates_, 
                 SEXP nSteps_, 
                 SEXP nSimul_, 
                 SEXP obsInterval_,
                 SEXP nObs_,
                 SEXP T0_,
                 SEXP splineSettings_) 
{                 
  using namespace Rcpp;
  
  try{
  
  /* Setting default values for splineSettings
  *  WARNING:  nBasis and degree are INTs, while period is a DOUBLE*/
  NumericMatrix procParams = as<NumericMatrix>(procParams_);
  NumericMatrix obsParams = as<NumericMatrix>(obsParams_);
  NumericMatrix splineParams = as<NumericMatrix>(splineParams_);
  NumericMatrix initStates = as<NumericMatrix>(initStates_);
  int nSteps = as<int>(nSteps_); 
  int nSimul = as<int>(nSimul_); 
  double obsInterval = as<double>(obsInterval_);
  int nObs = as<int>(nObs_);
  double T0 = as<double>(T0_);
  NumericVector splineSettings = as<NumericVector>(splineSettings_);
  
  
  if(splineSettings.size() == 0) splineSettings = NumericVector::create(_["nBasis"] = 3.0, _["degree"] = 2.0, _["period"] = 1.0);
  
  if(nSteps < 1) stop("nSteps cannot be < 1");
  if(nObs < 1) stop("nObs cannot be < 1");
  
  if(procParams.nrow() != nSimul) stop("procParams.nrow() != nSimul");
  if(obsParams.nrow() != nSimul) stop("obsParams.nrow() != nSimul");
  if(splineParams.nrow() != nSimul) stop("initStates.nrow() != nSimul");
  if(initStates.nrow() != nSimul) stop("splineParams.nrow() != nSimul");

  if(initStates.ncol() != 3) stop("initStates should have 3 columns: Susc, Inf and Rec");
  if(splineParams.ncol() != as<int>(splineSettings["nBasis"])) stop("ncol(splineParams) != seasSettings[[\"nBasis\"]]");
  
  NumericMatrix susceptibles(nSimul, nObs);
  NumericMatrix infected    (nSimul, nObs);
  NumericMatrix recovered   (nSimul, nObs);
  NumericMatrix hiddenCases (nSimul, nObs);
  NumericMatrix obsCases    (nSimul, nObs);
  
  NumericMatrix tmpProcResults(nSimul, 4);
  
  // Setting up spline coefficients for seasonal component
  NumericVector tBasis = seqCpp(T0, T0 + nObs*obsInterval, obsInterval);
  NumericMatrix splineBasis( my_periodic_bspline_basis_creator(tBasis, as<int>(splineSettings["nBasis"]), 
                                                               as<int>(splineSettings["degree"]), 
                                                               as<double>(splineSettings["period"])) );
  int obsCounter = 0; 
  
  /* Simulating the hidden states */
    tmpProcResults = stepSIR(procParams, splineParams, splineBasis, initStates, nSteps, nSimul, obsInterval, obsCounter);
  
  susceptibles(_ , 0) = tmpProcResults(_ , 0);
  infected    (_ , 0) = tmpProcResults(_ , 1);
  recovered   (_ , 0) = tmpProcResults(_ , 2);
  hiddenCases (_ , 0) = tmpProcResults(_ , 3);
  
  for(obsCounter = 1; obsCounter < nObs; obsCounter++)
  {
    tmpProcResults = stepSIR(procParams, splineParams, splineBasis, tmpProcResults, nSteps, nSimul, obsInterval, obsCounter);
    
    susceptibles(_ , obsCounter) = tmpProcResults(_ , 0);
    infected    (_ , obsCounter) = tmpProcResults(_ , 1);
    recovered   (_ , obsCounter) = tmpProcResults(_ , 2);
    hiddenCases (_ , obsCounter) = tmpProcResults(_ , 3);
  }
  
  /* Simulating the observed cases */
    binomObsNoise(obsParams, hiddenCases, obsCases);
  
  return List::create(_["S"] = susceptibles, 
                      _["I"] = infected, 
                      _["R"] = recovered, 
                      _["H_Cases"] = hiddenCases, 
                      _["Obs_Cases"] = obsCases);
                      
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  
} 










/*******************************************************************
 * *****************************************************************
 * **************  Malaria model of Bhadhra 2011      ************** 
 * ****************************************************************
 * ****************************************************************
 */
 
/*
* malStepSIR = simulated one step of a seasonal SIR malaria model of Bhadhra 2011 
*
* ARGS LIST
* params:      nSimul by nParams each row is the set of parameters used for one simulation (log scale)
* splineParams: nSimul by nBasis its a matrix
* splineBasis: (nObs * nSteps) by nBasis its a matrix of spline basis covariates.
* initStates:  nSimul by 7 each row is the set of initial values 
(susceptibles, infected and recovered) used for one simulation
* Pop:         size of the population
* nStep:       finite difference step in which one obsInterval is divided (mininum nStep = 1)
* nSimul:      number of simulated paths
* obsInterval: the interval between 2 observation (which will be subdivided in nSteps)
* T0:          initial time, useful if there is a periodic component in the system
*
* NOTES:  
* - susceptibles, infected and recovered in general don't sum up exactly to N but their sum
*   should oscillate around N
*/

NumericMatrix malStepSIR(const NumericMatrix & params, 
                         const NumericMatrix & splineParams,
                         const NumericMatrix & splineBasis,
                         const NumericMatrix & initStates,
                         const NumericVector & rainfall,
                         const double Pop, 
                         const int nSteps, 
                         const int nSimul, 
                         const double obsInterval,
                         const int currIteration)
{
  RNGScope scope;
  
  static int counter = 0;
  
  if(nSteps < 1) stop("nSteps cannot be < 1");

  /* Each simulated path should have it's own set of parameters */
  if(params.nrow() != nSimul)       stop("params.nrow() != nSimul");
  if(splineParams.nrow() != nSimul) stop("initStates.nrow() != nSimul");
  if(initStates.nrow() != nSimul)   stop("initStates.nrow() != nSimul");
    
  const int nParams = params.ncol();
  const double dt = obsInterval / nSteps;     // nSteps = 1 means no intermediate points
  double S1, S2, E, I1, I2, LAM1, LAM2;       // Susceptibles1, Susceptibles2, Exposed, Inf1, Inf2, lambda1 and lambda2
  double CumNewInf = 0.0;                     // Cumulative number of infections (E to I1) in obsInterval
  
  if(nParams != 10) stop("Wrong number of parameters");
  
  /* 
   *   Unpack the parameters (for legibility only) 
   */
  // Inter-class transition rates (mu_S1E missing because it's time dependent)
  double mu_EI1;      // rate from E to I1
  double mu_I1I2;     // rate from I1 to I2
  double mu_I2S2;     // rate from I2 to S2
  double mu_S2S1;     // rate from I2 to S2 
  double mu_I1S1;     // rate from I1 to S1
  double mu_S2I2;     // rate from S2 to I2
   
  // Death rate from any class
  const double mu_XD = 1.0/50.0;   
  
  // Constant params
  const double beta_bar = 1;  // Fixed to the value of the article
  const double kappa = 2.0;
  
  // Other params
  double q;
  double covarCoeff;          // Coefficient of the rainfall covariate
  double tau;                 // Mean of the gamma-distributed latency
  double sigma;               // Intensity of gamma noise
  double c;                   // Force of re-infection

  // Other temporary variables
  double betaCovar;
  double sigmaSquare;         // = sigma * sigma
  
  // Vector of spline coeffients and vector of spline covariates
  NumericVector splineCoeffs(splineParams.ncol()); 
  NumericVector splineCovar (splineBasis.ncol()); 
  
  /* Output S1, S2, E, I1, I2, LAM1, LAM2 and Cumulative new recoveries */
  NumericMatrix output(nSimul, initStates.ncol() + 1);  // +1 for CumNewRec
  
  // Rates of exit from each class
  double birthRate;
  NumericVector ratesS1(2);
  NumericVector ratesS2(3);
  NumericVector ratesE (2);
  NumericVector ratesI1(3);
  NumericVector ratesI2(2);
  
  // Number of transitions from each class
  double births;
  NumericVector perturbS1(2);
  NumericVector perturbS2(3);
  NumericVector perturbE (2);
  NumericVector perturbI1(3);
  NumericVector perturbI2(2);
    
  // Gamma noise vector 
  NumericVector gammaNoise(nSteps);
  NumericVector::iterator gammaIter;
  
  // Normalization used to make the compartments sum up to the total population Pop
  double normalization;
  
  /* Simulate all the nSimul paths */
  for(int simulIndex = 0; simulIndex < nSimul; simulIndex++)
  {
    
  /* Setting parameters for the simulIndex-th path */
  // Transition rates (mu)        
  mu_EI1  = exp(params(simulIndex, 0));     
  mu_I1I2 = exp(params(simulIndex, 1));    
  mu_I2S2 = exp(params(simulIndex, 2));  
  mu_S2S1 = exp(params(simulIndex, 3));  
  mu_I1S1 = exp(params(simulIndex, 4));       
  
  // Other parameters
  q          = exp(params(simulIndex, 5));
  covarCoeff = exp(params(simulIndex, 6));                       
  sigma      = exp(params(simulIndex, 7)); 
  c          = exp(params(simulIndex, 8));
  tau        = exp(params(simulIndex, 9));
  
  // Derived parameters
  sigmaSquare = sigma * sigma;
  birthRate = mu_XD * Pop * dt;
  
  /* Setting specific spline coefficients, taken from the simulIndex-th 
   * row of splineParams and stored in splineCoeffs                      */
  RowToVector(splineParams, splineCoeffs, simulIndex);
  
  // Setting the initial states  
  S1 = initStates(simulIndex, 0);
  S2 = initStates(simulIndex, 1);
  E  = initStates(simulIndex, 2);
  I1 = initStates(simulIndex, 3);
  I2 = initStates(simulIndex, 4);
  LAM1 = initStates(simulIndex, 5);
  LAM2 = initStates(simulIndex, 6);
  
  /* Setting specific rates of exit 
   * ratesS1[0] = mu_S1E and ratesS2[1] = mu_S2I2 = c * mu_S1E are time-varying and will be specified later
   */
  ratesS1[1] = mu_XD;
  
  ratesS2[0] = mu_S2S1; ratesS2[2] = mu_XD;
  
  ratesE [0] = mu_EI1;  ratesE [1] = mu_XD;
  
  ratesI1[0] = mu_I1S1; ratesI1[1] = mu_I1I2; ratesI1[2] = mu_XD;
  
  ratesI2[0] = mu_I2S2; ratesI2[1] = mu_XD;
   
  // Generating specific gamma noise and resetting cumulative new infections to 0
  CumNewInf = 0.0;
  gammaNoise = rgamma(nSteps, dt/sigmaSquare, sigmaSquare);
  gammaIter = gammaNoise.begin();
  
  // Inner loop with nSteps of multinomial sampling: here total population (Pop) is held constant 
   for(int stepIndex = 0; stepIndex < nSteps; stepIndex++, gammaIter++) 
  {
  
  // Calculating specific effect of rain and seasonal component 
  // as inner product of spline coefficient and spline basis 
  RowToVector(splineBasis, splineCovar, currIteration*nSteps + stepIndex);
  betaCovar  = exp( rainfall[currIteration*nSteps + stepIndex] * covarCoeff + 
                    std::inner_product(splineCoeffs.begin(), splineCoeffs.end(), splineCovar.begin(), 0.0) );
  
  // Formulas (8) and (9) in the paper
  LAM1 = std::max( LAM1 + (beta_bar * (I1 + q * I2) / Pop  * betaCovar - LAM1 )* kappa / tau * *gammaIter, 0.0);
  LAM2 = std::max( LAM2 + (LAM1 - LAM2) * kappa / tau * dt, 0.0 );
  
  ratesS1[0] = LAM2;
  ratesS2[1] = c * ratesS1[0];     // mu_S2I2 = c * mu_S1E 
  
  // Simulate the perturbations
  births = R::rpois(birthRate);                // New births
  myReulermulti(S1, ratesS1, dt, perturbS1);   // exits from S1 in perturbS1
  myReulermulti(S2, ratesS2, dt, perturbS2);   // exits from S2 in perturbS2
  myReulermulti( E,  ratesE, dt, perturbE );   // exits from  E in perturbE
  myReulermulti(I1, ratesI1, dt, perturbI1);   // exits from I1 in perturbI1
  myReulermulti(I2, ratesI2, dt, perturbI2);   // exits from I2 in perturbI2
  
  // Perturbing the compartments
  S1 += ( births + perturbS2[0] + perturbI1[0] - perturbS1[0] - perturbS1[1] );
  S2 += ( perturbI2[0] - perturbS2[0] - perturbS2[1] - perturbS2[2] );
   E += ( perturbS1[0] - perturbE[0]  - perturbE[1] );
  I1 += ( perturbE[0]  - perturbI1[0] - perturbI1[1] - perturbI1[2]);
  I2 += ( perturbI1[1] + perturbS2[1] - perturbI2[0] - perturbI2[1]);
  CumNewInf += perturbE[0];
  }
  
  //Normalizing to have all the compartments sum up to the total population
  normalization =  Pop / (S1 + S2 + E + I1 + I2);
  
  output(simulIndex, 0) = round(S1*normalization);
  output(simulIndex, 1) = round(S2*normalization);
  output(simulIndex, 2) = round( E*normalization);
  output(simulIndex, 3) = round(I1*normalization);
  output(simulIndex, 4) = round(I2*normalization);
  output(simulIndex, 5) = LAM1;
  output(simulIndex, 6) = LAM2;
  output(simulIndex, 7) = round(CumNewInf*normalization);
}
  
  return output;
}

//[1] Here you could put *gammaIter/dt and the gamma(zz, simulIndex)*dt in the formula for ratesS1[0],
//    but you don't get the same output unless you change also the initial values for lambda.





/* Simulator for the S2EI2 malaria model of Bhadhra et at. 2011.
 *
 * ARGS:
 * procParams_: matrix of size (nSimul_ by 10), the i-th row is the set of process parameters 
 *                    used for the i-th simulation. The columns represent, in order, the _LOG_ of parameters
 *                    mu_EI1, mu_I1I2, mu_I2S2, mu_S2S1, mu_I1S1, q, covarCoeff, sigma, c and tau.  
 * obsParams_: matrix of size (nSimul_ by 2), the i-th row is the set of parameters of the observational 
 *                   process used for the i-th simulation. The columns represent, in order, the _LOG_ of parameters
 *                   rho and psiSquare.                  
 * splineParams_: matrix of size (nSimul by splineSettings[["nBasis"]]), the i-th row is the set of spline coefficients 
 *                      used for the i-th simulation.
 * initStates_:  matrix of size (nSimul by 7), the i-th row is the set of initial values 
 *                     used for the i-th simulation. The columns represent, in order, the states
 *                     S1, S2, E, I1, I2, lambda1 and lambda2.
 *                     NB: S1, S2, E, I1 and I2 have to be integers that sum up to population_[0]
 * population_: vector of length ( nObs_ ) representing the total population at each month. In this model I'm assuming
 *              that the population stays constant between months (which is roughly true).
 * rainfall_: vector of length ( nObs_ * nSteps_) representing the rainfall at each step.
 * nSteps_: finite difference step in which one month is divided (mininum nStep = 1, i.e. you directly step from
 *          one month to the next).
 * nSimul_: number of simulated paths
 * obsInterval: the length of the interval between 2 observation (which will be subdivided in nSteps). Typically
                it is going to be 1/12 so each observation represents a month.
 * nObs_: number of observations (months) of each simulated trajectory.
 * T0_: initial time of the simulations. It is used to set up the spline basis.
 * splineSettings_: named list where [["nBasis"]], [["degree"]] and [["period"]] are the (self-explanatory)
                    settings used to set up the spline basis for the seasonal component of the model.
 *
 * RETURN:
 * A named list where S1, S2, E, I1, I2, lambda1 and lambda2 are matrices of size ( nsim by nObs ) each row
 * representing a simulated trajectory of a particular state.
 */
 
SEXP malSirSimulCpp( SEXP procParams_, 
                     SEXP obsParams_,
                     SEXP splineParams_,
                     SEXP initStates_,
                     SEXP population_,
                     SEXP rainfall_,
                     SEXP nSteps_, 
                     SEXP nSimul_, 
                     SEXP obsInterval_,
                     SEXP nObs_,
                     SEXP T0_,
                     SEXP splineSettings_) 
{ 
  using namespace Rcpp;
  
  try{
  
  NumericMatrix procParams = as<NumericMatrix>(procParams_);
  NumericMatrix obsParams = as<NumericMatrix>(obsParams_);
  NumericMatrix splineParams = as<NumericMatrix>(splineParams_);
  NumericMatrix initStates = as<NumericMatrix>(initStates_);
  NumericVector population = as<NumericVector>(population_);
  NumericVector rainfall = as<NumericVector>(rainfall_);
  int nSteps = as<int>(nSteps_);
  int nSimul = as<int>(nSimul_); 
  double obsInterval = as<double>(obsInterval_);
  int nObs = as<int>(nObs_);
  double T0 = as<double>(T0_);
  List splineSettings = as<List>(splineSettings_);
  
  if(rainfall.size() != nObs * nSteps) stop("rainfall.size() != nObs * nSteps");
  if(population.size() != nObs) stop("population.size() != nObs");
  
  if(splineParams.ncol() != as<int>(splineSettings["nBasis"])) stop("splineParams.ncol() != as<int>(splineSettings[\"nBasis\"])");
  
  NumericMatrix S1(nSimul, nObs);
  NumericMatrix S2(nSimul, nObs);
  NumericMatrix E (nSimul, nObs);
  NumericMatrix I1(nSimul, nObs);
  NumericMatrix I2(nSimul, nObs);
  NumericMatrix LAM1(nSimul, nObs);
  NumericMatrix LAM2(nSimul, nObs);
  NumericMatrix hiddenCases(nSimul, nObs);
  NumericMatrix obsCases   (nSimul, nObs);
    
  // Temporary matrix that holds the new states at each step
  NumericMatrix tmpProcResults(nSimul, initStates.ncol() + 1);

  // Setting up spline coefficients for seasonal component
  NumericVector tBasis = seqCpp(T0, T0 + nObs*obsInterval, obsInterval/nSteps);
  NumericMatrix splineBasis( my_periodic_bspline_basis_creator(tBasis, as<int>(splineSettings["nBasis"]), 
                                                                       as<int>(splineSettings["degree"]), 
                                                                       as<double>(splineSettings["period"])) );
  int obsCounter = 0; 

  /* Simulating the hidden states */
  tmpProcResults = malStepSIR(procParams, splineParams, splineBasis, initStates, rainfall,
                              population[obsCounter], nSteps, nSimul, obsInterval, obsCounter);

  S1(_ , 0)             = tmpProcResults(_ , 0);
  S2(_ , 0)             = tmpProcResults(_ , 1);
  E (_ , 0)             = tmpProcResults(_ , 2);
  I1 (_ , 0)            = tmpProcResults(_ , 3);
  I2 (_ , 0)            = tmpProcResults(_ , 4);
  LAM1 (_ , 0)          = tmpProcResults(_ , 5);
  LAM2 (_ , 0)          = tmpProcResults(_ , 6);
  hiddenCases (_ , 0)   = tmpProcResults(_ , 7);

  for(obsCounter = 1; obsCounter < nObs; obsCounter++)
  {
  tmpProcResults = malStepSIR(procParams, splineParams, splineBasis, tmpProcResults, rainfall,
                              population[obsCounter], nSteps, nSimul, obsInterval, obsCounter);
  S1(_ , obsCounter)                     = tmpProcResults(_ , 0);
  S2(_ , obsCounter)                     = tmpProcResults(_ , 1);
  E (_ , obsCounter)                     = tmpProcResults(_ , 2);
  I1 (_ , obsCounter)                    = tmpProcResults(_ , 3);
  I2 (_ , obsCounter)                    = tmpProcResults(_ , 4);
  LAM1 (_ , obsCounter)                  = tmpProcResults(_ , 5);
  LAM2 (_ , obsCounter)                  = tmpProcResults(_ , 6);
  hiddenCases (_ , obsCounter)           = tmpProcResults(_ , 7);
  }

  /* Simulating the observed cases using negative binomial*/
  NegBinomObsNoise( obsParams, hiddenCases, obsCases); 
  
  return List::create(_["S1"] = S1, 
                      _["S2"] = S2, 
                      _["E"]  = E,
                      _["I1"] = I1,
                      _["I2"] = I2, 
                      _["H_Cases"] = hiddenCases, 
                      _["Obs_Cases"] = obsCases,
                      _["LAM1"] = LAM1,
                      _["LAM2"] = LAM2
                      );
                      
 } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
   } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
}

} 






































/*  
             P O M P
*/

/*
* malStepSIR = simulated one step of a seasonal SIR malaria model of Bhadhra 2011 
*
* ARGS LIST
* params:      nSimul by nParams each row is the set of parameters used for one simulation
* splineParams: nSimul by nBasis its a matrix
* splineBasis: (nObs * nSteps) by nBasis its a matrix of spline basis covariates.
* initStates:  nSimul by 3 each row is the set of initial values 
(susceptibles, infected and recovered) used for one simulation
* Pop:         size of the population
* nStep:       finite difference step in which one obsInterval is divided (mininum nStep = 1)
* nSimul:      number of simulated paths
* obsInterval: the interval between 2 observation (which will be subdivided in nSteps)
* T0:          initial time, useful if there is a periodic component in the system
*
* NOTES:  
* - susceptibles, infected and recovered in general don't sum up exactly to N but their sum
*   should oscillate around N
*/

SEXP malStepPompCpp(SEXP params_, 
                    SEXP splineParams_,
                    SEXP splineBasis_,
                    SEXP initStates_,
                    SEXP rainfall_,
                    SEXP Pop_, 
                    SEXP nSteps_, 
                    SEXP nSimul_, 
                    SEXP obsInterval_,
                    SEXP currIteration_)
{
  using namespace Rcpp;
  
  try{
  RNGScope scope;
  
  const NumericMatrix params = as<NumericMatrix>(params_);
  const NumericMatrix splineParams = as<NumericMatrix>(splineParams_);
  const NumericMatrix splineBasis = as<NumericMatrix>(splineBasis_);
  const NumericMatrix initStates = as<NumericMatrix>(initStates_);
  const NumericVector rainfall = as<NumericVector>(rainfall_);
  const double Pop = as<double>(Pop_);
  const int nSteps = as<int>(nSteps_); 
  const int nSimul = as<int>(nSimul_); 
  const double obsInterval = as<double>(obsInterval_);
  const int currIteration = as<int>(currIteration_);
    
  if(nSteps < 1) stop("nSteps cannot be < 1");
  if(currIteration < 0) stop("currIteration < 0");

  /* Each simulated path should have it's own set of parameters */
  if(params.nrow() != nSimul)       stop("params.nrow() != nSimul");
  if(splineParams.nrow() != nSimul) stop("splineParams.nrow() != nSimul");
  if(initStates.nrow() != nSimul)   stop("initStates.nrow() != nSimul");
    
  const int nParams = params.ncol();
  const double dt = obsInterval / nSteps;     // nSteps = 1 means no intermediate points
  double S1, S2, E, I1, I2, LAM1, LAM2;       // Susceptibles1, Susceptibles2, Exposed, Inf1, Inf2, lambda1 and lambda2
  double CumNewInf = 0.0;                     // Cumulative number of infections (E to I1) in obsInterval
  
  if(nParams != 10) stop("Wrong number of parameters");
  
  // Inter-class transition rates (mu_S1E missing because it's time dependent)
  double mu_EI1;      // rate from E to I1
  double mu_I1I2;     // rate from I1 to I2
  double mu_I2S2;     // rate from I2 to S2
  double mu_S2S1;     // rate from I2 to S2 
  double mu_I1S1;     // rate from I1 to S1
  double mu_S2I2;     // rate from S2 to I2
   
  // Death rate from any class
  const double mu_XD = 1.0/50.0;   
  
  // Constant params
  const double beta_bar = 1;  // Fixed to the value of the article
  const double kappa = 2.0;
   
  // Other params
  double q;
  double covarCoeff;          // Coefficient of the rainfall covariate
  double tau;                 // Mean of the gamma-distributed latency
  double sigma;               // Intensity of gamma noise
  double c;                   // Force of re-infection

  // Other temporary variables
  double betaCovar;
  double sigmaSquare;         // = sigma * sigma
  
  // Vector of spline coeffients and vector of spline covariates
  NumericVector splineCoeffs(splineParams.ncol()); 
  NumericVector splineCovar (splineBasis.ncol()); 
  
  /* Output S1, S2, E, I1, I2, LAM1, LAM2 and Cumulative new recoveries */
  NumericMatrix output(nSimul, initStates.ncol()); 
  
  // Rates of exit from each class
  double birthRate;
  NumericVector ratesS1(2);
  NumericVector ratesS2(3);
  NumericVector ratesE (2);
  NumericVector ratesI1(3);
  NumericVector ratesI2(2);
  
  // Number of transitions from each class
  double births;
  NumericVector perturbS1(2);
  NumericVector perturbS2(3);
  NumericVector perturbE (2);
  NumericVector perturbI1(3);
  NumericVector perturbI2(2);
    
  // Gamma noise vector 
  NumericVector gammaNoise(nSteps);
  NumericVector::iterator gammaIter;
  
  // Normalization used to make the compartments sum up to the total population Pop
  double normalization;
  
  /* Simulate all the nSimul paths */
  for(int simulIndex = 0; simulIndex < nSimul; simulIndex++)
  {
    
  /* Setting parameters for the simulIndex-th path */
  // Transition rates (mu)        
  mu_EI1  = exp(params(simulIndex, 0));     
  mu_I1I2 = exp(params(simulIndex, 1));    
  mu_I2S2 = exp(params(simulIndex, 2));  
  mu_S2S1 = exp(params(simulIndex, 3));  
  mu_I1S1 = exp(params(simulIndex, 4));       
  
  // Other parameters
  q          = exp(params(simulIndex, 5));
  covarCoeff = exp(params(simulIndex, 6));                       
  sigma      = exp(params(simulIndex, 7)); 
  c          = exp(params(simulIndex, 8));
  tau        = exp(params(simulIndex, 9));
  
  // Derived parameters
  sigmaSquare = sigma * sigma;
  birthRate = mu_XD * Pop * dt;
  
  /* Setting specific spline coefficients, taken from the simulIndex-th 
   * row of splineParams and stored in splineCoeffs                      */
  RowToVector(splineParams, splineCoeffs, simulIndex);
  
  // Setting the initial states  
  S1 = initStates(simulIndex, 0);
  S2 = initStates(simulIndex, 1);
  E  = initStates(simulIndex, 2);
  I1 = initStates(simulIndex, 3);
  I2 = initStates(simulIndex, 4);
  LAM1 = initStates(simulIndex, 5);
  LAM2 = initStates(simulIndex, 6);
  
  /* Setting specific rates of exit 
   * ratesS1[0] = mu_S1E and ratesS2[1] = mu_S2I2 = c * mu_S1E are time-varying and will be specified later
   */
  ratesS1[1] = mu_XD;
  
  ratesS2[0] = mu_S2S1; ratesS2[2] = mu_XD;
  
  ratesE [0] = mu_EI1;  ratesE [1] = mu_XD;
  
  ratesI1[0] = mu_I1S1; ratesI1[1] = mu_I1I2; ratesI1[2] = mu_XD;
  
  ratesI2[0] = mu_I2S2; ratesI2[1] = mu_XD;
   
  // Generating specific gamma noise and resetting cumulative new infections to 0
  CumNewInf = 0.0;
  gammaNoise = rgamma(nSteps, dt/sigmaSquare, sigmaSquare);
  gammaIter = gammaNoise.begin();
  
  // Inner loop with nSteps of multinomial sampling: here total population (Pop) is held constant 
   for(int stepIndex = 0; stepIndex < nSteps; stepIndex++, gammaIter++) 
  {
  
  // Calculating specific effect of rain and seasonal component 
  // as inner product of spline coefficient and spline basis 
  RowToVector(splineBasis, splineCovar, currIteration*nSteps + stepIndex);
  betaCovar  = exp( rainfall[currIteration*nSteps + stepIndex] * covarCoeff + 
                    std::inner_product(splineCoeffs.begin(), splineCoeffs.end(), splineCovar.begin(), 0.0) );
  
  // Formulas (8) and (9) in the paper
  LAM1 = std::max( LAM1 + (beta_bar * (I1 + q * I2) / Pop  * betaCovar - LAM1 )* kappa / tau * *gammaIter, 0.0);
  LAM2 = std::max( LAM2 + (LAM1 - LAM2) * kappa / tau * dt, 0.0 );
  
  ratesS1[0] = LAM2;
  ratesS2[1] = c * ratesS1[0];     // mu_S2I2 = c * mu_S1E 
  
  // Simulate the perturbations
  births = R::rpois(birthRate);                // New births
  myReulermulti(S1, ratesS1, dt, perturbS1);   // exits from S1 in perturbS1
  myReulermulti(S2, ratesS2, dt, perturbS2);   // exits from S2 in perturbS2
  myReulermulti( E,  ratesE, dt, perturbE );   // exits from  E in perturbE
  myReulermulti(I1, ratesI1, dt, perturbI1);   // exits from I1 in perturbI1
  myReulermulti(I2, ratesI2, dt, perturbI2);   // exits from I2 in perturbI2
  
  // Perturbing the compartments
  S1 += ( births + perturbS2[0] + perturbI1[0] - perturbS1[0] - perturbS1[1] );
  S2 += ( perturbI2[0] - perturbS2[0] - perturbS2[1] - perturbS2[2] );
   E += ( perturbS1[0] - perturbE[0]  - perturbE[1] );
  I1 += ( perturbE[0]  - perturbI1[0] - perturbI1[1] - perturbI1[2]);
  I2 += ( perturbI1[1] + perturbS2[1] - perturbI2[0] - perturbI2[1]);
  CumNewInf += perturbE[0];
  }
  
  //Normalizing to have all the compartments sum up to the total population
  normalization =  Pop / (S1 + S2 + E + I1 + I2);
  
  output(simulIndex, 0) = round(S1*normalization);
  output(simulIndex, 1) = round(S2*normalization);
  output(simulIndex, 2) = round( E*normalization);
  output(simulIndex, 3) = round(I1*normalization);
  output(simulIndex, 4) = round(I2*normalization);
  output(simulIndex, 5) = LAM1;
  output(simulIndex, 6) = LAM2;
  output(simulIndex, 7) = round(CumNewInf*normalization);
}
  
  return wrap(output);
  
 } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
   } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
}

} 


//[1] Here you could put *gammaIter/dt and the gamma(zz, simulIndex)*dt in the formula for ratesS1[0],
//    but you don't get the same output unless you change also the initial values for lambda.

