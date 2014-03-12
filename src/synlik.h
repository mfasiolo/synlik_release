#ifndef _SYNLIK_H
#define _SYNLIK_H

#include <RcppArmadillo.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
 
RcppExport SEXP cleanStats(SEXP inMat);

RcppExport SEXP checkBoundsCpp(SEXP theMean_, SEXP cholFact_, SEXP toCheck_, SEXP upper_, SEXP lower_, SEXP output_);

RcppExport SEXP mahaCpp(SEXP X, SEXP mu, SEXP sigma, SEXP isChol);

RcppExport SEXP dmvnCpp(SEXP X_, SEXP mu_, SEXP sigma_, SEXP log_, SEXP isChol_);

/*
 * Simulators
 */

RcppExport SEXP simpleModelsWrap(SEXP model, SEXP days, SEXP nSimul, SEXP params, SEXP nBurn, SEXP randInit, SEXP initVal);

RcppExport SEXP loktaCpp(SEXP nObs_, SEXP nsim_, SEXP params_, SEXP nBurn_, SEXP totTime_, SEXP nSteps_, SEXP x0_, SEXP y0_, SEXP randInit_);

RcppExport SEXP volesFullCpp(SEXP nMon_, SEXP nSimul_, SEXP nBurn_, SEXP params_, SEXP nSteps_, SEXP T0_, 
                             SEXP randInit_, SEXP startVole_, SEXP startWeas_, SEXP addObsNoise_);

RcppExport SEXP volesStdCpp(SEXP nMon_, SEXP nSimul_, SEXP nBurn_, SEXP params_, SEXP nSteps_, SEXP T0_, 
                            SEXP randInit_, SEXP startVole_, SEXP startWeas_, SEXP addObsNoise_);
                            
RcppExport SEXP sirSimulCpp(SEXP procParams_, SEXP obsParams_, SEXP splineParams_, SEXP initStates_, SEXP nSteps_, SEXP nSimul_,  
                            SEXP obsInterval_, SEXP nObs_, SEXP T0_, SEXP splineSettings_);
                                                         
RcppExport SEXP malSirSimulCpp( SEXP procParams_, SEXP obsParams_, SEXP splineParams_, SEXP initStates_, SEXP population_, SEXP rainfall_,
                                SEXP nSteps_,  SEXP nSimul_,  SEXP obsInterval_, SEXP nObs_, SEXP T0_, SEXP splineSettings_ );
                                
RcppExport SEXP malStepPompCpp(SEXP params_, SEXP splineParams_, SEXP splineBasis_, SEXP initStates_,   SEXP rainfall_,  SEXP Pop_, 
                               SEXP nSteps_,  SEXP nSimul_,  SEXP obsInterval_, SEXP currIteration_);

#endif
