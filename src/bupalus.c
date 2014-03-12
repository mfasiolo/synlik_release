/* This file contains a model definition file for a stochastic parasitoid 
   interference model for bupalus. 
   The idea is to efficiently solve the  model for multiple 
   replicates, conditional on an unscaled noise vector.
*/
#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "summ_stats.h"


void bup_par(double *n,double *theta,double *e,int *burn_in,int *n_t, int *n_reps) {
/* Simulates `n_reps', length `n_t' replicates of bupalus-parasitoid model, 
   discarding an initial sequence of length `burn_in'.
   e_t must be length (burn_in+n_t)*n_reps
   n is n_t by n_reps
   This should be much faster than looping in R, even if all reps parallel.
*/
  int i,j;
  double r0,Ki,alpha,beta,sig_e,P,N,ner,pap;
  r0 = theta[0];
  Ki = 1/theta[1]; /* inverse carrying capacity */
  alpha = theta[2];
  beta = theta[3];
  sig_e = theta[4];

  /* following iterates the bupalus model... */
  for (j=0;j<*n_reps;j++) {
   P=.2;N=20.0;
   for (i=0;i< *burn_in;i++,e++) {
      ner = N*exp(r0*(1-N*Ki) + *e * sig_e); 
      pap = exp(-alpha * pow(P,beta));
      N = ner * pap;     /* hosts */
      P = ner * (1-pap); /* parasitoids */
    }
    for (i=0;i<*n_t;i++,e++,n++) {
      ner = N*exp(r0*(1-N*Ki) + *e * sig_e); 
      pap = exp(-alpha * pow(P,beta));
      *n = N = ner * pap;
      P = ner * (1-pap);
    }
  }
} /* end of bup_par */

