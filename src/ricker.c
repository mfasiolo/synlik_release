/* This file contains a model definition file for the Ricker model.
   The idea is to efficiently solve the Ricker model for multiple 
   replicates, conditional on an unscaled noise vector.
*/
#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "summ_stats.h"


void ricker(double *n,double *theta,double *e,int *burn_in,int *n_t, int *n_reps,double *n0) {
/* Simulates `n_reps', length `n_t' replicates of Ricker model, discarding an initial 
   sequence of length `burn_in'.
   e_t must be length (burn_in+n_t)*n_reps
   n is n_t by n_reps
   This should be much faster than looping in R, even if all reps parallel.
*/
  int i,j;
  double log_r,sig_e,x;
  log_r = theta[0];
  sig_e = theta[1];
  /* following iterates the Ricker map on log populations */
  for (j=0;j<*n_reps;j++) {
    for (x= *n0,i=0;i< *burn_in;i++,e++) x += log_r - exp(x) + *e * sig_e;     
    for (i=0;i<*n_t;i++,e++,n++) { x += log_r - exp(x) + *e * sig_e; *n = x; }
  }
} /* end of Ricker */

/* equivalent to...

rickR <- function(theta,e,n0=1,burn.in=50) {
  log.r <- theta[1]
  sig.e <- theta[2]
  n.reps <- ncol(e)
  n.t <- nrow(e)-burn.in
  x <- rep(n0,n.reps)
  for (i in 1:burn.in) x <- x + log.r - exp(x) + e[i,] * sig.e
  n <- matrix(0,n.t,n.reps);n[1,] <- x
  for (i in 2:n.t) n[i,] <- n[i-1,] + log.r - exp(n[i-1,]) + e[i+burn.in,] * sig.e
  n
}


*/
