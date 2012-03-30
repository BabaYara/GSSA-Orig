/* Productivity.cpp is a routine that draws random series of the productivity
   shocks and simulates the corresponding series of the productivity levels; 
   see "Numerically Stable and Accurate Stochastic Simulation Approaches for 
   Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia Maliar and 
   Serguei Maliar, (2011), Quantitative Economics 2/2, 173-210 (henceforth, 
   JMM, 2011). This is a translation of the original Matlab code provided
   by Lilia Maliar and Serguei Maliar.

   This version: 28 March 2012.
   -------------------------------------------------------------------------
   Inputs:  "T" is the simulation length; T>=1;
            "N" is the number of countries; N>=1;
	    "a" is the initial condition for the productivity levels of
	    N countries; 1-by-N;
	    "rho" and "sigma" are the parameters of the model

   Output:  "a" is a pointer to the time series of the productivity levels
            of N countries; T-by-N
   -------------------------------------------------------------------------
   Copyright © 2012 by Eric M. Aldrich. All rights reserved. The code may
   be used, modified and redistributed under the terms provided in the file
   "License_Agreement.txt".
   -------------------------------------------------------------------------*/

#include "global.h"
#include "ran.h"

void Productivity(int T, int N, REAL* a, REAL sigma, REAL rho)
{

  int ix, jx;

  // Initialize random number generator
  Normaldev randn(0,1,100125);

  // A random draw of common-for-all-countries productivity shocks for
  // T periods; T-by-1
  REAL* EPSI = new REAL[T];
  for(ix = 0 ; ix < T ; ++ix) EPSI[ix] = randn.dev();

  // A random draw of country-specific productivity shocks for T periods
  // and N countries; T-by-N
  // Compute the error terms in the process for productivity level using
  // condition (4) in JMM (2011); T-by-N
  REAL* epsi = new REAL[T*N];
  for(ix = 0 ; ix < T ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      epsi[ix*N + jx] = (randn.dev() + EPSI[ix])*sigma;
    }
  }

  // Compute the next-period productivity levels using condition (4) in
  // JMM (2011); 1-by-N 
  for(ix = 0 ; ix < (T-1) ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      a[(ix+1)*N + jx] = pow(a[ix*N + jx], rho)*exp(epsi[(ix+1)*N + jx]);
    }
  }

}
