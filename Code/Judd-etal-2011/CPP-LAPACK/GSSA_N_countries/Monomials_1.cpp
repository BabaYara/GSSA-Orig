// Monomials_1.cpp is a routine that constructs integration nodes and weights 
// under N-dimensional monomial (non-product) integration rule with 2N nodes; 
// see the Supplement to "Numerically Stable and Accurate Stochastic Simulation 
// Approaches for Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia 
// Maliar and Serguei Maliar, (2011), Quantitative Economics 2/2, 173-210 
// (henceforth, JMM, 2011). This is a translation of the original Matlab code
// provided by Lilia Maliar and Serguei Maliar.
//
// This version: 28 March 2012.
// ---------------------------------------------------------------------------
// Inputs:  "N" is the number of random variables; N>=1;
//          "vcv" is the variance-covariance matrix; N-by-N

// Outputs: "n_nodes" is the total number of integration nodes; 2*N;
//          "epsi_nodes" are the integration nodes; n_nodes-by-N;
//          "weight_nodes" are the integration weights; n_nodes-by-1
// ---------------------------------------------------------------------------
// Copyright © 2012 by Eric M. Aldrich. All rights reserved. The code may be
// used, modified and redistributed under the terms provided in the file
// "License_Agreement.txt".
// ---------------------------------------------------------------------------

#include "global.h"
#include <math.h>
#include <typeinfo>
#include <iostream>

using namespace std;

void Monomials_1(int N, REAL* vcv, int n_nodes, REAL* epsi_nodes,
		 REAL* weight_nodes)
{

  int ix, jx;

  //==========================================================================
  // 1. N-dimensional integration nodes for N uncorrelated random variables
  //    with zero mean and unit variance
  //==========================================================================

  // A supplementary matrix for integration nodes; n_nodes-by-N
  REAL z1[n_nodes*N];
  for(ix = 0 ; ix < n_nodes ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      z1[ix + jx*n_nodes] = 0.0;
    }
  }
                       
  // In each node, random variable i takes value either 1 or -1, and all
  // other variables take value 0
  for(ix = 0 ; ix < N ; ++ix){
    for(jx = (2*ix+1) ; jx < 2*(ix+1) ; ++jx){
      z1[2*ix + ix*n_nodes] = 1.0;
      z1[(2*ix+1) + ix*n_nodes] = -1.0;
    }
  }
  // For example, for N = 2, z1 = [1 0; -1 0; 0 1; 0 -1]

  //==========================================================================
  // 2. N-dimensional integration nodes and weights for N correlated random
  //    variables with zero mean and variance-covaraince matrix vcv
  //==========================================================================

  // Cholesky decomposition of the variance-covariance matrix
  int info;
  dpotrf("U", &N, vcv, &N, &info);
  for(ix = 1 ; ix < N ; ++ix){
    for(jx = 0 ; jx < ix ; ++jx){
      vcv[ix + jx*N] = 0.0;
    }
  }
                                 
  // Integration nodes; see condition (B.7) in the Supplement
  // to JMM (2011); n_nodes-by-N
  if(typeid(realtype) == typeid(singletype)){
    float f_alpha = sqrt(N);
    float f_beta = 0.0;
    sgemm("N", "N", &n_nodes, &N, &N, &f_alpha, (float*)z1, &n_nodes,
	  (float*)vcv, &N, &f_beta, (float*)epsi_nodes, &n_nodes);
  } else if(typeid(realtype) == typeid(doubletype)){
    double d_alpha = sqrt(N);
    double d_beta = 0.0;
    dgemm("N", "N", &n_nodes, &N, &N, &d_alpha, (double*)z1, &n_nodes,
	  (double*)vcv, &N, &d_beta, (double*)epsi_nodes, &n_nodes);
  }

  //==========================================================================
  // 3. Integration weights
  //==========================================================================

  // Integration weights are equal for all integration nodes; n_nodes-by-1;
  // the weights are the same for the cases of correlated and uncorrelated
  // random variables
  for(ix = 0 ; ix < n_nodes ; ++ix){
    weight_nodes[ix] = 1.0/n_nodes; 
  }

}
