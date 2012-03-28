// Monomials_2.cpp is a routine that constructs integration nodes and weights 
// under N-dimensional monomial (non-product) integration rule with 2N^2+1 
// nodes; see supplement to "Numerically Stable and Accurate Stochastic 
// Simulation Approaches for Solving Dynamic Economic Models" by Kenneth L. 
// Judd, Lilia Maliar and Serguei Maliar, (2011), Quantitative Economics 2/2, 
// 173-210 (henceforth, JMM, 2011). This is a translation of the original
// Matlab code provided by Lilia Maliar and Serguei Maliar.
//
// This version: 28 March 2012.
// ---------------------------------------------------------------------------
// Inputs:  "N" is the number of random variables; N>=1;
//          "vcv" is the variance-covariance matrix; N-by-N;

// Outputs: "n_nodes" is the total number of integration nodes; 2*N^2+1;
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
#include <gsl/gsl_cblas.h>
#include <iostream>
#include <iomanip>

using namespace std;

void Monomials_2(int N, REAL* vcv, int n_nodes, REAL* epsi_nodes,
		 REAL* weight_nodes)
{

  int ix, jx, px, qx;

  //==========================================================================
  // 1. N-dimensional integration nodes for N uncorrelated random variables
  //    with zero mean and unit variance
  //==========================================================================

  //==========================================================================
  // 1.1 The origin point
  //==========================================================================
  // A supplementary matrix for integration nodes: the origin point
  REAL z0[N];
  for(ix = 0 ; ix < N ; ++ix) z0[ix] = 0.0;

  //==========================================================================
  // 1.2 Deviations in one dimension
  //==========================================================================
  // A supplementary matrix for integration nodes; 2N-by-N
  REAL z1[2*N*N];
  for(ix = 0 ; ix < 2*N ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      z1[ix*N + jx] = 0.0;
    }
  }
                       
  // In each node, random variable i takes value either 1 or -1, and all
  // other variables take value 0
  for(ix = 0 ; ix < N ; ++ix){
    for(jx = (2*ix+1) ; jx < 2*(ix+1) ; ++jx){
      z1[2*ix*N + ix] = 1.0;
      z1[(2*ix+1)*N + ix] = -1.0;
    }
  }
  // For example, for N = 2, z1 = [1 0; -1 0; 0 1; 0 -1]

  //==========================================================================
  // 1.3 Deviations in two dimensions
  //==========================================================================
  // A supplementary matrix for integration nodes; 2N(N-1)-by-N
  REAL z2[2*N*(N-1)*N];
  for(ix = 0 ; ix < 2*N*(N-1) ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      z2[ix*N + jx] = 0.0;
    }
  }

  // In each node, a pair of random variables (p,q) takes either values (1,1)
  // or (1,-1) or (-1,1) or (-1,-1), and all other variables take value 0
  ix = 0;
  for(px = 0 ; px < (N-1) ; ++px){
    for(qx = (px+1) ; qx < N ; ++qx){
      z2[4*ix*N + px] = 1.0;
      z2[(4*ix+1)*N + px] = -1.0;
      z2[(4*ix+2)*N + px] = 1.0;
      z2[(4*ix+3)*N + px] = -1.0;
      z2[4*ix*N + qx] = 1.0;
      z2[(4*ix+1)*N + qx] = 1.0;
      z2[(4*ix+2)*N + qx] = -1.0;
      z2[(4*ix+3)*N + qx] = -1.0;
      ++ix;
    }
  }
  // For example, for N = 2, z2 = [1 1;1 -1;-1 1;-1 1]


  //==========================================================================
  // 2. N-dimensional integration nodes and weights for N correlated random
  //    variables with zero mean and variance-covaraince matrix vcv
  //==========================================================================

  // Cholesky decomposition of the variance-covariance matrix
  REAL* vcv_chol = new REAL[N*N];
  cholesky(vcv, N, vcv_chol);
                                 
  // Integration nodes; see condition (B.8) in the Supplement
  // to JMM (2011); n_nodes-by-N
  REAL z1_alt[2*N*N];
  REAL z2_alt[2*N*(N-1)*N];
  if(typeid(realtype) == typeid(singletype)){
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 2*N, N, N,
		sqrt(N+2), (float*)z1, N, (float*)vcv_chol, N, 0.0,
		(float*)z1_alt, N);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 2*N*(N-1), N, N,
		sqrt((N+2)/2), (float*)z2, N, (float*)vcv_chol, N, 0.0,
		(float*)z2_alt, N);
  } else if(typeid(realtype) == typeid(doubletype)){
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 2*N, N, N,
		sqrt(N+2), (double*)z1, N, (double*)vcv_chol, N, 0.0,
		(double*)z1_alt, N);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 2*N*(N-1), N, N,
		sqrt((N+2)/2), (double*)z2, N, (double*)vcv_chol, N, 0.0, 
		(double*)z2_alt, N);
  }
  for(ix = 0 ; ix < n_nodes ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      if(ix == 0){
	epsi_nodes[ix*N + jx] = z0[jx];
      } else if(ix >= 1 & ix <= 2*N){
	epsi_nodes[ix*N + jx] = z1_alt[(ix-1)*N+jx];
      } else {
	epsi_nodes[ix*N + jx] = z2_alt[(ix-2*N-1)*N+jx];
      }
    }
  }

  //==========================================================================
  // 3. Integration weights
  //==========================================================================

  // See condition in (B.8) in the Supplement to JMM (2011); n_nodes-by-1;
  // the weights are the same for the cases of correlated and uncorrelated
  // random variables
  weight_nodes[0] = 2.0/(N+2);
  for(ix = 1 ; ix <= 2*N ; ++ix){
    weight_nodes[ix] = (4-N)/(2*pow(N+2,2));
  }
  for(ix = (2*N+1) ; ix <= n_nodes ; ++ix){
    weight_nodes[ix] = 1/pow(N+2,2);
  }

}
