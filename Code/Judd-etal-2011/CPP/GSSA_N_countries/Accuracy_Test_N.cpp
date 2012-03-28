// Accuracy_Test_N.cpp is a routine for evaluating accuracy of solutions to the    
// multi-country model: it computes approximation errors in the optimality  
// conditions on a given set of points in the state space; see "Numerically 
// Stable and Accurate Stochastic Simulation Approaches for Solving Dynamic 
// Economic Models" by Kenneth L. Judd, Lilia Maliar and Serguei Maliar, 
// (2011), Quantitative Economics 2/2, 173-210 (henceforth, JMM, 2011)., and 
// Juillard and Villemot, (2011), "Multi-Country Real Business Cycle Models: 
// Accuracy Tests and Test Bench", Journal of Economic Dynamics and Control 
// 35, 178-185. This is a translation of the original Matlab code provided
// by Lilia Maliar and Serguei Maliar.
//
// This version: 28 March 2012.
// -------------------------------------------------------------------------
// Inputs:    "P" and "N" are the dimensions (rows and cols, respectively)
//            of matrices "k" and "a";
//            "k" and "a" are, respectively, current-period capital and 
//            productivity levels, in the given set of points on which the 
//            accuracy is tested; 
//            "bk" are the coefficients of the capital policy functions of N 
//            countries;
//            "IM" is the integration method for evaluating accuracy, 
//            IM=1,2,..,10=Gauss-Hermite quadrature rules with 1,2,...,10 
//            nodes in each dimension, respectively, 
//            IM=11=Monomial rule with 2N nodes,
//            IM=12=Monomial rule with 2N^2+1 nodes;
//            "alpha", "gam", "phi", "beta", "A", "tau", "rho" and "vcv"
//            are the parameters of the model;
//            "discard" is the number of data points to discard 
//
// Outputs:   "Errors_mean" and "Errors_max" are, respectively, the mean and
//            maximum approximation errors across all optimality conditions;
//            "Errors_max_EE", "Errors_max_MUC", "Errors_max_MUL", and 
//            "Errors_max_RC" are the maximum approximation errors disaggre- 
//            gated by optimality conditions 
// -------------------------------------------------------------------------
// Copyright © 2012 by Eric M. Aldrich. All rights reserved. The code may
// be used, modified and redistributed under the terms provided in the
// file "License_Agreement.txt".
// -------------------------------------------------------------------------

#include "global.h"
#include <math.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <typeinfo>
#include <gsl/gsl_cblas.h>

using namespace std;

void Accuracy_Test_N(int P, int N, REAL* k, REAL* a, REAL* bk, int D,
		     int IM, REAL alpha, REAL gam, REAL delta, REAL beta,
		     REAL A, REAL tau, REAL rho, REAL* vcv, int discard,
		     REAL& Errors_mean, REAL& Errors_max, REAL& time_test)
{

  int ix, jx, px;

  // Start counting time needed to run the test
  clock_t tic = clock();

  //==========================================================================
  // 2. Integration method for evaluating accuracy
  //==========================================================================
  REAL* epsi_nodes = NULL;
  REAL* weight_nodes = NULL;
  int n_nodes;
  if((IM>=1) && (IM<=10)){
    // Compute the number of integration nodes, n_nodes, integration
    // nodes, epsi_nodes, and integration weights, weight_nodes, for Gauss-
    // Hermite quadrature integration rule with IM nodes in each dimension
    n_nodes = pow(IM, N);
    epsi_nodes = (REAL*)realloc(epsi_nodes, n_nodes*N*sizeof(REAL));
    weight_nodes = (REAL*)realloc(weight_nodes, n_nodes*sizeof(REAL));
    GH_Quadrature(IM, N, vcv, n_nodes, epsi_nodes, weight_nodes);

  } else if(IM == 11){
    // Monomial integration rule with 2N nodes
    n_nodes = 2*N;
    epsi_nodes = (REAL*)realloc(epsi_nodes, n_nodes*N*sizeof(REAL));
    weight_nodes = (REAL*)realloc(weight_nodes, n_nodes*sizeof(REAL));
    Monomials_1(N, vcv, n_nodes, epsi_nodes, weight_nodes);

  } else if(IM == 12){
    // Monomial integration rule with 2N^2+1 nodes
    n_nodes = 2*pow(N, 2) + 1;
    epsi_nodes = (REAL*)realloc(epsi_nodes, n_nodes*N*sizeof(REAL));
    weight_nodes = (REAL*)realloc(weight_nodes, n_nodes*sizeof(REAL));
    Monomials_2(N, vcv, n_nodes, epsi_nodes, weight_nodes);

  }

  //==========================================================================
  // 3. Polynomial bases for the test
  //==========================================================================

  // Matrix which holds data for polynomial bases
  REAL X1[P*2*N];
  for(ix = 0 ; ix < P ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      X1[ix*2*N + jx] = k[ix*N + jx];
    }
    for(jx = N ; jx < 2*N ; ++jx){
      X1[ix*2*N + jx] = a[ix*N + jx - N];
    }
  }

  // Form a complete ordinary polynomial of degree D on given set of points
  int n_cols;
  REAL* poly_X = Ord_Polynomial_N(X1, P, 2*N, D, n_cols);

  //==========================================================================
  // 4. Given the solution for capital, compute consumption on the given set
  // of points
  //==========================================================================
  
  int nE_cols = 4*N+1;
  REAL k0[N], a0[N], poly_X0[n_cols], k1[N], C0, c0[N];
  REAL a1[n_nodes*N], k1_dupl[n_nodes*N], X1_dupl[n_nodes*2*N], k2[n_nodes*N];
  REAL C1[n_nodes], c1[n_nodes*N], MUC0j[N], MUC1j[n_nodes*N];
  REAL lambda0, lambda1[N], numerator, denominator, Errors[P*nE_cols];
  REAL* poly_X_dupl;
  for(px = 0 ; px < P ; ++px){ // For each given point, ...     
    
    // Display the point (with the purpose of monitoring the progress)
    //cout << "P = " << px << endl;
                                     
    //========================================================================
    // 4.1 Variables in point p
    //========================================================================
    
    // N capital stocks and productivity levels of period p
    for(ix = 0 ; ix < N ; ++ix){
      k0[ix] = k[px*N+ix];
      a0[ix] = a[px*N+ix];
    }
    // Complete (second-degree) polynomial bases at t
    for(ix = 0 ; ix < n_cols ; ++ix){
      poly_X0[ix] = poly_X[px*n_cols+ix];
    }          

    //========================================================================
    // 4.2 Capital and consumption choices at t
    //========================================================================

    // Compute a row-vector of capital of period t+1 (chosen at t) using
    // the corresponding capital policy functions; 1-by-N
    if(typeid(realtype) == typeid(singletype)){
      cblas_sgemv(CblasRowMajor, CblasTrans, n_cols, N, 1.0, (float*)bk,
		  N, (float*)poly_X0, 1, 0.0, (float*)k1, 1);
    } else if(typeid(realtype) == typeid(doubletype)){
      cblas_dgemv(CblasRowMajor, CblasTrans, n_cols, N, 1.0, (double*)bk,
		  N, (double*)poly_X0, 1, 0.0, (double*)k1, 1);
    }
        
    // C is computed by summing up individual consumption, which in turn, is 
    // found from the individual budget constraints; 1-by-1
    C0 = 0.0;
    for(ix = 0 ; ix < N ; ++ix){
      C0 += A*pow(k0[ix], alpha)*a0[ix] - k1[ix] + k0[ix]*(1-delta);
    }

    // Individual consumption is the same for all countries; 1-by-N
    for(ix = 0 ; ix < N ; ++ix){
      c0[ix] = C0/N;
    }

    //========================================================================
    // 4.3 Capital and consumption choices at t+1
    //========================================================================

    // Compute the next-period productivity levels in each integration node
    // using condition (?) in the online appendix; n_nodes-by-N
    for(ix = 0 ; ix < n_nodes ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	a1[ix*N+jx] = pow(a0[jx], rho)*exp(epsi_nodes[ix*N+jx]);
      }
    }

    // Duplicate k1 n_nodes times to create a matrix with n_nodes identical
    // rows; n_nodes-by-N 
    for(ix = 0 ; ix < n_nodes ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	k1_dupl[ix*N+jx] = k1[jx];
      }
    }

    // Matrix which holds data for polynomial bases
    for(ix = 0 ; ix < n_nodes ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	X1_dupl[ix*2*N + jx] = k1_dupl[ix*N + jx];
      }
      for(jx = N ; jx < 2*N ; ++jx){
	X1_dupl[ix*2*N + jx] = a1[ix*N + jx - N];
      }
    }

    // Form a complete polynomial of degree D (at t+1) in the given point 
    poly_X_dupl = Ord_Polynomial_N(X1_dupl, n_nodes, 2*N, D, n_cols);
    
    // Compute capital of period t+2 (chosen at t+1) using the second-
    // degree capital policy functions; n_nodes-by-N  
    if(typeid(realtype) == typeid(singletype)){
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_nodes, N,
		  n_cols, 1.0, (float*)poly_X_dupl, n_cols, (float*)bk,
		  N, 0.0, (float*)k2, N);
    } else if(typeid(realtype) == typeid(doubletype)){
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_nodes, N,
		  n_cols, 1.0, (double*)poly_X_dupl, n_cols, (double*)bk,
		  N, 0.0, (double*)k2, N);
    }
    
    // Aggregate consumption is computed by summing up individual consumption, 
    // which in turn, is found from the individual budget constraints; 
    // n_nodes-by-1 
    for(ix = 0 ; ix < n_nodes ; ++ix){
      C1[ix] = 0.0;
      for(jx = 0 ; jx < N ; ++jx){
	C1[ix] += A*pow(k1_dupl[ix*N+jx], alpha)*a1[ix*N+jx] - k2[ix*N+jx] + k1_dupl[ix*N+jx]*(1-delta);
      }
    }

    // Individual consumption is the same for all countries; n_nodes-by-N
    for(ix = 0 ; ix < n_nodes ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	c1[ix*N + jx] = C1[ix]/N;
      }
    }

  //==========================================================================
  // 5. Approximation errors in point p
  //==========================================================================

    //========================================================================
    // 5.1 Lagrange multiplier associated with the aggregate resource
    // constraint
    //========================================================================

    // Compute a country's marginal utility of consumption multiplied 
    // by its welfare weight
    lambda0 = 0.0;
    for(ix = 0 ; ix < N ; ++ix){
      MUC0j[ix] = tau*pow(c0[ix], -gam); 
      lambda0 += MUC0j[ix];
    }

    // An optimality condition w.r.t. consumption of period t equates 
    // the Lagrange multiplier of the aggregate resource constraint of 
    // period t and each country's marginal utility of consumption 
    // multiplied by its welfare weight; to infer the Lagrange multiplier,  
    // we average across N countries; 1-by-1
    lambda0 /= N;
        
    for(ix = 0 ; ix < n_nodes ; ++ix){
      lambda1[ix] = 0.0;
      for(jx = 0 ; jx < N ; ++jx){

	// Compute a country's marginal utility of consumption multiplied 
	// by its welfare weight
	MUC1j[ix*N+jx] = tau*pow(c1[ix*N+jx], -gam);
	lambda1[ix] += MUC1j[ix*N+jx];

      }

      // Similarly, the Lagrange multiplier of the aggregate resource 
      // constraint of period t+1 is equal to a country's marginal utility 
      // of consumption multiplied by its welfare weight; to infer the 
      // Lagrange multiplier, we average across N countries; 1-by-n_nodes
      lambda1[ix] /= N;

    }  

    //========================================================================
    // 5.2 Unit-free Euler-equation errors
    //========================================================================

    for(jx = 0 ; jx < N ; ++jx){
      Errors[px*nE_cols+jx] = 1.0;
      for(ix = 0 ; ix < n_nodes ; ++ix){
        // A unit-free Euler-equation approximation error of country j
	Errors[px*nE_cols+jx] -= weight_nodes[ix]*(beta*(lambda1[ix]/lambda0)*(1-delta+alpha*A*pow(k1[jx], alpha-1)*a1[ix*N+jx]));
      }
    }
        
    //========================================================================
    // 5.2 Unit-free errors in the optimality conditions w.r.t. consumption
    //======================================================================== 
    for(jx = 0 ; jx < N ; ++jx){
      Errors[px*nE_cols+N+jx] = 1-lambda0/(tau*pow(c0[jx], -gam)); 
      // A unit-free approximation error in the optimality condition w.r.t. 
      // consumption of country j (this condition equates marginal utility 
      // of consumption, multiplied by the welfare weight, and the 
      // Lagrange multiplier of the aggregate resource constraint)
    }
        
    //========================================================================
    // 5.3 Unit-free errors in the optimality conditions w.r.t. labor 
    //======================================================================== 
    for(jx = 0 ; jx < N ; ++jx){    
      Errors[px*nE_cols + 2*N + jx] = 0.0;
      // These errors  are zero by construction 
    }
        
    //========================================================================
    // 5.4 Unit-free approximation error in the aggregate resource constraint
    //======================================================================== 
    numerator = 0.0;
    denominator = 0.0;
    for(ix = 0 ; ix < N ; ++ix){
      numerator += c0[ix] + k1[ix] - k0[ix]*(1-delta);
      denominator += A*pow(k0[ix], alpha)*a0[ix];
    }
    Errors[px*nE_cols+3*N] = 1.0 - numerator/denominator;
    // This error is a unit-free expression of the resource constraint  
    // (?) in the online appendix

    //========================================================================
    // 5.5 Approximation errors in the capital-accumulation equation
    //======================================================================== 
    for(jx = 0 ; jx < N ; ++jx){    
      Errors[px*nE_cols + 3*N + 1 + jx] = 0.0;
      // These errors are always zero by construction 
    }
        
    // For this model, GSSA produces zero errors (with machine accuracy)
    // in all the optimality conditions except of the Euler equation. 
    // Here, the errors in all the optimality conditions are introduced 
    // to make our accuracy measures comparable to those in the February 
    // 2011 special issue of the Journal of Economic Dynamics and Control 

  }

  //==========================================================================
  // 6. Mean and maximum approximation errors computed after discarding the
  // first "discard" observations
  //==========================================================================

  //==========================================================================
  // 6.1 Approximation errors across all the optimality conditions
  //========================================================================== 
  
  // Average and Maximum absolute approximation error
  int ind_max;
  if(typeid(realtype) == typeid(singletype)){
    Errors_mean = cblas_sasum((P-discard)*nE_cols, (float*)Errors+discard*nE_cols, 1);
    ind_max = cblas_isamax((P-discard)*nE_cols, (float*)Errors+discard*nE_cols, 1);
  } else if(typeid(realtype) == typeid(doubletype)){
    Errors_mean = cblas_dasum((P-discard)*nE_cols, (double*)Errors+discard*nE_cols, 1);
    ind_max = cblas_idamax((P-discard)*nE_cols, (double*)Errors+discard*nE_cols, 1);
  }
  Errors_mean = log10(Errors_mean/((P-discard)*nE_cols));
  Errors_max = log10(fabs(Errors[discard*nE_cols + ind_max]));
 
  //==========================================================================
  // 6.2 Maximum approximation errors disaggregated by the optimality
  // conditions
  // This is currently not supported.
  //========================================================================== 

  //Errors_max_EE = log10(max(max(abs(Errors(1+discard:end,1:N)))));
  // Across N Euler equations

  //Errors_max_MUC = log10(max(max(abs(Errors(1+discard:end,N+1:2*N)))));    
  // Across N optimality conditions w.r.t. consumption (conditions on   
  // marginal utility of consumption, MUC)

  //Errors_max_MUL = log10(max(max(abs(Errors(1+discard:end,2*N+1:3*N)))));    
  // Across N optimality conditions w.r.t. labor (conditions on marginal 
  // utility of labor, MUL)
  
  //Errors_max_RC = log10(max(max(abs(Errors(1+discard:end,3*N+1)))));    
  // In the aggregate resource constraint 
  
  //==========================================================================
  // 7. Time needed to run the test
  //==========================================================================
  time_test = (clock() - tic)/(REAL)CLOCKS_PER_SEC;

}
