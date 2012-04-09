// Num_Stab_Approx.cpp is a routine that implements the approximation methods 
// described in "Numerically Stable and Accurate Stochastic Simulation 
// Approaches for Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia 
// Maliar and Serguei Maliar, (2011), Quantitative Economics 2/2, 173-210 
// (henceforth, JMM, 2011). This is a translation of the original Matlab code
// provided by Lilia Maliar and Serguei Maliar.
//
// This version: 28 March 2012.
// ---------------------------------------------------------------------------
// Inputs:  "T" and "n" are the dimensions (rows and cols, repsectively) of
//          matrix "X"
//          "X" is a matrix of dependent variables in a regression, T-by-n,
//          where n corresponds to the total number of coefficients in the 
//          original regression (i.e. with unnormalized data);
//          "Y" is a matrix of independent variables, T-by-N; 
//          "RM" is the regression (approximation) method, RM=1,...,8:  
//          1=OLS,          2=LS-SVD,    3=LAD-PP,    4=LAD-DP, 
//          5=RLS-Tikhonov, 6=RLS-TSVD,  7=RLAD-PP,   8=RLAD-DP;
//          "penalty"  is a parameter determining the value of the regulari-
//          zation parameter for a regularization methods, RM=5,6,7,8;
//          "normalize"  is the option of normalizing the data, 
//          0=unnormalized data,  1=normalized data                

// Outputs: "B" is a matrix of the regression coefficients 
// ---------------------------------------------------------------------------
// Copyright © 2012 by Eric M. Aldrich. All rights reserved. The code may be
// used, modified and redistributed under the terms provided in the file
// "License_Agreement.txt".
// ---------------------------------------------------------------------------

#include "global.h"
#include <math.h>
#include <typeinfo>
#include <iostream>
#include <iomanip>

using namespace std;

REAL* Num_Stab_Approx(int T, int n, REAL* X, int N, REAL* Y, int RM,
		      int penalty, int normalize)
{

  int ix, jx;

  // Storage for regression coefficients
  REAL* B = new REAL[n*N];

  //==========================================================================
  // 2. Normalize the data
  //==========================================================================

  int Tplus1 = T+1;
  int n1;
  if(normalize == 1){
    n1 = n-1;
  } else {
    n = n;
  }

  REAL* X1 = new REAL[Tplus1*n1];
  REAL* Y1 = new REAL[T*N];
  REAL* B1 = new REAL[n1*N];
  REAL* X_means = new REAL[n1];
  REAL* X_sds = new REAL[n1];
  REAL* Y_means = new REAL[N];
  REAL* Y_sds = new REAL[N];

  if((normalize == 1) | (RM >= 5)){ 

    // If we use normalized data or consider a regularization method, ...

    // eliminate the first column of X
    REAL X_nocol[Tplus1*n1];
    for(ix = 0 ; ix < Tplus1 ; ++ix){
      for(jx = 0 ; jx < n1 ; ++jx){
	X_nocol[ix + jx*Tplus1] = X[ix + (jx+1)*Tplus1];
      }
    }

    // Center and scale X1
    matrix_normalize(T, n1, X_nocol, T+1, X1, X_means, X_sds);

    // Center and scale Y
    matrix_normalize(T, N, Y, T, Y1, Y_means, Y_sds);

  } else {

    // If we use unnormalized data, ...
    X1 = X;    // Leave X without changes
    Y1 = Y;    // Leave Y without changes

  }

  //==========================================================================
  // 3. Regression methods
  //==========================================================================

  //==========================================================================
  // 3.1 OLS
  //==========================================================================

  if(RM == 1){

    // If the regression method is OLS, ...

    // storage for linear regression variables
    float f_wkopt;
    double d_wkopt;
    REAL* work = NULL;
    int info;
    int lwork;
    
    // Compute B using the formula of the OLS estimator; note that
    // all the regressions, j=1,...,N, are run at once
    lwork = -1;
    if(typeid(realtype) == typeid(singletype)){
      sgels("N", &T, &n1, &N, (float*)X1, &Tplus1, (float*)Y1, &T, &f_wkopt, &lwork, &info);
      lwork = (int)f_wkopt;
    } else if(typeid(realtype) == typeid(doubletype)){
      dgels("N", &T, &n1, &N, (double*)X1, &Tplus1, (double*)Y1, &T, &d_wkopt, &lwork, &info);
      lwork = (int)d_wkopt;
    }
    work = (REAL*)realloc(work, lwork*sizeof(REAL));
    if(typeid(realtype) == typeid(singletype)){
      sgels("N", &T, &n1, &N, (float*)X1, &Tplus1, (float*)Y1, &T, (float*)work, &lwork, &info);
    } else if(typeid(realtype) == typeid(doubletype)){
      dgels("N", &T, &n1, &N, (double*)X1, &Tplus1, (double*)Y1, &T, (double*)work, &lwork, &info);
    }
    for(ix = 0 ; ix < n1 ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	B1[ix+jx*n1] = Y1[ix+jx*T];
      }
    }

  }

  //==========================================================================
  // 3.2 LS-SVD
  //==========================================================================
  else if(RM == 2){

    // If the regression method is LS-SVD, ...

    // storage for linear regression variables
    REAL S[n];
    float f_rcond = -1.0;
    double d_rcond = -1.0;
    int rank;
    float f_wkopt;
    double d_wkopt;
    REAL* work = NULL;
    int info;
    int lwork;
    
    // Compute B using the formula of the OLS estimator; note that
    // all the regressions, j=1,...,N, are run at once
    lwork = -1;
    if(typeid(realtype) == typeid(singletype)){
      sgelss(&T, &n1, &N, (float*)X1, &Tplus1, (float*)Y1, &T, (float*)S, &f_rcond, &rank, &f_wkopt, &lwork, &info);
      lwork = (int)f_wkopt;
    } else if(typeid(realtype) == typeid(doubletype)){
      dgelss(&T, &n1, &N, (double*)X1, &Tplus1, (double*)Y1, &T, (double*)S, &d_rcond, &rank, &d_wkopt, &lwork, &info);
      lwork = (int)d_wkopt;
    }
    work = (REAL*)realloc(work, lwork*sizeof(REAL));
    if(typeid(realtype) == typeid(singletype)){
      sgelss(&T, &n1, &N, (float*)X1, &Tplus1, (float*)Y1, &T, (float*)S, &f_rcond, &rank, (float*)work, &lwork, &info);
    } else if(typeid(realtype) == typeid(doubletype)){
      dgelss(&T, &n1, &N, (double*)X1, &Tplus1, (double*)Y1, &T, (double*)S, &d_rcond, &rank, (double*)work, &lwork, &info);
    }
    for(ix = 0 ; ix < n1 ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	B1[ix+jx*n1] = Y1[ix+jx*T];
      }
    }

  }
  // For each regression j, RM 3, 4, 7, 8, we  solve a linear programming 
  // problem written in MATLAB:
  //                       min  f'*xlp 
  //                       s.t. Aeq*xlp=beq       (*)
  //                            Aineq*xlp<=bineq  (**) 
  //                            LB=xlp<=UB
  // where xlp is the vector of unknowns; Aeq is the matrix of coefficients 
  // in the equality restriction (*); Aineq is the matrix of coefficients in the
  // inequality restriction (**); beq and bineq are vectors of constants in the 
  // right sides of (*) and (**), respectively; LB means "lower bound", UB  
  // means "upper bound"

  //==========================================================================
  // 3.3 LAD-PP
  //==========================================================================
  /*  
  else if (RM == 3) {          // If the regression method is LAD-PP, ...
    // We solve the linear programming problem (27)-(29) in JMM (2011). 
    // xlp=[B(:,j); ups_plus; ups_minus] where B(:,j) is the vector of coeffi-
    // cients in the regression j, ups_plus and ups_minus are the deviations; 
    // f'=[0,...,0,1,...,1] where zeros correspond to the coefficients in the 
    // objective function on B and ones correspond to the coefficients on 
    // ups_plus and ups_minus
    
    // Initialize tableau matrix
    int nA = 2*(n1+T)+1; // number of columns in tableau matrix
    REAL* A = new REAL[(T+2)*nA];

    // First row of tableau matrix has a zero followed by the negative of
    // the objective function coefficients f'=[0,...,0,1,...,1]
    A[0] = 0.0;
    for(ix = 1 ; ix < (2*n1+1) ; ++ix) A[ix] = 0.0;
    for(ix = (2*n1+1) ; ix < nA ; ++ix) A[ix] = -1.0;

    // The first column of the remaining rows will be the constraint values
    // Y[,i], filled in later. The other columns consist of the matrix
    // [X1, -X1, I, -I], where I is the TxT identity matrix
    for(ix = 1 ; ix < Tplus1 ; ++ix){
      for(jx = 1 ; jx < (n1+1) ; ++jx) A[ix*nA+jx] = X1[(ix-1)*n1+jx-1];
      for(jx = (n1+1) ; jx < (2*n1+1) ; ++jx) A[ix*nA+jx] = -X1[(ix-1)*n1+jx-1-n1];
      for(jx = (2*n1+1) ; jx < nA ; ++jx){
	if((jx - 2*n1) == ix){
	  A[ix*nA+jx] = 1.0;
	} else if((jx - 2*n1 - T) == ix) {
	  A[ix*nA+jx] = -1.0;
	} else {
	  A[ix*nA+jx] = 0.0;
	}
      }
    }

    // Storage for simplex algorithm
    int icase;
    int izrov[nA-1];
    int iposv[T];

    // Run simplex for each column of Y
    for(jx = 0 ; jx < N ; ++jx){ // For each regression j, ...

      // fill in the first column of the tableau matrix for rows 1 to T
      for(ix = 1 ; ix < Tplus1 ; ++ix){
	A[ix*nA] = Y1[(ix-1)*N+jx];
      }
      
      // Simplex
      simplx(A, T, nA-1, 0, 0, T, icase, izrov, iposv);

      // Store the regression coefficients for the regression j,
      // xlp(1:n1,1), into the matrix B 
      for(ix = 0 ; ix < T ; ++ix){
	if(iposv[ix] < n1) {
	  B[iposv[ix]*N+jx + N] += A[(ix+1)*nA];
	} else if(iposv[ix] < 2*n1) {
	  B[(iposv[ix]-n1)*N+jx + N] += -A[(ix+1)*nA];
	}
      }
    }
    for(jx = 0 ; jx < T ; ++jx){
      if(iposv[jx] < 2*n1){
	cout << "Variable " << iposv[jx] << " is equal to " << A[(jx+1)*nA] << endl; 
      }
    }
  }

  
 // 3.4 LAD-DP
 //-----------
 elseif RM == 4;          // If the regression method is LAD-DP, ...
 // We solve the linear programming problem (30)-(32) in JMM (2011). 
 // xlp=[q] where q is a vector of unknowns in (30)-(32) of JMM (2011)
   LB = -ones(1,T);      // Lower bound on q is -1
   UB = ones(1, T);      // Upper bound on q is 1
   Aeq =  X1';           // n1-by-T    
   beq = zeros(n1,1);    // n1-by-1
   B = zeros(size(X1,2),N);
                         // Allocate memory for the regression coefficients;
                         // n1-by-N
   for j = 1:N           // For each regression j, ...
       f = -Y1(:,j);     // Our objective function is -Y1'(:,j)*q (minus sign 
                         // appears because (30)-(32) in JMM (2011) is a 
                         // minimization problem)
       [xlp,fval,exitflag,output,lambda] = linprog(f,[],[],Aeq,beq,LB,UB,[]);
                         // Find the values of the Lagrange multipliers, 
                         // lambda, on all the constraints
        B(:,j) = lambda.eqlin;
                         // Coefficients for the regression j, B(:,j),
                         // are equal to the Lagrange multipliers on the  
                         // equality constraint (*)
   end
  */

  //==========================================================================
  // 3.5 RLS-Tikhonov
  //==========================================================================

  else if(RM == 5){

    // If the regression method is RLS-Tikhonov, ...
    // storage for linear regression variables
    REAL XX[n1*n1];
    REAL XY[n1*N];
    
    // Use the formula (22) in JMM (2011) where the regularization parameter
    // is T/n1*10^-penalty
    // First compute X'X and X'Y.
    if(typeid(realtype) == typeid(singletype)){
      float f_alpha = 1.0;
      float f_beta = 0.0;
      sgemm("T", "N", &n1, &n1, &T, &f_alpha, (float*)X1, &Tplus1, (float*)X1, &Tplus1, &f_beta, (float*)XX, &n1);
      sgemm("T", "N", &n1, &N, &T, &f_alpha, (float*)X1, &Tplus1, (float*)Y1, &T, &f_beta, (float*)XY, &n1);
    } else if(typeid(realtype) == typeid(doubletype)){
      double d_alpha = 1.0;
      double d_beta = 0.0;
      dgemm("T", "N", &n1, &n1, &T, &d_alpha, (double*)X1, &Tplus1, (double*)X1, &Tplus1, &d_beta, (double*)XX, &n1);
      dgemm("T", "N", &n1, &N, &T, &d_alpha, (double*)X1, &Tplus1, (double*)Y1, &T, &d_beta, (double*)XY, &n1);
    }

    // Add regularization matrix
    REAL eta = (T/n1)*pow(10, penalty);
    for(ix = 0 ; ix < n1 ; ++ix) XX[ix + ix*n1] += eta;

    // Solve the system of linear equations
    int ipiv[n1];
    int info;
    if(typeid(realtype) == typeid(singletype)){
      sgesv(&n1, &N, (float*)XX, &n1, ipiv, (float*)XY, &n1, &info);
    } else if(typeid(realtype) == typeid(doubletype)){
      dgesv(&n1, &N, (double*)XX, &n1, ipiv, (double*)XY, &n1, &info);
    }
    for(ix = 0 ; ix < n1 ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	B1[ix+jx*n1] = XY[ix+jx*n1];
      }
    }

  }

  //==========================================================================
  // 3.6 RLS-TSVD
  //==========================================================================

  else if(RM == 6){

    // If the regression method is RLS-TSVD, ...

    // storage for linear regression variables
    REAL S[n1];
    float f_rcond = 1/pow(10, penalty);
    double d_rcond = 1/pow(10, penalty);
    int rank;
    float f_wkopt;
    double d_wkopt;
    REAL* work = NULL;
    int info;
    int lwork;
    
    // Compute B using the formula of the OLS estimator; note that
    // all the regressions, j=1,...,N, are run at once
    lwork = -1;
    if(typeid(realtype) == typeid(singletype)){
      sgelss(&T, &n1, &N, (float*)X1, &Tplus1, (float*)Y1, &T, (float*)S, &f_rcond, &rank, &f_wkopt, &lwork, &info);
      lwork = (int)f_wkopt;
    } else if(typeid(realtype) == typeid(doubletype)){
      dgelss(&T, &n1, &N, (double*)X1, &Tplus1, (double*)Y1, &T, (double*)S, &d_rcond, &rank, &d_wkopt, &lwork, &info);
      lwork = (int)d_wkopt;
    }
    work = (REAL*)realloc(work, lwork*sizeof(REAL));
    if(typeid(realtype) == typeid(singletype)){
      sgelss(&T, &n1, &N, (float*)X1, &Tplus1, (float*)Y1, &T, (float*)S, &f_rcond, &rank, (float*)work, &lwork, &info);
    } else if(typeid(realtype) == typeid(doubletype)){
      dgelss(&T, &n1, &N, (double*)X1, &Tplus1, (double*)Y1, &T, (double*)S, &d_rcond, &rank, (double*)work, &lwork, &info);
    }
    for(ix = 0 ; ix < n1 ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	B1[ix+jx*n1] = Y1[ix+jx*T];
      }
    }

  }

  /*
// 3.7 RLAD-PP
//------------
elseif RM == 7;          // If the regression method is RLAD-PP, ...  
// We solve the linear programming problem (34)-(37) in JMM (2011). 
// xlp=[phi_plus; phi_minus; ups_plus; ups_minus] where phi_plus, phi_minus, ups_plus, ups_minus are defined in 
// (34)-(37) of JMM (2011) 
   LB = [zeros(2*n1,1); zeros(2*T,1)];
                         // Lower bounds on phi_plus, phi_minus, ups_plus, ups_minus are 0
   UB = [];              // No upper bounds at all
   f = [10^penalty*ones(n1*2,1)*T/n1; ones(2*T,1)];
                         // Regularization parameter is T/n1*10^-penalty; 
                         // (2*n1+2T)-by-1
   Aeq =  [X1 -X1 eye(T,T) -eye(T,T)];
                         // T-by-(2*n1+2*T)
   B = zeros(size(X1,2),N);
                         // Allocate memory for the regression coefficients;
                         // n1-by-N
   for j = 1:N           // For each regression j, ...
       beq =  Y1(:,j);   // T-by-1
       [xlp,fval,exitflag,output,lambda] = linprog(f,[],[],Aeq,beq,LB,UB,[]);
                         // Find xlp in which xlp(1:n1,1) corresponds to phi_plus 
                         // and xlp(n1+1:2*n1,1) corresponds to phi_minus
       B(:,j) = xlp(1:n1,1)-xlp(n1+1:2*n1,1);
                         // Coefficients for the regression j, B(:,j),
                         // are given by the difference phi_plus - phi_minus;   
                         // store these coefficients into the matrix B 
   end

// 3.8 RLAD-DP
//------------
elseif RM == 8;          // If the regression method is RLAD-DP, ...
// We solve the linear programming problem (38)-(41) in JMM (2011). 
// xlp=[q] where q is a vector of unknowns in (38)-(41) of JMM (2011) 
   LB = -ones(1,T);      // Lower bound on q is -1
   UB = ones(1,T);       // Upper bound on q is 1
   Aineq = [X1'; -X1'];    // 2*n1-by-T
   bineq = 10^penalty*ones(n1*2,1)*T/n1;
                         // *T/n1*10^-penalty is a regularization parameter;
                         // 2*n1-by-1
   B = zeros(size(X1,2),N);
                         // Allocate memory for the regression coefficients;
                         // n1-by-N
   for j = 1:N           // For each regression j, ...
       f = -Y1(:,j);     // Our objective function is -Y1'(:,j)*q (minus sign 
                         // appears because (38)-(41) in JMM (2011) is a 
                         // minimization problem)
       [xlp,fval,exitflag,output,lambda] = linprog(f,Aineq,bineq,[],[],LB,UB,[]);
                         // Find the values of the Lagrange multipliers, 
                         // lambda, on all the constraints; phi_plus and phi_minus 
                         // (defined in (38)-(41) in JMM, 2011) are equal to  
                         // the Lagrange multipliers lambda.ineqlin(1:n1) 
                         // and lambda.ineqlin(n1+1:2*n1), respectively
       B(:,j) = lambda.ineqlin(1:n1)-lambda.ineqlin(n1+1:2*n1);
                         // Coefficients for the regression j, B(:,j),
                         // are given by the difference phi_plus - phi_minus
   end
   
end
  */

// 10. Infer the regression coefficients in the original regression with
// unnormalized data
//----------------------------------------------------------------------

  if((normalize == 1) | (RM >= 5)){

    // If data were normalized, ... 
 
    // Infer all the regression coefficients except the intercept
    for(ix = 1 ; ix < n ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	B[ix + jx*n] = (Y_sds[jx]/X_sds[ix-1])*B1[(ix-1)+jx*n1];
      }
    }

    // Infer the intercept
    for(jx = 0 ; jx < N ; ++jx){
      B[jx*n] = Y_means[jx];
      for(ix = 1 ; ix < n ; ++ix){
	B[jx*n] -= X_means[ix-1]*B[ix+jx*n];
      }
    }
  }

  free(X1);
  free(Y1);
  free(B1);
  free(X_means);
  free(X_sds);
  free(Y_means);
  //free(Y_sds);
  return B;
}
