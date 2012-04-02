/* C++ software that solves the multi-country model using the generalized 
   stochastic-simulation algorithm (GSSA), as described in the article  
   "Numerically Stable and Accurate Stochastic Simulation Approaches for 
   Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia Maliar and 
   Serguei Maliar, (2011), Quantitative Economics 2/2, 173-210 (henceforth, 
   JMM, 2011). This is a translation of the original Matlab code provided
   by Lilia Maliar and Serguei Maliar.
 
   This version: 28 March 2012.
 
  -----------------------------------------------------------------------------
  The software uses the following files: 
  -----------------------------------------------------------------------------
  1. "Main_GSSA_N.cpp"      computes GSSA solutions to the N-country model
  2. "Accuracy_Test_N.cpp"  computes approximation errors in the optimality 
                            conditions on a given set of points in the state 
                            space, for the N-country model 
  3. "Productivity.cpp"     generates random draws of the productivity shocks  
                            and simulates the corresponding series of the  
                            productivity levels
  4. "Num_Stab_Approx.cpp"  implements the numerically stable LS and LAD 
                            approximation methods
  5. "Ord_Polynomial_N.cpp" constructs the sets of basis functions for ordinary
                            polynomials of the degrees from one to five, for
                            the N-country model
  6. "Monomials_1.cpp"      constructs integration nodes and weights for an N-
                            dimensional monomial (non-product) integration rule 
                            with 2N nodes 
  7. "Monomials_2.cpp"      constructs integration nodes and weights for an N-
                            dimensional monomial (non-product) integration rule 
                            with 2N^2+1 nodes
  8. "GH_Quadrature.cpp"    constructs integration nodes and weights for the 
                            Gauss-Hermite rules with the number of nodes in 
                            each dimension ranging from one to ten                     
  9. "aT20200N10.dat"       contains the series of the productivity levels of 
                            length 20,200 for 10 countries that are used for 
                            computing solutions and for evaluating accuracy 
  10. "ran.h"               generates a normal random variable
  11. "cholesky.cpp"        computes a Cholesky decomposition of a matrix
  12. "qrdcmp.cpp"          computes the QR decomposition of square matrix and
                            provides a function to solve a linear system with
                            the resulting decomposition
  13. "simplx.cpp"          implements the simplex method to solve a linear
                            programming problem in standard form
  14. "svdcmp.cpp"          computes the singular value decomposition of a
                            rectangular matrix and provides a function to
                            solve a linear system with the resulting
                            decomposition 
  15. "auxfuncs.cpp"        various helper functions
  -----------------------------------------------------------------------------
  Copyright © 2012 by Eric M. Aldrich. All rights reserved. The code may be
  used, modified and redistributed under the terms provided in the file
  "License_Agreement.txt".
  ---------------------------------------------------------------------------*/

#include "global.h"
#include "svdcmp.h"
#include "qrdcmp.h"
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <typeinfo>

using namespace std;

int main()
{

  int hx, ix, jx, kx, tx;

  //==========================================================================
  // 1. Choose the number of countries and simulation length
  //==========================================================================
  int N     = 2;    // Choose the number of countries, 1<=N<=10 (note that the
                    // code also works for the one-country case, N=1)
  int Ndat;         // Either number of countries in existing data set or
                    // number of countries in created data set
  int T     = 200; // Choose the simulation length for the solution procedure,
                    // T<=10,000   

  // To solve models with N>10 or T>10,000, one needs to simulate new series
  // of the productivity levels by enabling the code in paragraph 5.

  //==========================================================================
  // 2. Model's parameters
  //==========================================================================
  REAL gam     = 1;      // Utility-function parameter
  REAL alpha   = 0.36;   // Capital share in output
  REAL beta    = 0.99;   // Discount factor
  REAL delta   = 0.025;  // Depreciation rate 
  REAL rho     = 0.95;   // Persistence of the log of the productivity level
  REAL sigma   = 0.01;   // Standard deviation of shocks to the log of the 
                         // productivity level

  // Variance-covariance matrix of the countries' productivity 
  // shocks in which diagonal terms are equal to 2*sigma^2   
  // and in which off-diagonal terms are equal to sigma^2; 
  // this vcv follows from the assumption that a country's 
  // shock has both common-for-all-countries and country-
  // specific components; N-by-N, column-major format
  REAL* vcv = new REAL[N*N];
  for(jx = 0 ; jx < N ; ++jx){
    for(ix = 0 ; ix < N ; ++ix){
      if(jx == ix){
	vcv[ix + jx*N] = 2*pow(sigma,2);
      } else {
	vcv[ix + jx*N] = pow(sigma,2);
      }
    }
  }

                                     
  //==========================================================================
  // 3. The normalizing constant, A, and welfare weight, tau
  //==========================================================================

  // The normalizing constant in output  
  REAL A       = (1-beta+beta*delta)/alpha/beta;

  // The welfare weight of country j 
  REAL tau     = 1;

  // The above normalization ensures that steady state of capital of all 
  // countries is equal to one 

  //==========================================================================
  // 4. Initial condition
  //==========================================================================

  // Initial conditions for capital and productivity (equal to steady state)
  REAL* k = new REAL[(T+1)*N];
  REAL* a = new REAL[T*N];
  for(jx = 0 ; jx < N ; ++jx){
    k[jx*(T+1)] = 1;
    a[jx*T] = 1;
  }

  //==========================================================================
  // 5. Construct the productivity levels, a 
  //==========================================================================

  // simulate or load existing data
  bool simulate = 0;
  const char* fileName = "aT20200N10.dat";

  if(simulate){

    // Generate a random draw of the productivity shocks and simulate the
    // corresponding series of the productivity levels of length T periods
    // for N countries 
    Productivity(T, N, a, sigma, rho);

    // Save the series of the productivity levels into a file "aT20200N10.dat" 
    // Save in row-major format (althought variable is in column-major format)
    ofstream fileOut;
    fileOut.open(fileName);
    for(ix = 0 ; ix < T ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	fileOut << a[ix + jx*T] << "\n";
      }
    }
    fileOut.close();

    // number of countries in data set
    Ndat = N;

  } else {
    
    // Load the previously saved series of the productivity levels of length
    // 20,200 for 10 countries (the first 10,000 observations are used for
    // finding a solution, and the remaining 10,200 observations are used for
    // evaluating accuracy)
    // Restrict the series of the productivity levels for the solution
    // procedure to the given T<=10,000 and N<=10
    Ndat = 10; // the number of countries in the exisiting dataset
    ifstream fileIn;
    REAL trash;
    fileIn.open(fileName);
    for(ix = 0 ; ix < T ; ++ix){
      for(jx = 0 ; jx < Ndat ; ++jx){
	if(jx < N){
	  fileIn >> a[ix + jx+T];
	} else {
	  fileIn >> trash;
	}
      }
    }
    fileIn.close();
  }

  //==========================================================================
  // Compute a first-degree polynomial solution using the one-node Monte Carlo  
  // integration method (this solution will be used as an initial guess for
  // the other cases)
  //==========================================================================
  clock_t tic = clock(); // Start counting time needed to compute the solution
                            
  //==========================================================================
  // 6. The GSSA parameters
  //==========================================================================

  // Damping parameter for (fixed-point) iteration on the coefficients of
  // the capital policy functions
  REAL kdamp = 0.1; 

  // Set the initial difference between the series from two iterations in
  // the convergence criterion (condition (10) in JMM, 2011) to a very large
  // number
  // To achieve convergence under N>10, one may need to modify the values of 
  // the damping parameter kdamp or refine the initial guess 

  REAL dif_1d = 1e+10;

  //==========================================================================
  // 7. Initialize the first-degree capital policy functions of N countries 
  //==========================================================================

  REAL* bk_1d = new REAL[(1+2*N)*N];
  for(ix = 0 ; ix < (1+2*N) ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      if(ix == (jx+1)){
	bk_1d[ix + jx*(1+2*N)] = 0.9;
      } else if(ix == (jx+1+N)){
	bk_1d[ix + jx*(1+2*N)] = 0.1;
      } else {
	bk_1d[ix + jx*(1+2*N)] = 0;
      }
    }
  }
  // Matrix of polynomial coefficients of size (1+2N)-by-N: for each country  
  // (each column), 1+2N rows correspond to a constant, N coefficients on the
  // capital stocks, k(t,1),...,k(t,N), and N coefficients on the productivity
  // levels, a(t,1),...,a(t,N)

  // As an initial guess, assume that a country's j capital depends only on 
  // its own capital and productivity level as k(t+1,j)=0.9*k(t,j)+0.1*a(t,j); 
  // (note that in the steady state, we have k(t+1,j)=0.9*k(t,j)+0.1*a(t,j)=1)
  
  // Note that diag(ones(1,N)*q) delivers an N-by-N matrix with  diagonal
  // entries equal to q.
    
  //==========================================================================
  // 8. Initialize the capital series
  //==========================================================================

  // Initialize the series of next-period capital of N countries; these
  // series are used to check the convergence on the subsequent iteration
  // (initially, capital can take any value); (T+1)-by-N
  REAL* k_old = new REAL[(T+1)*N];
  for(ix = 0 ; ix < (T+1) ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      k_old[ix + jx*(T+1)] = 1.0;
    }
  }

  //==========================================================================
  // 9. The main iterative cycle of GSSA
  //==========================================================================

  // auxiliary variables and storage
  REAL* k_temp;                    // damping parameter
  REAL* x = new REAL[T*(1+2*N)];   // storage for basis functions
  REAL* C = new REAL[T];           // storage for aggregate consumption
  REAL* c = new REAL[T*N];         // storage for individual consumption
  REAL* Y = new REAL[(T-1)*N];     // storage for income
  REAL a1[T*N], k2[T*N], c1[T*N];  // variables for euler equation errors
  
  // storage for linear regression variables (including SVD variables)
  int count = 0;
  REAL wkopt;
  REAL* work = NULL;
  int info;
  int lwork = -1;

  // loop
  while(dif_1d > 1e-4*kdamp){

    //========================================================================
    // 9.1 Generate time series of capital
    //========================================================================

    for(ix = 0 ; ix < T ; ++ix){
      for(jx = 0 ; jx < (1+2*N) ; ++jx){

	// The basis functions of the first-degree polynomial at time t
	if(jx == 0){
	  x[ix + jx*T] = 1.0;
	} else if((jx >= 1) & (jx <= N)){
	  x[ix + jx*T] = k[ix + (jx-1)*T];
	} else {
	  x[ix + jx*T] = a[ix + (jx-N-1)*T];
	}

	// Compute next-period capital using bk_1d
	REAL alph = 1.0;
	REAL bet = 0.0;
	int nrows = 1+2*N;
	int inc = 1;
	char trans = 'T';
	//if(typeid(realtype) == typeid(singletype)){
	  //SGEMV("T", &nrows, &N, &(float)alph, (float*)bk_1d, &nrows, (float*)x+ix, &inc,
	  //	&(float)bet, (float*)k+(ix+1), &inc);
	//} else if(typeid(realtype) == typeid(doubletype)){
	  //dgemv("T", &(1+2*N), &N, &1.0, (double*)bk_1d, &(1+2*N), (double*)x+ix, &1,
	  //	&0.0, (double*)k+(ix+1), &1);
	  dgemv(&trans, &nrows, &N, &alph, bk_1d, &nrows, x+ix, &inc, &bet, k+(ix+1), &inc);
	  //}
      }
    }

    //========================================================================
    // 9.2 Compute consumption series 
    //========================================================================

    for(ix = 0 ; ix < T ; ++ix){
      C[ix] = 0.0;
      for(jx = 0 ; jx < N ; ++jx){
	C[ix] += A*pow(k[ix+jx*T], alpha)*a[ix+jx*T] - k[(ix+1)+jx*T] + k[ix+jx*T]*(1-delta);
      }
    }

    // Aggregate consumption is computed by summing up individual consumption, 
    // which in turn, is found from the individual budget constraints; T-by-1

    // Individual consumption is the same for all countries; T-by-N
    for(ix = 0 ; ix < T ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	c[ix + jx*T] = C[ix]/N;
      }
    }

    //========================================================================
    // 9.3 Evaluate the percentage (unit-free) difference between the series  
    //========================================================================

    // from the previous and current iterations
    dif_1d = 0.0;
    for(ix = 0 ; ix < (T+1) ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	dif_1d += fabs(1 - k[ix + jx*(T+1)]/k_old[ix + jx*(T+1)]);
      }
    }
    dif_1d /= (T+1)*N;
    // Compute a unit-free difference between the capital series 
    // from two iterations; see condition (10) in JMM (2011)

    //========================================================================
    // 9.4 Monte Carlo realizations of the right side of the Euler equation,
    // Y, in condition (C4) in the online Appendix C
    //========================================================================

    for(ix = 0 ; ix < (T-1) ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	Y[ix+jx*(T-1)] = beta*pow(c[(ix+1)+jx*(T-1)]/c[ix+jx*(T-1)], -gam)*(1-delta+A*alpha*pow(k[(ix+1)+jx*(T-1)], alpha-1)*a[(ix+1)+jx*(T-1)])*k[(ix+1)+jx*(T-1)];
      }
    }

    //========================================================================
    // 9.5 Compute and update the coefficients of the capital policy functions
    //========================================================================

    // Compute new coefficients of capital policy functions using OLS (QR).
    if(typeid(realtype) == typeid(singletype)){
      sgels('N', &(T-1), &(1+2*N), &N, (float*)x, &(T-1), (float*)Y, &(T-1), &wkopt, &lwork, &info);
    } else if(typeid(realtype) == typeid(doubletype)){
      dgels('N', &(T-1), &(1+2*N), &N, (double*)x, &(T-1), (double*)Y, &(T-1), &wkopt, &lwork, &info);
    }
    lwork = (int)wkopt;
    work = (REAL*)realloc(work, lwork*sizeof(REAL));
    if(typeid(realtype) == typeid(singletype)){
      sgels('N', &(T-1), &(1+2*N), &N, (float*)x, &(T-1), (float*)Y, &(T-1), work, &lwork, &info);
    } else if(typeid(realtype) == typeid(doubletype)){
      dgels('N', &(T-1), &(1+2*N), &N, (double*)x, &(T-1), (double*)Y, &(T-1), work, &lwork, &info);
    }

    // Update the coefficients of the capital policy functions using damping 
    for(ix = 0 ; ix < (1+2*N) ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	bk_1d[ix+jx*(1+2*N)] = kdamp*Y[ix+jx*(1+2*N)] + (1-kdamp)*bk_1d[ix+jx*(1+2*N)];
      }
    }

                                     
    //========================================================================
    // 9.6 Store the capital series 
    //========================================================================

    // The stored capital series will be used for checking the convergence
    // on the subsequent iteration
    k_temp = k_old;
    k_old = k;
    k = k_temp;
    ++count;
    cout<< "Iteration: " << count <<endl;
    cout<< "dif_1d: " << dif_1d <<endl;
  }

  //==========================================================================
  // 10. Time needed to compute the initial guess
  //==========================================================================
  clock_t toc = clock();
  REAL time_GSSA_1d = (toc - tic)/(REAL)CLOCKS_PER_SEC; 
  cout << time_GSSA_1d << endl;


  //==========================================================================
  // Compute polynomial solutions of the degrees from one to D_max using one 
  // of the following integration methods: Monte Carlo, Gauss-Hermite product  
  // and monomial non-product methods
  //==========================================================================
  tic = clock();      // Start counting time needed to compute the solution
                            
  //==========================================================================
  // 11. The GSSA parameters
  //==========================================================================

  // Damping parameter for (fixed-point) iteration on the coefficients of
  // the capital policy functions
  kdamp = 0.1;

  // Set the initial difference between the series from two iterations in
  // the convergence criterion (condition (10) in JMM, 2011) to a very
  // large number
  REAL dif_GSSA_D  = 1e+10;

  // To achieve convergence under N>10, one may need to modify the values of 
  // the damping parameter kdamp or refine the initial guess 

  
  //==========================================================================
  // 12. The matrix of the polynomial coefficients
  //==========================================================================
  // Maximum degree of a polynomial: the program computes polynomial
  // solutions of the degrees from one to D_max; (D_max can be from 1 to 5)
  int D_max  = 5;

  // For the polynomial degrees from one to D_max compute the number of
  // polynomial bases (this is needed for finding the number of the
  // polynomial coefficients in the policy functions)
  int* npol = new int[D_max];
  npol[0] = 1 + 2*N;
  int n_high_order = 2*N;
  for(ix = 1 ; ix < D_max ; ++ix){
    n_high_order = (n_high_order*(2*N+ix))/(ix+1);
    npol[ix] = npol[ix-1] + n_high_order;
  }

  // Matrix of polynomial coefficients of the capital policy functions for
  // the polynomial solutions of the degrees from one to D_max;
  // npol[D_max-1]-by-N-by-D_max
  REAL* BK = new REAL[npol[D_max-1]*N*D_max];
  for(ix = 0 ; ix < npol[D_max-1] ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      for(kx = 0 ; kx < D_max ; ++kx){
	BK[ix*N*D_max + jx*D_max + kx] = 0.0;
      }
    }
  }

  //==========================================================================
  // 13. Choose an integration method 
  //==========================================================================
                                 
  // 0 = a one-node Monte Carlo method(default);
  // 1,2,..,10 = Gauss-Hermite quadrature rules with 1,2,...,10 nodes in each
  // dimension, respectively;
  // 11 = Monomial rule with 2N nodes;
  // 12 = Monomial rule with 2N^2+1 nodes
  int IM = 12;
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

  // Under the one-node Gauss-Hermite quadrature rule, the conditional 
  // expectation (integral) is approximated by the value of the integrand 
  // evaluated in one integration node in which the next-period productivity 
  // shock is zero, i.e., the next-period productivity level is 
  // a(t+1,:)=a(t,:).^rho*exp(0)=a(t,:).^rho

  //==========================================================================
  // 14. Choose an regression method 
  //==========================================================================

  int RM = 6;        // Choose a regression method: 
                     // 1=OLS,          2=LS-SVD,   3=LAD-PP,  4=LAD-DP, 
                     // 5=RLS-Tikhonov, 6=RLS-TSVD, 7=RLAD-PP, 8=RLAD-DP
  int normalize = 1; // Option of normalizing the data; 0=unnormalized data; 
                     // 1=normalized data                    
  int penalty = 7;   // Degree of regularization for a regularization methods, 
                     // RM=5,6,7,8 (must be negative, e.g., -7 for RM=5,7,8 
                     // and must be positive, e.g., 7, for RM=6)

  //==========================================================================
  // 15. Compute the polynomial solutions of the degrees from one to D_max
  //==========================================================================

  // Admin
  int n_cols, D;
  REAL* poly_X;
  REAL* poly_X_integral;
  REAL* poly_X_row;
  REAL* bk_D;
  REAL* bk_hat_D;
  REAL time_GSSA[D_max];

  // Matrix which holds data for polynomial bases - begin with only shocks
  REAL X1[T*2*N], X1_integral[T*2*N];
  for(ix = 0 ; ix < T ; ++ix){
    for(jx = N ; jx < 2*N ; ++jx){
      X1[ix*2*N + jx] = a[ix*N + jx - N];
    }
  }

  for(D = 1 ; D <= D_max ; ++D){

    //========================================================================
    // 15.1 Using the previously computed capital series, compute the initial 
    // guess on the coefficients under the  selected approximation method
    //========================================================================

    // Complete the underlying data matrix, filling with most recent capital
    for(ix = 0 ; ix < T ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	X1[ix*2*N + jx] = k[ix*N + jx];
      }
    }

    // Construct the polynomial bases on the series of state variables from the
    // previously computed time-series solution
    // eliminate the first column of X
    poly_X = Ord_Polynomial_N(X1, T, 2*N, D, n_cols);
   
    // Compute the initial guess on the coefficients  using the chosen
    // regression method
    bk_D = Num_Stab_Approx(T-1, n_cols, poly_X, N, Y, RM, penalty, normalize);

    // Initialize the series of next-period capital of N countries; these
    // series are used to check the convergence on the subsequent iteration
    // (initially, capital can take any value); (T+1)-by-N
    for(ix = 0 ; ix < (T+1) ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	k_old[ix*N + jx] = 1.0;
      }
    }

    // Convergence criterion (initially is not satisfied)
    dif_GSSA_D  = 1e+10;

    //========================================================================
    // 15.2 The main iterative cycle of GSSA
    //========================================================================
    
    // 10^(-4-D)*kdamp is a convergence parameter, adjusted to the polynomial
    // degree D and the damping parameter kdamp; see the discussion in
    // JMM (2011)
    count = 0;
    while(dif_GSSA_D > pow(10,-4-D)*kdamp){
    
      //======================================================================
      // 15.2.1 Generate time series of capital
      //======================================================================

      for(ix = 0 ; ix < T ; ++ix){

	// The basis functions of a polynomial of degree D at time t	
	for(jx = 0 ; jx < 2*N ; ++jx){
	  if(jx < N){
	    X1[ix*2*N + jx] = k[ix*N + jx];
	  } else {
	    X1[ix*2*N + jx] = a[ix*N + jx - N];
	  }
	}

	// Polynomial bases
	poly_X_row = Ord_Polynomial_N(X1+ix*2*N, 1, 2*N, D, n_cols);
	for(jx = 0 ; jx < n_cols ; ++jx){
	  poly_X[ix*n_cols+jx] = poly_X_row[jx];
	}

	// Compute next-period capital using bk_D	
	if(typeid(realtype) == typeid(singletype)){
	  cblas_sgemv(CblasRowMajor, CblasTrans, n_cols, N, 1.0, (float*)bk_D,
		      N, (float*)poly_X_row, 1, 0.0, (float*)k+(ix+1)*N, 1);
	} else if(typeid(realtype) == typeid(doubletype)){
	  cblas_dgemv(CblasRowMajor, CblasTrans, n_cols, N, 1.0, (double*)bk_D,
		      N, (double*)poly_X_row, 1, 0.0, (double*)k+(ix+1)*N, 1);
	}

      }

      //======================================================================
      // 15.2.2 Compute consumption series of all countries 
      //======================================================================

      // Aggregate consumption is computed by summing up individual consumption, 
      // which in turn, is found from the individual budget constraints; T-by-1
      for(ix = 0 ; ix < T ; ++ix){
	C[ix] = 0.0;
	for(jx = 0 ; jx < N ; ++jx){
	  C[ix] += A*pow(k[ix*N+jx], alpha)*a[ix*N+jx] - k[(ix+1)*N+jx] + k[ix*N+jx]*(1-delta);
	}
      }

      // Individual consumption is the same for all countries; T-by-N
      for(ix = 0 ; ix < T ; ++ix){
	for(jx = 0 ; jx < N ; ++jx){
	  c[ix*N + jx] = C[ix]/N;
	}
      }
      
      //======================================================================
      // 15.2.3 Approximate the conditional expectations for t=1,...T-1 using 
      // the integration method chosen 
      //======================================================================

      //======================================================================
      // 15.2.3.1 The one-node Monte Carlo integration method approximates the
      // values of the conditional expectations, Y, in the Euler equation with 
      // the realization of the integrand in the next period 
      //======================================================================
      if(IM == 0){
	for(ix = 0 ; ix < (T-1) ; ++ix){
	  for(jx = 0 ; jx < N ; ++jx){
	    Y[ix*N+jx] = beta*pow(c[(ix+1)*N+jx]/c[ix*N+jx], -gam)*(1-delta+A*alpha*pow(k[(ix+1)*N+jx], alpha-1)*a[(ix+1)*N+jx])*k[(ix+1)*N+jx];
	  }
	}
      }

      //======================================================================
      // 15.2.3.2 Deterministic integration methods approximate the values of 
      // conditional expectations, Y, in the Euler equation as a weighted
      // average of the values of the integrand in the given nodes with the
      // given weights 
      //======================================================================
      else{

	// Initialize variable Y
	for(ix = 0 ; ix < T ; ++ix){
	  for(jx = 0 ; jx < N ; ++jx){
	    Y[ix*N+jx] = 0.0;
	  }
	}

	// Compute Euler Equation Errors via integration
	for(hx = 0 ; hx < n_nodes ; ++hx){
	  for(ix = 0 ; ix < T ; ++ix){
	    for(jx = 0 ; jx < N ; ++jx){

	      // Compute the next-period productivity levels for each integration       
	      // node using condition (C3) in the online Appendix C; n_nodes-by-N
	      a1[ix*N+jx] = pow(a[ix*N+jx],rho)*exp(epsi_nodes[hx*N+jx]);

	      // The basis functions of a polynomial of degree D at time t	
	      X1_integral[ix*2*N + jx] = k[(ix+1)*N + jx];
	      X1_integral[ix*2*N + jx + N] = a1[ix*N + jx];
	    }
	  }
	  
	  // Polynomial bases
	  poly_X_integral = Ord_Polynomial_N(X1_integral, T, 2*N, D, n_cols);

	  // Compute capital of period t+2 (chosen at t+1) using the
	  // capital policy functions; n_nodes-by-N 
	  if(typeid(realtype) == typeid(singletype)){
	    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, T, N, n_cols,
			1.0, (float*)poly_X_integral, n_cols, (float*)bk_D, N, 0.0,
			(float*)k2, N);
	  } else if(typeid(realtype) == typeid(doubletype)){
	    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, T, N, n_cols,
			1.0, (double*)poly_X_integral, n_cols, (double*)bk_D, N, 0.0,
			(double*)k2, N);
	  }

	  // C is computed by summing up individual consumption, which in
	  // turn, is found from the individual budget constraints; T-by-1
	  for(ix = 0 ; ix < T ; ++ix){
	    C[ix] = 0.0;
	    for(jx = 0 ; jx < N ; ++jx){
	      C[ix] += A*pow(k[(ix+1)*N+jx], alpha)*a1[ix*N+jx] - k2[ix*N+jx] + k[(ix+1)*N+jx]*(1-delta);
	    }
	  }

	  // Compute next-period consumption for N countries; n_nodes-by-N
	  for(ix = 0 ; ix < T ; ++ix){
	    for(jx = 0 ; jx < N ; ++jx){
	      c1[ix*N + jx] = C[ix]/N;
	    }
	  }

	  // Sum Euler Equation errors over integration weights
	  for(ix = 0 ; ix < T ; ++ix){
	    for(jx = 0 ; jx < N ; ++jx){
	      Y[ix*N+jx] = Y[ix*N+jx]+weight_nodes[hx]*beta*pow(c1[ix*N+jx]/c[ix*N+jx],-gam)*(1-delta+alpha*A*pow(k[(ix+1)*N+jx],alpha-1)*a1[ix*N+jx])*k[(ix+1)*N+jx];
	    } // T-by-N
	  }
	}
      }
 
      //======================================================================
      // 15.2.4 Evaluate the percentage (unit-free) difference between the 
      // capital series from the previous and current iterations
      //======================================================================

      // Compute a unit-free difference between the capital series from two
      // iterations; see condition (10) in JMM (2011)
      dif_GSSA_D = 0.0;
      for(ix = 0 ; ix < (T+1) ; ++ix){
	for(jx = 0 ; jx < N ; ++jx){
	  dif_GSSA_D += fabs(1 - k[ix*N + jx]/k_old[ix*N + jx]);
	}
      }
      dif_GSSA_D /= (T+1)*N;

      //======================================================================
      // 15.2.5 Compute and update the coefficients of the capital policy 
      // functions
      //======================================================================

      // Compute new coefficients of the capital policy functions using the
      // chosen approximation method
      bk_hat_D = Num_Stab_Approx(T-1, n_cols, poly_X, N, Y, RM, penalty, normalize);

      // Update the coefficients of the capital policy functions using damping 
      for(ix = 0 ; ix < n_cols ; ++ix){
	for(jx = 0 ; jx < N ; ++jx){
	  bk_D[ix*N+jx] = kdamp*bk_hat_D[ix*N+jx] + (1-kdamp)*bk_D[ix*N+jx];
	}
      }

      //======================================================================
      // 15.2.6 Store the capital series 
      //======================================================================

      // The stored capital series will be used for checking the convergence
      // on the subsequent iteration
      k_temp = k_old;
      k_old = k;
      k = k_temp;
      ++count;
      cout << "Order of Integration: " << D << endl;
      cout<< "Iteration: " << count <<endl;
      cout<< "dif_GSSA_D: " << dif_GSSA_D <<endl;
      
    }

    //========================================================================
    // 15.2.7 The GSSA output for the polynomial solution of degree D
    //========================================================================
    // Store the coefficients of the polynomial of degree D that approximates
    // capital policy functions of N countries 
    for(ix = 0 ; ix < npol[D-1] ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	BK[ix*N + jx + (D-1)*N*npol[D-1]] = bk_D[ix*N+jx];
      }
    }

    // Time needed to compute the polynomial solution of degree D
    toc = clock();
    time_GSSA[D-1] = (toc - tic)/(REAL)CLOCKS_PER_SEC; 
    cout << time_GSSA[D-1] << endl;

  }

  //===========================================================================
  // 16. Accuracy test of the GSSA solutions: errors on a stochastic simulation
  //===========================================================================

  //===========================================================================
  // 16.1 Specify a set of points on which the accuracy is evaluated
  //===========================================================================

  //Choose the simulation length for the test on a stochastic simulation, 
  // T_test<=10,200
  int T_test = 10200;

  // Restrict the series of the productivity levels for the test on a
  // stochastic simulation to the given T_test<=10,200 and N<=10
  REAL* a_test = new REAL[T_test*N];
  REAL trash;
  ifstream fileIn;
  fileIn.open(fileName);
  for(ix = 0 ; ix < (T+T_test) ; ++ix){
    for(jx = 0 ; jx < Ndat ; ++jx){
      if(ix < T){
	fileIn >> a_test[0];
      } else {
	if(jx < N){
	  fileIn >> a_test[(ix-T)*N + jx];
	} else {
	  fileIn >> trash;
	}
      }
    }
  }
  fileIn.close();
  

  // Initial condition for capital (equal to steady state)
  REAL* k_test = new REAL[T_test*N];
  for(ix = 0 ; ix < T_test ; ++ix){
    for(jx = 0 ; jx < N ; ++jx){
      k_test[ix*N + jx] = 1.0;
    }
  }

  //===========================================================================
  // 16.2 Choose an integration method for evaluating accuracy of solutions
  //===========================================================================

  // See paragraph 13 for the integration options
  int IM_test = 11;

  // To implement the test on a stochastic simulation with T_test>10,200, one
  // needs to simulate new series of the productivity levels with larger T_test 
  // by enabling the code in paragraph 6.

  //===========================================================================
  // 16.3 Compute errors on a stochastic simulation for the GSSA polynomial
  // solution of degrees D=1,...,D_max
  //===========================================================================

  n_cols = npol[0];
  REAL* bk = new REAL[n_cols*N];
  REAL* X1_test = new REAL[T_test*2*N];
  REAL* poly_X_row_test = new REAL[n_cols];
  int discard = 200; // discard the first 200 observations
  REAL Errors_mean[D_max], Errors_max[D_max], time_test[D_max];
  for(D = 1 ; D <= D_max ; ++D){
    
    //=========================================================================
    // 16.3.1 Simulate the time series solution under the given capital-
    // policy-function coefficients, BK(:,:,D) with D=1,...,D_max
    //=========================================================================

    // Proper memory allocation
    if(D > 1){
      n_cols = npol[D-1];
      bk = (REAL*)realloc(bk, n_cols*N*sizeof(REAL));
      poly_X_row_test = (REAL*)realloc(poly_X_row_test, n_cols*sizeof(REAL));
    }

    // The vector of coefficients of the polynomial of degree D
    for(ix = 0 ; ix < n_cols ; ++ix){
      for(jx = 0 ; jx < N ; ++jx){
	bk[ix*N+jx] = BK[ix*N + jx + (D-1)*N*n_cols];
      }
    }

    for(ix = 0 ; ix < (T_test-1) ; ++ix){

      // The basis functions of a polynomial of degree D at time t	
      for(jx = 0 ; jx < 2*N ; ++jx){
	if(jx < N){
	  X1_test[ix*2*N + jx] = k_test[ix*N + jx];
	} else {
	  X1_test[ix*2*N + jx] = a_test[ix*N + jx - N];
	}
      }

      // Polynomial bases
      poly_X_row_test = Ord_Polynomial_N(X1_test+ix*2*N, 1, 2*N, D, n_cols);

      // Compute next-period capital using bk
      if(typeid(realtype) == typeid(singletype)){
	cblas_sgemv(CblasRowMajor, CblasTrans, n_cols, N, 1.0, (float*)bk, N,
		    (float*)poly_X_row_test, 1, 0.0, (float*)k_test+(ix+1)*N,
		    1);
      } else if(typeid(realtype) == typeid(doubletype)){
	cblas_dgemv(CblasRowMajor, CblasTrans, n_cols, N, 1.0, (double*)bk, N,
		    (double*)poly_X_row_test, 1, 0.0, (double*)k_test+(ix+1)*N,
		    1);
      }
    }

    //=========================================================================
    // 16.3.2 Errors across 10,000 points on a stochastic simulation
    //=========================================================================

    Accuracy_Test_N(T_test, N, k_test, a_test, bk, D, IM_test, alpha, gam,
		    delta, beta, A, tau, rho, vcv, discard, Errors_mean[D-1],
		    Errors_max[D-1], time_test[D-1]);

    // Errors_mean    is the unit-free average absolute approximation error  
    //                across 4N+1 optimality conditions (in log10) 
    // Errors_max     is the unit-free maximum absolute approximation error   
    //                across 4N+1 optimality conditions (in log10) 
  }

  //===========================================================================
  // 17. Display the results for the polynomial solutions of the degrees from
  // one to D_max
  //===========================================================================

  cout << "GSSA OUTPUT:\n" << endl;
  cout << "RUNNING TIME (in seconds):\n" << endl;

  cout << "a) for computing the solution\n" << endl; 
  for(D = 1 ; D <= D_max ; ++D) cout << setw(10) << D << "";
  cout << endl;
  for(D = 1 ; D <= D_max ; ++D) cout << setw(10) << time_GSSA[D-1] << "";
  cout << endl;

  cout << "\nb) for implementing the accuracy test\n" << endl;
  for(D = 1 ; D <= D_max ; ++D) cout << setw(10) << D << "";
  cout << endl;
  for(D = 1 ; D <= D_max ; ++D) cout << setw(10) << time_test[D-1] << "";
  cout << endl;

  cout << "\nAPPROXIMATION ERRORS (log10):\n" << endl;
  cout << "a) mean error across 4N+1 optimality conditions\n" << endl; 
  for(D = 1 ; D <= D_max ; ++D) cout << setw(10) << D << "";
  cout << endl;
  for(D = 1 ; D <= D_max ; ++D) cout << setw(10) << Errors_mean[D-1] << "";
  cout << endl;

  cout << "\nb) max error across 4N+1 optimality conditions\n" << endl; 
  for(D = 1 ; D <= D_max ; ++D) cout << setw(10) << D << "";
  cout << endl;
  for(D = 1 ; D <= D_max ; ++D) cout << setw(10) << Errors_max[D-1] << "";
  cout << endl;

  //save Results_N time_GSSA time_test Errors_mean Errors_max kdamp RM IM N T BK k_test a_test IM_test alpha gam delta beta A tau rho vcv discard npol D_max T_test ;

}
