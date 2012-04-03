// Ord_Polynomial_N.cpp is a routine that constructs the basis functions of 
// complete ordinary polynomial of the degrees from one to five for the 
// multi-dimensional case; see "Numerically Stable and Accurate Stochastic 
// Simulation Approaches for Solving Dynamic Economic Models" by Kenneth L. 
// Judd, Lilia Maliar and Serguei Maliar, (2011), Quantitative Economics 2/2, 
// 173-210 (henceforth, JMM, 2011). This is a translation of the original
// Matlab code provided by Lilia Maliar and Serguei Maliar.
//
// This version: 28 March 2012
// -------------------------------------------------------------------------
// Inputs:  "z" is the data points on which the polynomial basis functions  
//               must be constructed; n_rows-by-dimen; 
//          "n_rows" and "dimen" are the dimensions of matrix "z"
//          "D" is the degree of the polynomial whose basis functions must 
//               be constructed; (can be 1,2,3,4 or 5)
//          "n_cols" is the number of columns in output matrix "basis_fs"
//
// Output:  "basis_fs" is the matrix of basis functions of a complete 
//           polynomial of the given degree 
// -------------------------------------------------------------------------
// Copyright © 2012 by Eric M. Aldrich. All rights reserved. The code may 
// be used, modified and redistributed under the terms provided in the file
// "License_Agreement.pdf".
// -------------------------------------------------------------------------

#include "global.h"
#include <math.h>
#include <iostream>

using namespace std;

REAL* Ord_Polynomial_N(REAL* z, int n_rows, int dimen, int D, int& n_cols)
{

  int ix, jx, j1x, j2x, j3x, j4x, j5x;

  // A polynomial is given by the sum of polynomial basis functions, phi(i),
  // multiplied by the coefficients; see condition (13) in JMM (2011). By 
  // convention, the first basis function is one. 

  // allocate storage for the matrix of basis function
  int higher_order = 1;
  for(ix = 1 ; ix <= D ; ++ix){
    higher_order *= (dimen+ix);
    higher_order /= ix;
  }
  n_cols = higher_order;
  REAL* basis_fs = new REAL[n_rows*n_cols];

  //==========================================================================
  // The matrix of the basis functions - 1st degree
  //==========================================================================
  for(ix = 0 ; ix < n_rows ; ++ix){
    basis_fs[ix] = 1.0;
    for(jx = 1 ; jx < (dimen+1) ; ++jx){
      basis_fs[ix + jx*n_rows] = z[ix + (jx-1)*n_rows];
    }

    //========================================================================
    // 2. Second-degree columns
    //========================================================================

    if(D == 2){
      for(j1x = 0 ; j1x < dimen ; ++j1x){
	for(j2x = j1x ; j2x < dimen ; ++j2x){
	  basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows];
	  ++jx;
	}
      }

    //========================================================================
    // 3. Third-degree columns
    //========================================================================

    } else if(D == 3){
      for(j1x = 0 ; j1x < dimen ; ++j1x){
	for(j2x = j1x ; j2x < dimen ; ++j2x){
	  basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows];
	  ++jx;
	  for(j3x = j2x ; j3x < dimen ; ++j3x){
	    basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows];
	    ++jx;
	  }
	}
      }

    //========================================================================
    // 4. Fourth-degree columns
    //========================================================================

    } else if(D == 4){
      for(j1x = 0 ; j1x < dimen ; ++j1x){
	for(j2x = j1x ; j2x < dimen ; ++j2x){
	  basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows];
	  ++jx;
	  for(j3x = j2x ; j3x < dimen ; ++j3x){
	    basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows];
	    ++jx;
	    for(j4x = j3x ; j4x < dimen ; ++j4x){
	      basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows]*z[ix + j4x*n_rows];
	      ++jx;
	    }
	  }
	}
      }

    //========================================================================
    // 5. Fifth-degree columns
    //========================================================================

    } else if(D == 5){
      for(j1x = 0 ; j1x < dimen ; ++j1x){
	for(j2x = j1x ; j2x < dimen ; ++j2x){
	  basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows];
	  ++jx;
	  for(j3x = j2x ; j3x < dimen ; ++j3x){
	    basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows];
	    ++jx;
	    for(j4x = j3x ; j4x < dimen ; ++j4x){
	      basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows]*z[ix + j4x*n_rows];
	      ++jx;
	      for(j5x = j4x ; j5x < dimen ; ++j5x){
		basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows]*z[ix + j4x*n_rows]*z[ix + j5x*n_rows];
		++jx;
	      }
	    }
	  }
	}
      }
    }
  }
  return basis_fs;
}
