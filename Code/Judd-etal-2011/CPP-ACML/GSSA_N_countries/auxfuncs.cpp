#include "global.h"
#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;

REAL mean(REAL* x, const int N)
{
  int ix;
  REAL mean = 0.0;
  for(ix = 0 ; ix < N ; ++ix) mean += x[ix]/N;
  return mean;
}

REAL var(REAL* x, const int N)
{
  int ix;
  const REAL x_bar = mean(x,N);
  REAL var = 0.0;
  for(ix = 0 ; ix < N ; ++ix) var += pow(x[ix]-x_bar,2)/(N-1);
  return var;
}

REAL sd(REAL* x, const int N)
{
  return sqrt(var(x,N));
}

REAL norm(REAL* x, const int N)
{
  int ix;
  REAL norm = 0.0;
  for(ix = 0 ; ix < N ; ++ix) norm += fabs(x[ix]);
  return norm;
}

void matrix_normalize(const int M, const int N, REAL* X, const int ldx,
		      REAL* X_norm, REAL* X_means, REAL* X_sds)
{

  // M-by-N matrix stored in flat array, row-major format
  // normalize value with column mean and standard deviation
  int ix, jx;
  for(jx = 0 ; jx < N ; ++jx){
    X_means[jx] = mean(X+jx*ldx, M);
    X_sds[jx] = sd(X+jx*ldx, M);
    for(ix = 0 ; ix < M ; ++ix){
      X_norm[ix + jx*ldx] = (X[ix + jx*ldx] - X_means[jx])/X_sds[jx];
    }
  }
}


void print_matrix(const int M, const int N, REAL* X,
		  const int printrows, const int printcols)
{
  cout.precision(5);
  for(int ix = 0 ; ix < printrows ; ++ix){
    for(int jx = 0 ; jx < printcols ; ++jx){
      cout << setw(10) << X[ix+jx*M] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

void print_vector(const int N, REAL* X)
{
  cout.precision(5);
  for(int ix = 0 ; ix < N ; ++ix){
    cout << X[ix] << endl;
  }
}
