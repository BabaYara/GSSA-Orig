#include "global.h"
#include <math.h>

void cholesky(REAL* A, int N, REAL* A_chol)
{
  int ix, jx, kx;
  REAL sum;
  for(ix = 0 ; ix < N ; ++ix){

    // fill the upper part with zeros
    for(jx = (ix+1) ; jx < N ; ++jx){
      A_chol[ix*N + jx] = 0.0;
    }

    // compute the lower part
    for(jx = ix ; jx < N ; ++jx){
      sum = A[ix*N + jx];
      for (kx = ix-1 ; kx >= 0 ; --kx) sum -= A_chol[ix*N + kx]*A_chol[jx*N + kx];
      if (ix == jx) {
	A_chol[ix*N + ix] = sqrt(sum);
      } else A_chol[jx*N + ix] = sum/A_chol[ix*N + ix];
    }

  }
}
