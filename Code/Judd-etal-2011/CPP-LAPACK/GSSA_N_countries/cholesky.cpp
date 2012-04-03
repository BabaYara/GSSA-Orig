#include "global.h"
#include <math.h>

void cholesky(REAL* A, int N, REAL* A_chol)
{
  int ix, jx, kx;
  REAL sum;
  for(ix = 0 ; ix < N ; ++ix){

    // fill the upper part with zeros
    for(jx = (ix+1) ; jx < N ; ++jx){
      A_chol[ix + jx*N] = 0.0;
    }

    // compute the lower part
    for(jx = ix ; jx < N ; ++jx){
      sum = A[ix + jx*N];
      for (kx = ix-1 ; kx >= 0 ; --kx) sum -= A_chol[ix + kx*N]*A_chol[jx + kx*N];
      if (ix == jx) {
	A_chol[ix + ix*N] = sqrt(sum);
      } else A_chol[jx + ix*N] = sum/A_chol[ix + ix*N];
    }

  }
}
