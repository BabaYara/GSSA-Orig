#include "qrdcmp.h"
#include <cmath>

using namespace std;

template<class T>
inline const T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

template<class T>
inline const T MAX(const T &a, const T &b) {return b > a ? (b) : (a);}

template<class T>
inline const T MIN(const T &a, const T &b) {return b < a ? (b) : (a);}

template<class T> inline const T SQR(const T a) {return a*a;}

void qrdcmp(REAL* a, const int n, REAL* c, REAL* d, bool &sing)
{
  int i,j,k;
  REAL scale,sigma,sum,tau;

  sing=false;
  for (k=0;k<n-1;k++) {
    scale=0.0;
    for (i=k;i<n;i++) scale=MAX(scale,fabs(a[i*n+k]));
    if (scale == 0.0) {
      sing=true;
      c[k]=d[k]=0.0;
    } else {
      for (i=k;i<n;i++) a[i*n+k] /= scale;
      for (sum=0.0,i=k;i<n;i++) sum += SQR(a[i*n+k]);
      sigma=SIGN(sqrt(sum),a[k*n+k]);
      a[k*n+k] += sigma;
      c[k]=sigma*a[k*n+k];
      d[k] = -scale*sigma;
      for (j=k+1;j<n;j++) {
	for (sum=0.0,i=k;i<n;i++) sum += a[i*n+k]*a[i*n+j];
	tau=sum/c[k];
	for (i=k;i<n;i++) a[i*n+j] -= tau*a[i*n+k];
      }
    }
  }
  d[n-1]=a[(n-1)*n+n-1];
  if (d[n-1] == 0.0) sing=true;
}


void rsolv(REAL* a, const int n, REAL* d, REAL* b, const int ncol)
{
  int h,i,j;
  REAL sum;
  
  for(h=0;h<ncol;h++){
    b[(n-1)*ncol+h] /= d[n-1];
    for (i=n-2;i>=0;i--) {
      for (sum=0.0,j=i+1;j<n;j++) sum += a[i*n+j]*b[j*ncol+h];
      b[i*ncol+h]=(b[i*ncol+h]-sum)/d[i];
    }
  }
}

void qrsolv(REAL* a, const int n, REAL* c, REAL* d, REAL* b, const int ncol)
{
  int h,i,j;
  REAL sum,tau;
  
  for(h=0;h<ncol;h++){
    for (j=0;j<n-1;j++) {
      for (sum=0.0,i=j;i<n;i++) sum += a[i*n+j]*b[i*ncol+h];
      tau=sum/c[j];
      for (i=j;i<n;i++) b[i*ncol+h] -= tau*a[i*n+j];
    }
  }
  rsolv(a,n,d,b,ncol);
}
