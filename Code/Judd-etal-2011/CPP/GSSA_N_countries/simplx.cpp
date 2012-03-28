#include "simplx.h"
#include <iostream>
#include <cmath>

using namespace std;

template<class T> inline void SWAP(T &a, T &b) {T dum=a; a=b; b=dum;}

void simp1(REAL* a, const int mm, const int n, int* ll, const int nll,
	   const int iabf, int &kp, REAL &bmax)
{
  int k;
  REAL test;
  
  if (nll <= 0)
    bmax=0.0;
  else {
    kp=ll[0];
    bmax=a[mm*(n+1)+kp];
    for (k=1;k<nll;k++) {
      if (iabf == 0)
	test=a[mm*(n+1)+ll[k]]-bmax;
      else
	test=fabs(a[mm*(n+1)+ll[k]])-fabs(bmax);
      if (test > 0.0) {
	bmax=a[mm*(n+1)+ll[k]];
	kp=ll[k];
      }
    }
  }
}

void simp2(REAL* a, const int m, const int n, int &ip, const int kp)
{
  const REAL EPS=1.0e-14;
  int k,i;
  REAL qp,q0,q,q1;
  
  ip=0;
  for (i=0;i<m;i++)
    if (a[(i+1)*(n+1)+kp] < -EPS) break;
  if (i+1>m) return;
  q1 = -a[(i+1)*(n+1)+0]/a[(i+1)*(n+1)+kp];
  ip=i+1;
  for (i=ip;i<m;i++) {
    if (a[(i+1)*(n+1)+kp] < -EPS) {
      q = -a[(i+1)*(n+1)+0]/a[(i+1)*(n+1)+kp];
      if (q < q1) {
	ip=i+1;
	q1=q;
      } else if (q == q1) {
	for (k=0;k<n;k++) {
	  qp = -a[ip*(n+1)+k+1]/a[ip*(n+1)+kp];
	  q0 = -a[i*(n+1)+k+1]/a[i*(n+1)+kp];
	  if (q0 != qp) break;
	}
	if (q0 < qp) ip=i+1;
      }
    }
  }
}

void simp3(REAL* a, const int n, const int i1, const int k1,
	   const int ip, const int kp)
{
  int ii,kk;
  REAL piv;
  
  piv=1.0/a[ip*(n+1)+kp];
  for (ii=0;ii<i1+1;ii++)
    if (ii != ip) {
      a[ii*(n+1)+kp] *= piv;
      for (kk=0;kk<k1+1;kk++)
	if (kk != kp)
	  a[ii*(n+1)+kk] -= a[ip*(n+1)+kk]*a[ii*(n+1)+kp];
    }
  for (kk=0;kk<k1+1;kk++)
    if (kk != kp) a[ip*(n+1)+kk] *= -piv;
  a[ip*(n+1)+kp]=piv;
}

void simplx(REAL* a, const int m, const int n, const int m1, const int m2,
	    const int m3, int &icase, int* izrov, int* iposv)
{
  const REAL EPS=1.0e-14;
  int i,k,ip,is,kh,kp,nl1;
  REAL q1,bmax;
  
  if (m != (m1+m2+m3)) {
    cout << "Bad input constraint counts in simplx" << endl;
    return;
  }
  int l1[n+1],l3[m];
  nl1=n;
  for (k=0;k<n;k++) {
    l1[k]=k+1;
    izrov[k]=k;
  }
  for (i=1;i<=m;i++) {
    iposv[i-1]=n+i-1;
  }
  if (m2+m3 != 0) {
    for (i=0;i<m2;i++) l3[i]=1;
    for (k=0;k<(n+1);k++) {
      q1=0.0;
      for (i=m1+1;i<m+1;i++) q1 += a[i*(n+1)+k];
      a[(m+1)*(n+1)+k] = -q1;
    }
    for (;;) {
      simp1(a,m+1,n,l1,nl1,0,kp,bmax);
      if (bmax <= EPS && a[(m+1)*(n+1)+0] < -EPS) {
	icase = -1;
	return;
      } else if (bmax <= EPS && a[(m+1)*(n+1)+0] <= EPS) {
	for (ip=m1+m2+1;ip<m+1;ip++) {
	  if (iposv[ip-1] == (ip+n-1)) {
	    simp1(a,ip,n,l1,nl1,1,kp,bmax);
	    if (bmax > EPS)
	      goto one;
	  }
	}
	for (i=m1+1;i<=m1+m2;i++)
	  if (l3[i-m1-1] == 1)
	    for (k=0;k<n+1;k++)
	      a[i*(n+1)+k]= -a[i*(n+1)+k];
	break;
      }
      simp2(a,m,n,ip,kp);
      if (ip == 0) {
	icase = -1;
	return;
      }
    one:	simp3(a,n,m+1,n,ip,kp);
      if (iposv[ip-1] >= (n+m1+m2)) {
	for (k=0;k<nl1;k++)
	  if (l1[k] == kp) break;
	--nl1;
	for (is=k;is<nl1;is++) l1[is]=l1[is+1];
      } else {
	kh=iposv[ip-1]-m1-n+1;
	if (kh >= 1 && l3[kh-1] != 0) {
	  l3[kh-1]=0;
	  ++a[(m+1)*(n+1)+kp];
	  for (i=0;i<m+2;i++)
	    a[i*(n+1)+kp]= -a[i*(n+1)+kp];
	}
      }
      SWAP(izrov[kp-1],iposv[ip-1]);
    }
  }
  for (;;) {
    simp1(a,0,n,l1,nl1,0,kp,bmax);
    if (bmax <= EPS) {
      icase=0;
      return;
    }
    simp2(a,m,n,ip,kp);
    if (ip == 0) {
      icase=1;
      return;
    }
    simp3(a,n,m,n,ip,kp);
    SWAP(izrov[kp-1],iposv[ip-1]);
  }
}
