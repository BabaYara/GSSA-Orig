#include "simplx.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{

  int i, j;
  int m = 5;
  int n = 14;
  int m1 = 0;
  int m2 = 0;
  int m3 = 5;
  REAL a[(m+2)*(n+1)];

  // equality coefficients
  a[16] = 4;
  a[17] = 3;
  a[18] = -4;
  a[19] = -3;
  a[31] = 7;
  a[32] = 1;
  a[33] = -7;
  a[34] = -1;
  a[46] = 5;
  a[47] = 6;
  a[48] = -5;
  a[49] = -6;
  a[61] = 4;
  a[62] = 8;
  a[63] = -4;
  a[64] = -8;
  a[76] = 9;
  a[77] = 2;
  a[78] = -9;
  a[79] = -2;
  for(i = 1 ; i < m+1 ; ++i){
    for(j = 5 ; j < n+1 ; ++j){
      if((j-4) == i){
	a[i*(n+1)+j] = 1.0;
      } else if((j-9) == i){
	a[i*(n+1)+j] = -1.0;
      } else {
	a[i*(n+1)+j] = 0.0;
      }
    }
  }

  // equality values
  a[15] = 19;
  a[30] = 27;
  a[45] = 12;
  a[60] = 15;
  a[75] = 8;

  // objective function coefficients
  a[0] = 0.0;
  a[1] = 0.0;
  a[2] = 0.0;
  a[3] = 0.0;
  a[4] = 0.0;
  for(j = 5 ; j < n+1 ; ++j) a[j] = -1.0;

  // check the matrix (print it)
  cout << endl;
  cout.precision(4);
  for(i = 0 ; i < m+1 ; ++i){
    for(j = 0 ; j < n+1 ; ++j){
      cout << setw(7) << a[i*(n+1)+j] << " ";
    }
    cout << endl;
  }
  cout << endl;  

  // simplex
  int icase;
  int izrov[n];
  int iposv[m];
  simplx(a,m,n,m1,m2,m3,icase,izrov,iposv);

  // check the output
  cout << endl;
  cout.precision(4);
  for(i = 0 ; i < m+1 ; ++i){
    for(j = 0 ; j < n+1 ; ++j){
      cout << setw(7) << a[i*(n+1)+j] << " ";
    }
    cout << endl;
  }
  cout << endl;  

  cout << endl;
  for(i = 0 ; i < m ; ++i) cout << iposv[i] << endl;
  cout << endl;

  return 0;
}
