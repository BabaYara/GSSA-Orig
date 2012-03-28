#ifndef _SIMPLX_H_
#define _SIMPLX_H_

typedef double REAL;

void simplx(REAL* a, const int m, const int n, const int m1, const int m2,
	    const int m3, int &icase, int* izrov, int* iposv);

#endif /* _SIMPLX_H_ */
