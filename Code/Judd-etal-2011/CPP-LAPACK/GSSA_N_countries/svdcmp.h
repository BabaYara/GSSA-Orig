#ifndef _SVDCMP_H_
#define _SVDCMP_H_

#include "global.h"

void svdcmp(REAL* a, int m, int n, REAL* w, REAL* v);

void svbksb(REAL* u, const int m, const int n, REAL* w,
	    REAL* v, REAL* b, REAL* x, const int ncol);

#endif /* _SVDCMP_H_ */
