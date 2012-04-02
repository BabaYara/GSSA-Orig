#ifndef _QRDCMP_H_
#define _QRDCMP_H_

typedef double REAL;

void qrdcmp(REAL* a, const int n, REAL* c, REAL* d, bool &sing);

void rsolv(REAL* a, const int n, REAL* d, REAL* b, const int ncol);

void qrsolv(REAL* a, const int n, REAL* c, REAL* d, REAL* b, const int ncol);

#endif /* _QRDCMP_H_ */
