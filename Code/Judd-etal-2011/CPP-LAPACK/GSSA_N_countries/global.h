#ifndef __FILE_GLOBAL_H_SEEN__
#define __FILE_GLOBAL_H_SEEN__

//#include <Accelerate/Accelerate.h>
#include "mkl.h"
#include "mkl_cblas.h"

typedef double REAL;
typedef unsigned long int Ullong;
typedef unsigned int Uint;

// to determine whether single or double precision is being used
extern const float singletype;
extern const double doubletype;
extern const REAL realtype;

void Productivity(int T, int N, REAL* a_init, REAL sigma, REAL rho);
REAL* Ord_Polynomial_N(REAL* z, int n_rows, int dimen, int D, int& n_cols);
void GH_Quadrature(int Qn, int N, REAL* vcv, int n_nodes, REAL* epsi_nodes,
		   REAL* weight_nodes);
void Monomials_1(int N, REAL* vcv, int n_nodes, REAL* epsi_nodes,
		 REAL* weight_nodes);
void Monomials_2(int N, REAL* vcv, int n_nodes, REAL* epsi_nodes,
		 REAL* weight_nodes);
REAL* Num_Stab_Approx(int T, int n, REAL* X, int N, REAL* Y, int RM,
		      int penalty, int normalize);
void Accuracy_Test_N(int P, int N, REAL* k, REAL* a, REAL* bk, int D,
		     int IM, REAL alpha, REAL gam, REAL delta, REAL beta,
		     REAL A, REAL tau, REAL rho, REAL* vcv, int discard,
		     REAL& Errors_mean, REAL& Errors_max, REAL& time_test);

// from auxfuncs.cpp
REAL mean(REAL* x, const int N);
REAL var(REAL* x, const int N);
REAL sd(REAL* x, const int N);
REAL norm(REAL* x, const int N);
void matrix_normalize(const int M, const int N, REAL* X, const int ldx,
		      REAL* X_norm, REAL* X_means, REAL* X_sds);
void print_matrix(const int M, const int N, REAL* X,
		  const int printrows, const int printcols);
void print_vector(const int N, REAL* X);

#endif
