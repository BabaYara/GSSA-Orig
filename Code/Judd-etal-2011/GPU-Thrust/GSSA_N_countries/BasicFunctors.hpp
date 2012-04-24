#ifndef _BASIC_FUNCTORS_H_
#define _BASIC_FUNCTORS_H_

#include <thrust/iterator/zip_iterator.h>
#include <thrust/for_each.h>
#include <thrust/host_vector.h>

template <typename T>
struct productivity_functor
{
  const T rho;
  const T eps;
  productivity_functor(T _rho, T _eps) : rho(_rho), eps(_eps) {}

  __host__ __device__
  T operator()(const T& a) const 
  { 
    return ((T)std::pow(a, rho))*((T)std::exp(eps));
  }
};

template <typename T>
struct consumption_functor
{
  const T A;
  const T alpha;
  const T delta;
  consumption_functor(T _A, T _alpha, T _delta) : A(_A), alpha(_alpha), delta(_delta) {}

  template <typename Tuple>
  __host__ __device__
  void operator()(Tuple t)
  { 
    // For this the zip order is:
    // 0: k
    // 1: a1
    // 2: k2
    // 3: c
    T k = thrust::get<0>(t);
    T a1 = thrust::get<1>(t);
    T k2 = thrust::get<2>(t);
    thrust::get<3>(t) = A*((T)std::pow(k, alpha))*a1 - k2 + k*(1-delta);
  }
};

template <typename T>
struct euler_functor_mc
{
  const T A;
  const T alpha;
  const T beta;
  const T delta;
  const T gam;
  euler_functor_mc(T _A, T _alpha, T _beta, T _delta, T _gam)
    : A(_A), alpha(_alpha), beta(_beta), delta(_delta), gam(_gam) {}

  template <typename Tuple>
  __host__ __device__
  void operator()(Tuple t)
  { 
    // For this the zip order is:
    // 0: k
    // 1: a1
    // 2: c
    // 3: c1
    // 4: Y
    T k = thrust::get<0>(t);
    T a1 = thrust::get<1>(t);
    T c = thrust::get<2>(t);
    T c1 = thrust::get<3>(t);
    thrust::get<4>(t) = beta*((T)std::pow(c1/c,-gam))*(1-delta+alpha*A*((T)std::pow(k,alpha-1))*a1)*k;
  }
};

template <typename T>
struct euler_functor
{
  const T A;
  const T alpha;
  const T beta;
  const T delta;
  const T gam;
  const T weight;
  euler_functor(T _A, T _alpha, T _beta, T _delta, T _gam, T _weight)
    : A(_A), alpha(_alpha), beta(_beta), delta(_delta), gam(_gam), weight(_weight) {}

  template <typename Tuple>
  __host__ __device__
  void operator()(Tuple t)
  { 
    // For this the zip order is:
    // 0: k
    // 1: a1
    // 2: c
    // 3: c1
    // 4: Y
    T k = thrust::get<0>(t);
    T a1 = thrust::get<1>(t);
    T c = thrust::get<2>(t);
    T c1 = thrust::get<3>(t);
    thrust::get<4>(t) = thrust::get<4>(t) + weight*beta*((T)std::pow(c1/c,-gam))*(1-delta+alpha*A*((T)std::pow(k,alpha-1))*a1)*k;
  }
};

template <typename T>
struct ind_consumption_functor
{
  const int M;
  const int N;
  T* c;
  ind_consumption_functor(int _M, int _N, T* _c) : M(_M), N(_N), c(_c) {}

  __host__ __device__
  void operator()(const int& row)
  { 
    int jx;
    T AggCons = 0;
    for(jx = 0 ; jx < N ; ++jx) AggCons += c[row+jx*M];
    for(jx = 0 ; jx < N ; ++jx) c[row+jx*M] = AggCons/N;    
  }
};

template <typename T>
struct Ord_Polynomial_N_functor
{
  const int n_rows;
  const int dimen;
  const int D;
  T* z;
  T* basis_fs;
  Ord_Polynomial_N_functor(int _n_rows, int _dimen, int _D, T* _z, T* _basis_fs)
    : n_rows(_n_rows), dimen(_dimen), D(_D), z(_z), basis_fs(_basis_fs) {}

  __host__ __device__
  void operator()(const int& ix)
  {

    int jx, j1x, j2x, j3x, j4x, j5x;

    // A polynomial is given by the sum of polynomial basis functions, phi(i),
    // multiplied by the coefficients; see condition (13) in JMM (2011). By 
    // convention, the first basis function is one. 
    
    //==========================================================================
    // The matrix of the basis functions - 1st degree
    //==========================================================================
    basis_fs[ix] = 1.0;
    for(jx = 1 ; jx < (dimen+1) ; ++jx){
      basis_fs[ix + jx*n_rows] = z[ix + (jx-1)*n_rows];
    }

    //========================================================================
    // 2. Second-degree columns
    //========================================================================
    
    if(D == 2){
      for(j1x = 0 ; j1x < dimen ; ++j1x){
	for(j2x = j1x ; j2x < dimen ; ++j2x){
	  basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows];
	  ++jx;
	}
      }
      
      //========================================================================
      // 3. Third-degree columns
      //========================================================================
      
    } else if(D == 3){
      for(j1x = 0 ; j1x < dimen ; ++j1x){
	for(j2x = j1x ; j2x < dimen ; ++j2x){
	  basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows];
	  ++jx;
	  for(j3x = j2x ; j3x < dimen ; ++j3x){
	    basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows];
	    ++jx;
	  }
	}
      }
      
      //========================================================================
      // 4. Fourth-degree columns
      //========================================================================
      
    } else if(D == 4){
      for(j1x = 0 ; j1x < dimen ; ++j1x){
	for(j2x = j1x ; j2x < dimen ; ++j2x){
	  basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows];
	  ++jx;
	  for(j3x = j2x ; j3x < dimen ; ++j3x){
	    basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows];
	    ++jx;
	    for(j4x = j3x ; j4x < dimen ; ++j4x){
	      basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows]*z[ix + j4x*n_rows];
	      ++jx;
	    }
	  }
	}
      }
      
      //========================================================================
      // 5. Fifth-degree columns
      //========================================================================
      
    } else if(D == 5){
      for(j1x = 0 ; j1x < dimen ; ++j1x){
	for(j2x = j1x ; j2x < dimen ; ++j2x){
	  basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows];
	  ++jx;
	  for(j3x = j2x ; j3x < dimen ; ++j3x){
	    basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows];
	    ++jx;
	    for(j4x = j3x ; j4x < dimen ; ++j4x){
	      basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows]*z[ix + j4x*n_rows];
	      ++jx;
	      for(j5x = j4x ; j5x < dimen ; ++j5x){
		basis_fs[ix + jx*n_rows] = z[ix + j1x*n_rows]*z[ix + j2x*n_rows]*z[ix + j3x*n_rows]*z[ix + j4x*n_rows]*z[ix + j5x*n_rows];
		++jx;
	      }
	    }
	  }
	}
      }
    }
    //printf("The value of basis_fs[%d] = %f\n", ix, basis_fs[ix]);
  }
};

#endif //BASICFUNCTORS_H
