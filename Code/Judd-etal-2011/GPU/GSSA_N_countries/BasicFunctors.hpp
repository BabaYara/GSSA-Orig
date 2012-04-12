#ifndef _BASIC_FUNCTORS_H_
#define _BASIC_FUNCTORS_H_

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

#endif //BASICFUNCTORS_H
