#ifndef LIB_GAUSS_QUADRATURE
#define LIB_GAUSS_QUADRATURE

#include "legendre.h"

// 1d integrand
typedef double (*integrand)(double,void*);


// basic implementation of 1d and 2d gaussian
// quadrature
class GaussQuadrature
{
  size_t N;
  Legendre<double> * legendre;

 public:
 GaussQuadrature(size_t order) : N(order), legendre(0)
    {
      legendre = new Legendre<double>(N);
    }
  ~GaussQuadrature()
    {
      if (legendre) { delete legendre; legendre = 0; }
    }
  bool integrate(integrand func, double xlow, double xhi, void* params, double& res, int evals);

  // 2d error norms
  double L2norm();
};

#endif
