#ifndef INTEGRATE_GAUSS_QUADRATURE
#define INTEGRATE_GAUSS_QUADRATURE

#include "legendre.h"

// 1d integrand
typedef double (*integrand)(double, void*);

// basic implementation of 1d and 2d gaussian quadrature
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
  double integrate(integrand func, double xlow, double xhi, void* params, int evals);
  // 2d error norms
  double L2norm();
};

#endif
