#include "interpolate.h"
//#include "legendre.h"
#include "gauss_quadrature.h"
#include <math.h>

double func(double x, void * params)
{
  return 2*cos(x*x*x)*x*sin(x) + x*x*sin(x)*x + 3*x*x*x*x;
}

int main(int argc, char ** argv)
{

  // interpolation tests
  /*
  double * xdata = new double[3];
  double * ydata = new double [3];

  for (int i =0; i<3; ++i)
    {
      xdata[i] = i+1;
      ydata[i] = (i+1)*(i+1)*(i+1);
    }

  LagrangeInterpolate<double> lagrange(xdata,ydata,3);
  double yval = lagrange.interpolate(2.5);
  printf("..dbg x(2.5) = %f\n",yval);
  */


  // legendre roots/weights tests
  /*
  int nroot = 8;
  Legendre<double> legendre(nroot);
  */
  /*
  double x1 = legendre.compute(1);
  printf("val1=%f\n",x1);
  double x2 = legendre.compute(2);
  printf("val2=%f\n",x2);
  double x3 = legendre.compute(3);
  printf("val3=%f\n",x3);
  */
  /*
  double roots[nroot];
  double weights[nroot];

  legendre.roots_weights(roots, weights, 1e-15);
  for (int j=0; j<nroot; ++j)
    printf("j=%d, r=%3.15f, w=%3.15f\n",j,roots[j],weights[j]);
  
  delete [] xdata; xdata = 0;
  delete [] ydata; ydata = 0;  
  */

  // gaussian quadrature tests
  GaussQuadrature gquad(512);
  double res=0;
  gquad.integrate((integrand)&func, -2.0, 2.0, NULL, res, 512);
  printf("..integral val = %3.8f\n",res);

  return 0;
}
