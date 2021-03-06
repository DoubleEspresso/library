#include "gauss_quadrature.h";
#include <math.h>
#include <stdio.h>


// 1d integration by gaussian quadrature
double GaussQuadrature::integrate(integrand func, double xlow, double xhi, void* params, int evals)
{
	double x = 0;
	double A, B, Ax, y;
	int pow2 = (evals + 1) >> 1;
	bool odd = evals & 1;
	// limits
	A = 0.5 * (xhi - xlow);
	B = 0.5 * (xhi + xlow);

	double * roots = new double[N]; double *  weights = new double[N];
	legendre->roots_weights(roots, weights, 1e-15);

	y = odd ? weights[0] * func(B, params) : 0;

	for (int i = (odd ? 1 : 0); i < pow2; ++i)
	{
		Ax = A*roots[i];
		y += weights[i] * (func(B + Ax, params) + func(B - Ax, params));
	}
	if (roots) { delete[] roots; roots = 0; }
	if (weights) { delete[] weights; weights = 0; }
	return A*y;
}
