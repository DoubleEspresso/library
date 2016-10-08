#include "definitions.h"

/*math exports*/
#include "../math/integrate/gauss_quadrature.h"
extern "C"
{
	/*integration exports*/
	DllExport GaussQuadrature * DECL GaussQuadratureInstance(size_t N) { return new GaussQuadrature(N); }
	DllExport double DECL GaussQuadratureIntegrate(
		void * func,
		double xlow,
		double xhi,
		int evals,
		GaussQuadrature * handle) {
		return handle->integrate(((integrand)func), xlow, xhi, 0, evals);
	}
}