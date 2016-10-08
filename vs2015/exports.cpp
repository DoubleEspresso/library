#include "../utils/definitions.h"

/*math exports*/
#include "../math/integrate/gauss_quadrature.h"
extern "C"
{
	/*integration exports*/
	DllExport GaussQuadrature * GaussQuadratureInstance(size_t N) { return new GaussQuadrature(N); }
	DllExport bool DECL GaussQuadratureIntegrate(
		void * func,
		double xlow,
		double xhi,
		void * params,
		double& res,
		int evals,
		GaussQuadrature * handle) {
		return handle->integrate(*((integrand*)func), xlow, xhi, params, res, evals);
	}
}