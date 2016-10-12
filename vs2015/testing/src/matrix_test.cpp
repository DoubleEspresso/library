#include <stdio.h>
#include <iostream>

#include "../../../math/base/constants.h"
#include "../../../math/base/complex.h"
#include "../../../math/matrix/matrix.h"
#include "../../../math/matrix/linalg.h"

int main()
{
	/*Test LU decomposition*/
	int nb_passed = 0; int tests = 500;
	for (int j = 2; j < tests; ++j)
	{
		uint N = j;
		Matrix<double> M(N, N);
		for (int r = 0; r < N; ++r)
		for (int c = 0; c < N; ++c)
			M.set(r, c, log(1+r*c) + sin(r) * sin(r) + 0.5 + 2 * exp(- r - c) + cos(c) * cos(c) ); // expression likely to be invertible

		Matrix<double> D(M); Matrix<double> L(N, N); Matrix<double> U(N, N); Matrix<double> R(N, N); Matrix<double> I(N, N); I = I.identity();
		LinearAlgebra::LU(M, D);
		D.upper_triangle(U);
		D.lower_triangle(L); L = L + I;
		R = L * U;
		Matrix<double> S = M - R;
		bool isOK = true; int idx = 0;
		for (int j = 0; j < N*N; ++j) if (S.data_at(j) >= 1e-6) { isOK = false; idx = j;  break; }

		if (isOK) ++nb_passed;
		else printf("..test(%dx%d) failed: M(%d)=%.3f R(%d)=%.3f\n", j, j, idx, M.data_at(idx), idx, R.data_at(idx));
		if (j % 50 == 0) printf("..test size (%dx%d)\n", j, j);
	}
	printf("..%d/%d passed\n", nb_passed, tests-2);
	std::cin.get();
	return 0;
}