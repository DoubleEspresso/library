#include <stdio.h>
#include <iostream>

#include "../../../math/base/constants.h"
#include "../../../math/base/complex.h"
#include "../../../math/fft/fft.h"


int main()
{
	fft f;

	// fft data
	size_t N = 2048;
	double * data = new double[N];
	double freq = 3;

	for (unsigned int j = 0; j<N; ++j)
	{
		data[j] = sin(PI2 * freq * ((double)((double)(j + 1) / (double)N)));
	}

	/*simple fft*/
	Complex<double> * coeffs = new Complex<double>[N];
	f.forward_plan_1d(N); // pre-computes phase factors (cheating really..)
	double start = clock();
	f.forward_1d(data, coeffs, N);
	printf("..fft time=%3.3fms\n", (clock() - start) / CLOCKS_PER_SEC * 1000.0);
	
	/*dbg - naive method comparison*/
	start = clock();
	f.forward_naive(data, coeffs, N);
	printf("..naive time=%3.3fms\n", (clock() - start) / CLOCKS_PER_SEC * 1000.0);
	
	std::cin.get();
	return 0;
}