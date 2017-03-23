#include "fft.h"

/*vanilla 1d fft following Cooley-Tookey algorithm*/
void fft::forward_1d(const double * data, Complex<double> * coeffs, const size_t N)
{
	/*assume n divisible by 2*/
	size_t sz = N / 2;

	Complex<double> * E = new Complex<double>[sz];
	Complex<double> * O = new Complex<double>[sz];

	for (unsigned int k = 0; k < sz; ++k)
	{
		Complex<double> xk;
		Complex<double> ek;
		Complex<double> ok;

		for (unsigned int i = 0, j = 0; i < sz; ++i, j += 2)  // i in [0,N/2-1]
		{
			Complex<double> r; r = r.root(-i * k, sz);

			ek += r * data[j];   // even indices
			ok += r * data[j + 1]; // odd indices	  
		}
		E[k] = ek;
		O[k] = ok;
	}

	for (unsigned int i = 0, j = sz; i < sz; ++j, ++i)
	{
		Complex<double> r; r = r.root(-i, N);
		coeffs[i] = E[i] + r * O[i];
		coeffs[j] = E[i] - r * O[i];
	}

	if (E) { delete[] E; E = 0; }
	if (O) { delete[] O; O = 0; }
}

void fft::forward_naive(const double * data, Complex<double> * coeffs, const size_t N)
{
	for (unsigned int k = 0; k < N; ++k)
	{
		Complex<double> xk;
		for (unsigned int i = 0; i < N; ++i)
		{
			Complex<double> r; r = r.root(-i * k, N);
			xk += r * data[i];
		}
		coeffs[k] = xk;
	}
}


