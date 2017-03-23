#pragma once
#ifndef FFT_FFT_H
#define FFT_FFT_H

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <cmath> 
#include <time.h>

#include "../base/complex.h"
#include "../base/exp.h"


class fft
{
	Complex<double> * phases;
	Complex<double> * twiddles;
public:
	fft() {}
	~fft() 
	{
		if (phases) { delete[] phases; phases = 0; }
		if (twiddles) { delete[] twiddles; twiddles = 0; }
	}

	/*naive dft for comparison*/
	void forward_naive(const double * data, Complex<double> * coeffs, const size_t N);

	/*fft-algorithms*/
	virtual void forward_plan_1d(const size_t N);
	virtual void forward_1d(const double * data, Complex<double> * coeffs, const size_t N);
	//virtual void backward_1d(complex * coeffs, double * data, size_t N); 
};

/*follow fftw "plan" notation to pre-compute twiddle factors for 1d fft*/
void fft::forward_plan_1d(const size_t N)
{
	size_t sz = N / 2;
	Complex<double> r;

	//if (phases) { delete[] phases; phases = 0; }
	//if (twiddles) { delete[] twiddles; twiddles = 0; }

	phases = new Complex<double>[sz * sz];
	twiddles = new Complex<double>[sz];

	for (unsigned int k = 0, idx = 0; k < sz; ++k)
	{
		for (unsigned int i = 0; i < sz; ++i, ++idx)  // i in [0,N/2-1]
		{
			phases[idx] = r.root(-i * k, sz); // exp[ i * k / (N/2) ]
		}
	}

	// we also need exp[ i / N ] .. not straightforward to re-use phases[..] for that
	for (unsigned int i = 0; i < sz; ++i)
	{
		twiddles[i] = r.root(-i, N); // exp[ i * k / (N/2) ]
	}
}


/*vanilla 1d fft following Cooley-Tookey algorithm*/
void fft::forward_1d(const double * data, Complex<double> * coeffs, const size_t N)
{
	/*assume n divisible by 2*/
	size_t sz = N / 2;

	Complex<double> * E = new Complex<double>[sz];
	Complex<double> * O = new Complex<double>[sz];

	for (unsigned int k = 0, idx = 0; k < sz; ++k)
	{
		Complex<double> xk;
		Complex<double> ek;
		Complex<double> ok;

		for (unsigned int i = 0, j = 0; i < sz; ++i, j += 2, ++idx)  // i in [0,N/2-1]
		{
			ek += phases[idx] * data[j];   // even indices
			ok += phases[idx] * data[j + 1]; // odd indices	  
		}
		E[k] = ek;
		O[k] = ok;
	}

	for (unsigned int i = 0, j = sz; i < sz; ++j, ++i) // i in [0,N/2-1]
	{
		coeffs[i] = E[i] + twiddles[i] * O[i];
		coeffs[j] = E[i] - twiddles[i] * O[i];
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

#endif
