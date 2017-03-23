#include <iostream>
#include <string.h>
#include <stdio.h>
#include <cmath> 
#include <time.h>

#include "../base/complex.h"
#include "../base/exp.h"


class fft
{
public:
	fft() {}
	~fft() {}

	/*naive dft for comparison*/
	void forward_naive(const double * data, Complex<double> * coeffs, const size_t N);

	/*fft-algorithms*/
	virtual void forward_1d(const double * data, Complex<double> * coeffs, const size_t N);
	//virtual void backward_1d(complex * coeffs, double * data, size_t N); 
};
