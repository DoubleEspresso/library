#ifndef BASE_EXP_H
#define BASE_EXP_H

#pragma once 

#include <cmath>
#include "complex.h"
#include "constants.h"

#define C1 6.93145751953125E-1
#define C2 1.42860682030941723212E-6

typedef union
{
	double d;
	unsigned short s[4];
} ieee754;

// source: http://www.chokkan.org/software/dist/fastexp.c.html
// simple benching on cpp shows exp2 ~1.5x's faster than standard exp function
namespace QEXP
{
	inline double exp(double x)
	{
		int n;
		double a, xx, px, qx;
		ieee754 u;

		/* n = round(x / log 2) */
		a = INV_LOG_2 * x + 0.5;
		n = (int)a;
		n -= (a < 0);

		/* x -= n * log2 */
		px = (double)n;
		x -= px * C1;
		x -= px * C2;
		xx = x * x;

		/* px = x * P(x**2). */
		px = 1.26177193074810590878E-4;
		px *= xx;
		px += 3.02994407707441961300E-2;
		px *= xx;
		px += 9.99999999999999999910E-1;
		px *= x;

		/* Evaluate Q(x**2). */
		qx = 3.00198505138664455042E-6;
		qx *= xx;
		qx += 2.52448340349684104192E-3;
		qx *= xx;
		qx += 2.27265548208155028766E-1;
		qx *= xx;
		qx += 2.00000000000000000009E0;

		/* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
		x = px / (qx - px);
		x = 1.0 + 2.0 * x;

		/* Build 2^n in double. */
		u.d = 0;
		n += 1023;
		u.s[3] = (unsigned short)((n << 4) & 0x7FF0);

		return x * u.d;
	}
	inline Complex<double> expi(double x) 
	{ 
		return NULL;
	}
}
#endif
