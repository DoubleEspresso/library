//#pragma once

#ifndef LIB_LEGENDRE_H
#define LIB_LIEGENDRE_H

#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <stdio.h>

// *unoptimized* legendre polynomial class includes
// generator formulae and weight/abcissa table computation
// for gaussian quadrature

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define PI 3.1415926535897932384626433832795028841971693993751

template<class T>
class Legendre
{
  size_t N;
  T P_0;
  T P_1;
  T P_N;
  T * rts;
  T * wts;

 public:
 Legendre<T>(size_t order) : N(order), P_0(1), P_N(0), rts(0), wts(0) 
    {
      rts = new T[N];
      wts = new T[N];
    };
  ~Legendre<T>() 
    {
      if (rts) { delete [] rts; rts = 0; }
      if (wts) { delete [] wts; wts = 0; }
    };
  
  T compute(T x); // return L_n(x) polynomial of degree n, evaluated at x.
  T derivative(T x);
  T* roots(T eps);
  void roots_weights(T* roots, T* weights, T eps);

};


template <class T>
T Legendre<T>::compute(T x)
{
  P_0 = T(1.0);
  P_1 = T(x);
  for (int j=1; j<N; ++j)
    {
      P_N = T(T(2.0*j+1.0)/T(j+1.0)*x*P_1) - T(T(j)/T(j+1.0)*P_0);
      P_0 = T(P_1);
      P_1 = T(P_N);
    }
  if (N == 0) return 1.0;
  else if (N == 1) return P_1;
  else return P_N;
}


template <class T>
T Legendre<T>::derivative(T x)
{
  if (P_N == 0) compute(x);  
  return  T (N*x*P_N - N*P_0) / T (x*x-1);
}


// unoptimized newton method to compute roots
// note this method is quite slow for eps ~ 1e-16
template <class T>
T * Legendre<T>::roots(T eps)
{
  T x1 = 0;
  T t0 = T (1.0 - T (1.0 - 1.0/N)/ T(8.0 * N* N));
  T t1 = T ( 1.0 / T(4.0 * N + 2.0) );
 
  for (int j=1; j<=N; ++j)
    {
      T x0 = cos(PI * T( ((j<<2)-1) * t1 ))*t0;
      T dx = 1;
      int iter = 0;

      do
	{
	  x0 = iter > 0 ? x1 : x0;
	  compute(x0);

	  T dpdx = derivative(x0);
	  
	  // newton update
	  x1 = x0 - P_N/dpdx;

	  dx = x0 - x1;
	  
	  iter++;
	} while ( ABS(dx) > eps);
      
      rts[(N-1)-(j-1)] = x1;
    }
  return rts;
}

// compute/return both weights and roots for gaussian quadrature applications
template <class T>
void Legendre<T>::roots_weights(T* roots, T* weights, T eps)
{
  T w0, w1;
  T x1 = 0;
  T t0 = T (1.0 - T (1.0 - 1.0/N)/ T(8.0 * N* N));
  T t1 = T ( 1.0 / T(4.0 * N + 2.0) );
 
  for (int j=1; j<=N; ++j)
    {
      T x0 = cos(PI * T( ((j<<2)-1) * t1 ))*t0;
      T dx = 1;
      T dw = 1;
      int iter = 0;

      do
	{
	  x0 = iter > 0 ? x1 : x0;
	  compute(x0);

	  T dpdx = derivative(x0);
	  
	  // newton update
	  x1 = x0 - P_N/dpdx;

	  dx = x0 - x1;
	  
	  // weights
	  w0 = (iter == 0) ? T(2.0) / T((1.0 - x0 * x0)*dpdx * dpdx) : w1;
	  w1 = T(2.0) / T((1.0 - x1*x1)*dpdx*dpdx);	  
	  dw = w0-w1;

	  iter++;

	} while ( (ABS(dx) > eps || ABS(dw) > eps) && j < 1000);
      
      rts[(N-1)-(j-1)] = x1;
      wts[(N-1)-(j-1)] = w1;

    }

  memcpy(roots,rts,N*sizeof(T));
  memcpy(weights,wts,N*sizeof(T));
}

#endif










