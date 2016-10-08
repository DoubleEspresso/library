#ifndef INTEGRATE_INTERPOLATE_H
#define INTEGRATE_INTERPOLATE_H

#include <stdio.h>
#include <string.h>

// base class for interpolation  
template<class T>
class InterpolateBase
{
 public:
  InterpolateBase() {};
  virtual ~InterpolateBase() {};
  virtual T interpolate(T x)=0;  
};


// Basic lagrange point interpolation (useful for x-data which is not evenly spaced) 
template<class T>
class LagrangeInterpolate : public InterpolateBase<T>
{
 public :
  LagrangeInterpolate<T>(T*xd, T*yd, size_t size) : xdata(xd),
    ydata(yd),
    N(size) 
    {}; 
  
  ~LagrangeInterpolate<T>() { };
  virtual T interpolate(T x);

 private: 
  T*xdata;
  T*ydata;
  size_t N;
  
};


//////////////////////////////////////////////////////////////////
// base class implementation
template<class T>
T InterpolateBase<T>::interpolate(T x) {};


//////////////////////////////////////////////////////////////////
// lagrange interpolation implementation
template<class T>
T LagrangeInterpolate<T>::interpolate(T x)
{
  T result = 0;
  for (int j=0; j < N; ++j)
    {
      T product = 1;
      for (int  m=0; m < N; ++m)
	{
	  if (m == j) continue;
	  else
	    {
	      product *= T ( T(x-xdata[m]) / T( xdata[j]-xdata[m] ) );	      
	    }
	}
      result += ydata[j] * product;
    }
  return result;
}
#endif
