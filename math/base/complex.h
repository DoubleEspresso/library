#ifndef LIB_MATHBASE_COMPLEX
#define LIB_MATHBASE_COMPLEX

#include <cmath>

template<typename T>
struct Complex
{
  T real; T imag;
  Complex(T r, T i) : real(r), imag(i) {};
  Complex() : real(0), imag(0) {};
  Complex(const Complex<T>& other) : real(other.real), imag(other.imag) {};
  Complex(const T& other) : real(other), imag(0) {};
  ~Complex() {};
  
  Complex operator=(const Complex& other) { this->imag = other.imag; this->real = other.real; }

  // basic operations
    Complex conj() { return Complex<T>(real, -imag); }
  T norm() { return T(real*real+imag*imag); }
  T abs() { return T(sqrt( norm() )); }
  Complex sqrt();
  Complex operator+(const Complex& other);
  Complex operator-(const Complex& other);
  Complex operator*(const Complex& other);
  Complex operator/(const Complex& other);
  
  // increment/decrement operators..
  Complex operator+=(const Complex& other);
  Complex operator-=(const Complex& other);
  Complex operator*=(const Complex& other);
  Complex operator/=(const Complex& other);

  // operators for complex types and other scalar types
  Complex operator+(const T& other);
  Complex operator-(const T& other);
  Complex operator*(const T& other);
  Complex operator/(const T& other);
  Complex operator+=(const T& other);
  Complex operator-=(const T& other);
  Complex operator*=(const T& other);
  Complex operator/=(const T& other);

};

// The principal root of a complex number -- 
// using polar form z^0.5 = sqrt(r)*exp{i*arctan(theta)*0.5} --> euler's formula exp{i*theta} = cos + i sin and half angle formulas 
// will also work for real types .. FIXME
template<typename T>
inline Complex<T> Complex<T>::sqrt()
{
  if (real == 0) return Complex<T>(0,std::sqrt(imag));
  
  double r = std::sqrt(real*real+imag*imag);
  double rsqrt = std::sqrt(r);
  double theta = atan(imag / real);
  return Complex<T>(rsqrt * cos(theta * 0.5), rsqrt * sin(theta * 0.5));
  /*
  double cos = double( real / double(r) );
  double sin = double( imag / double(r) );

  double cos_half = std::sqrt(0.5 * (cos+1));
  double sin_half = std::sqrt(0.5 * (-cos+1));
  return Complex<T>(rsqrt * cos_half, rsqrt * sin_half);
  */
  
}

// complex addition
template<typename T>
inline Complex<T> Complex<T>::operator+(const Complex<T>& other)
{
  return Complex<T>(this->real+other.real, this->imag+other.imag);
}

// complex subtraction
template<typename T>
inline Complex<T> Complex<T>::operator-(const Complex<T>& other)
{
  return Complex<T>(this->real-other.real, this->imag-other.imag);
}

// complex multiplication
template<typename T>
inline Complex<T> Complex<T>::operator*(const Complex<T>& other)
{
  return Complex<T>(this->real*other.real-this->imag*other.imag, 
		    this->imag*other.real + this->real*other.imag);
}

// complex division
template<typename T>
inline Complex<T> Complex<T>::operator/(const Complex<T>& other)
{
  Complex<T> c(other);
  T n = c.norm();
  return Complex<T>( T(this->real*other.real + this->imag*other.imag) / n, 
		     T(this->imag*other.real - this->real*other.imag) / n);
}

template<typename T>
inline Complex<T> Complex<T>::operator+=(const Complex& other)
{
  *(this) = (*this) + other;
  return *(this);
}
template<typename T>
inline Complex<T> Complex<T>::operator-=(const Complex& other)
{
  *(this) = (*this) - other;
  return *(this);
}
template<typename T>
inline Complex<T> Complex<T>::operator*=(const Complex& other)
{
  *(this) = (*this) * other;
  return *(this);
}
template<typename T>
inline Complex<T> Complex<T>::operator/=(const Complex& other)
{
  *(this) = (*this) / other;
  return *(this);
}

// ======================================================================= //
// overloads for other types with Complex class (double, float, int..etc)
template<typename T>
inline Complex<T> Complex<T>::operator/=(const T& other)
{
  this->real = T(this->real / T(other));
  this->imag = T(this->imag / T(other));
  return *(this);
}
template<typename T>
inline Complex<T> Complex<T>::operator+(const T& other)
{
  return Complex<T>(this->real + other, this->imag);
}
template<typename T>
inline Complex<T> Complex<T>::operator-(const T& other)
{
  return Complex<T>(this->real - other, this->imag);
}
template<typename T>
inline Complex<T> Complex<T>::operator*(const T& other)
{
  return Complex<T>(this->real * other, this->imag * other);
}
template<typename T>
inline Complex<T> Complex<T>::operator/(const T& other)
{
  return Complex<T>( T(this->real / T(other)), T(this->imag / T(other)));
}
template<typename T>
inline Complex<T> Complex<T>::operator+=(const T& other)
{
  this->real = this->real + other;
  return *(this);
}
template<typename T>
inline Complex<T> Complex<T>::operator-=(const T& other)
{
  this->real = this->real - other;
  return *(this);
}
template<typename T>
inline Complex<T> Complex<T>::operator*=(const T& other)
{
  this->real = this->real * other;
  this->imag = this->imag * other;
  return *(this);
}
// ======================================================================= //
// the commuting operations
template<typename T>
inline Complex<T> operator/(const T& x, const Complex<T>& other)
{
  return other.conj()*x/other.norm();
}
template<typename T>
inline Complex<T> operator*(const T& x, const Complex<T>& other)
{
  return Complex<T>(other) * T(x);
}
template<typename T>
inline Complex<T> operator+(const T& x, const Complex<T>& other)
{
  return Complex<T>(other) + T(x);
}
template<typename T>
inline Complex<T> operator-(const T& x, const Complex<T>& other)
{
  return Complex<T>(other) - T(x);
}

// typedefs
typedef Complex<double> Complex_d;
typedef Complex<float> Complex_f;
typedef Complex<int> Complex_i;


#endif
