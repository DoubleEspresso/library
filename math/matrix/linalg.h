#pragma once

#ifndef MATH_LINALG_H
#define MATH_LINALG_H

#include <cassert>
#include <string.h>
#include <vector>

#include "matrix.h"

namespace LinearAlgebra
{
  // modified gram-schmidt orthogonalization
  // returns a matrix whose columns are the orthonormal-ized
  // column vectors of the input matrix
  template<typename T> Matrix<T> gram_schmidt(const Matrix<T>& input);
  
  // qr-factorization
  template<typename T> void QR(const Matrix<T>& input, Matrix<T>& Q, Matrix<T>& R);

  // hessenberg matrix, for symmetric matrices, this will return the householder transformation
  // matrix, which can be used to convert symmetric matrices to tri-diagonal form
  template<typename T> Matrix<T> hessenberg_form(const Matrix<T>& input, const int col);

  // given an nxn symmetric matrix A, this method reduces A to tridiagonal form using n-2
  // orthogonal transmformations (householder method)
  template<typename T> Matrix<T> tridiagonal_householder(const Matrix<T>& input);
  
  // svd, ..
};


// modified gram-schmidt orthonormalization for numerical stability
// input: mxn matrix 
// output: mxn matrix whose columns are orthonormal
// todo -- parallelize over cols for large matrices
template<typename T>
inline Matrix<T> LinearAlgebra::gram_schmidt(const Matrix<T>& input)
{
  Matrix<T> res(input);  
  for (int c=0; c<input.nb_cols(); ++c)
    {
      Vector<T> vi = res.column(c);
      vi = vi.normalize();
      res.set_column(c, vi);
      
      for(int j=c+1; j<input.nb_cols(); ++j)
	{
	  Vector<T> vj = res.column(j);
	  vj = vj - vj.dot(vi) * vi; // vi is unit, so norm(vi) is not included in projection
	  res.set_column(j,vj);
	}
    }
  return res;
}

// qr-decomposition of nxn matrix 
// uses modified gram-schmidt to compute matrix Q
template<typename T>
inline void LinearAlgebra::QR(const Matrix<T>& input, Matrix<T>& Q, Matrix<T>& R)
{
  assert(input.nb_rows() == input.nb_cols());
  Matrix<T> res(input);  
  
  Q = gram_schmidt(input);
  Matrix<T> Qt = Q.transpose().conj();
  R = Qt * input;
}

// the hessenberg form of a square matrix .. for symmetric matrices
// this method returns the householder transformation matrix
template<typename T>
inline Matrix<T> LinearAlgebra::hessenberg_form(const Matrix<T>& input, const int col)
{
  assert(input.nb_rows() == input.nb_cols() && col < input.nb_cols());

  // Remarks/notes: To compute the Hessenberg form (or householder matrix) P .. we must find x,y pair such that, P.x = y,
  // x is given as the nth-col of the input matrix, and y can be computed from the property P.x = y. This argument assumes one knows
  // that P generally takes the form of a reflection operator about the hyperplane defined by unit vector w ... ie P = I - 2.0 * w * wtranspose
  Vector<T> x = input.column(col);
  Vector<T> y(x.nb_rows());

  T ynorm = T(0);
  for (int j=0; j<col+1;++j)
    {
      y.set(j,x(j));
      ynorm += x(j)*x(j);
    }
  T x2 = x.dot(x);
  T remainder = x2 - ynorm; remainder = sqrt(remainder); //remainder.sqrt();
  y.set(col+1,remainder);
  
  // since P is the reflection operator about the hyperplane whose orientation is given by w,
  // and y-defined ..  w can be computed from ==> w = (x-y)/norm(x-y), and P = I - 2*w.wtranspose()
  Vector<T> w2 = x - y;
  w2 = w2.normalize();

  Matrix<T> w(input.nb_rows(),1);
  w.set_column(0,w2);

  Matrix<T> P(input.nb_rows(), input.nb_cols());
  P = P.identity();
  Matrix<T> wt = w.transpose();

  // TODO: floating point exception without the braces around (w * wt) ??
  // For symmetric input matrices, P is the householder transformation, for general input, P is the hessenberg matrix
  // define w to be a matrix of column nb 1 -- technically a vector
  P = P - 2.0 * (w * wt);  
  return P;
}

// TODO: speedups as outlined in NR 11.2.10
template<typename T>
Matrix<T> LinearAlgebra::tridiagonal_householder(const Matrix<T>& input)
{
  assert(input.nb_rows() == input.nb_cols());
  int n = input.nb_rows();

  Matrix<T> A(input);
  Matrix<T> prev(input);
  
  for(int i=0; i<n-2; ++i)
    {
      A = hessenberg_form(A,i);
      A = A * prev * A;
      prev = A;
    }
  return A;
}

#endif
