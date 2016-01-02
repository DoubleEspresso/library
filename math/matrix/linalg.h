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
	  vj = vj - vj.dot(vi.conj()) * vi; // vi is unit, so norm(vi) is not included in projection
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


#endif
