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
  // returns a matrix whos columns are the orthonormal-ized
  // column vectors of the input matrix
  template<typename T> Matrix<T> gram_schmidt(const Matrix<T>& input);
  
  // qr-factorization, svd, etc..
};


// modified gram-schmidt orthonormalization
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

#endif
