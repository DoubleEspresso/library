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
  Matrix<T> res(input);//input.nb_rows(), input.nb_cols());
  
  for (int c=0; c<input.nb_cols(); ++c)
    {
      Vector<T> vi = res.column(c);

      //      printf("\n\n--before normalization--");
      //for(int xx=0; xx<3; ++xx) printf(" (%g,%g) ", vi(xx).real, vi(xx).imag);
      // printf("\n\n");

      vi = vi.normalize(); // in place normalization ?

      // printf("\n\n--after normalization--");
      //for(int xx=0; xx<3; ++xx) printf(" (%g,%g) ", vi(xx).real, vi(xx).imag);
      //printf("\n\n");
      //printf("%d (%g,%g)\n",c,vi.norm().real,vi.norm().imag);
      // store vi
      res.set_column(c, vi);
      
      for(int j=c+1; j<input.nb_cols(); ++j)
	{
	  //printf("    %d\n",j);
	  Vector<T> vj = res.column(j);

	  vj = vj - vi.dot(vj.conj()) * vi; // vi is unit, so norm(vi) is not included in projection
	  //for(int xx=0; xx<3; ++xx) printf(" (%g,%g) ", vj(xx).real, vj(xx).imag);
	  //printf("\n");
	  //for(int xx=0; xx<3; ++xx) printf(" (%g,%g) ", vj(xx).conj().real, vj(xx).conj().imag);
	  //printf("\n");
	  res.set_column(j,vj);
	}
      if (c < input.nb_cols()-1) 
	{
	  Vector<T> vj = res.column(c+1);
	  T dot_res = vi.conj().dot(vj);
	  printf("..c(%d,%d), (%g,%g)\n",c,c+1, dot_res.real, dot_res.imag);
	}
      

    }
  return res;
}

#endif
