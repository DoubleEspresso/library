#pragma once

#ifndef MATH_MATRIX_H
#define MATH_MATRIX_H

#include <cassert>
#include <string.h>
#include <vector>
#include <math.h>

#include "../../concurrent/threads.h"
#include "../../utils/timer.h"
#include "../../system/hardware.h"

template <typename T> class Matrix;

template<typename T>
class Vector 
{
 public:
 Vector(int rws) :  rows(rws), size(rws), parallel(false), nb_threads(1)
    {
      data = new T[size];
      memset(data,0,sizeof(T)*size);
      parallel = (rows >= 275);
      if (parallel)
	{
	  bool hyperthreading = false;
	  nb_threads = Hardware::cpu_count(hyperthreading);
	  nb_threads *= (hyperthreading ? 2 : 1);	  
	}      
    }
  
  Vector(const Matrix<T>& other, int r, int c, bool is_col=true)
    {
      data = 0;
      rows = (is_col ? other.nb_rows() : other.nb_cols());
      size = rows;
      data = new T[size];
      for (int j=0; j<size; ++j) data[j] = other(r+(is_col?j:0),c+(!is_col?j:0));
    }
  Vector<T>(const Vector& other)
    {
      data = 0;
      rows = other.rows;
      size = rows;
      data = new T[size];
      memcpy(this->data, other.data, sizeof(T) * size);      
    }
  ~Vector() { if (data) { delete[] data; data = 0; }}

  T operator()(int r) const { return data[r]; }  
  void operator=(const Vector& other) { memcpy(this->data, other.data, sizeof(T) * size); }
  void set(int j, T val) { data[j] = T(val); }
  int nb_rows() const { return rows; }
  int nb_cols() const { return rows; }
  void set_threads(int j) { nb_threads=j; }
  int get_threads() { return nb_threads; }
  T dot(const Vector& other);
  T norm();
  Vector normalize();
  Vector cross(const Vector& other);
  Vector conj();

  // vector-vector operations
  Vector operator+(const Vector& other);
  Vector operator-(const Vector& other);
  Vector operator+=(const Vector& other);
  Vector operator-=(const Vector& other);
  Vector operator*=(const T& other);
  Vector operator*(const T& other);
  Vector operator/(const T& other);

  // parallel operator overloads
  static void vec_parallel_multiply(void * data);
  static void vec_parallel_add(void * data);
  static void vec_parallel_sub(void * data);

 private:
  int rows;
  size_t size;
  T * data;
  bool parallel;
  int nb_threads;
};

// vector-vector operators -- todo parallelize
template<typename T>
T Vector<T>::dot(const Vector& other)
{
  assert( this->rows == other.nb_rows());
  T res = 0;
  for(int j=0; j<this->rows; ++j) res += (*this)(j) * other(j);
  return res;
}

template<typename T>
Vector<T> Vector<T>::operator/(const T& other)
{
  Vector res(*this);
  for(int j=0; j<this->rows; ++j) res.set(j, data[j] / other);
  return res;
}

template<typename T>
Vector<T> Vector<T>::operator*(const T& other)
{
  Vector res(*this);
  for(int j=0; j<rows; ++j) res.set(j, data[j] * other);
  return res;
}

template<typename T>
Vector<T> Vector<T>::operator-(const Vector& other)
{
  Vector res(*this);
  for (int j=0; j<rows; ++j) res.set(j, data[j] - other(j));
  return res;
}

template<typename T>
inline Vector<T> operator*(const T& other, Vector<T>& vec)
{
  return vec * other;
}

template<typename T>
inline Vector<T> Vector<T>::conj()
{
  Vector res(*this);
  for(int j=0; j<rows; ++j) res(j) = (*this)(j).conj();
  return res;
}

// generic implementation for complex types
template<typename T>
T Vector<T>::norm()
{
  T tmp = (*this).conj().dot(*this);
  return tmp.sqrt();
}

// specialized for int,float, and double data types
template<>
inline int Vector<int>::norm()
{
  int tmp = this->dot(*this);
  return sqrt(tmp);
}
template<>
inline float Vector<float>::norm()
{
  float tmp = this->dot(*this);
  return sqrt(tmp);
}
template<>
inline double Vector<double>::norm()
{
  double tmp = this->dot(*this);
  return sqrt(tmp);
}

template<typename T>
Vector<T> Vector<T>::normalize()
{
  //T tmp = this->dot(*this);
  return (*this) / this->norm();
}

template<typename T>
class Matrix
{
 public:
  Matrix() {};
  Matrix(const Vector<T>& other)
    {
      size = other.nb_rows();
      cols = 0;
      rows = other.nb_rows();
      data = new T[size];
      memset(data,0,sizeof(T)*size);
      parallel = (size >= 275);
      if (parallel)
	{
	  bool hyperthreading = false;
	  nb_threads = Hardware::cpu_count(hyperthreading);
	  nb_threads *= (hyperthreading ? 2 : 1);	  
	}
    }
 Matrix(int r, int c) : rows(r), cols(c), size(rows*cols), parallel(false), nb_threads(1)
    {
      data = new T[size];
      memset(data,0,sizeof(T)*size);
      scan();
      if (parallel)
	{
	  bool hyperthreading = false;
	  nb_threads = Hardware::cpu_count(hyperthreading);
	  nb_threads *= (hyperthreading ? 2 : 1);	  
	}
    }
  Matrix(const Matrix& other)
    {
      data = 0;
      rows = other.rows;
      cols = other.cols;
      size = rows * cols;
      data = new T[size];
      memcpy(this->data, other.data, sizeof(T) * size);      
    }

  ~Matrix() { if (data) { delete[] data; data = 0; }}

  T operator()(int r, int c) const { return data[r*cols + c]; }
  
  void operator=(const Matrix& other) { memcpy(this->data, other.data, sizeof(T) * size); }

  T data_at(int j) { return T(data[j]); } 
  void set(int r, int c, T val) { data[r*cols + c] = T(val); }
  void set(int j, T val) { data[j] = val; }
  T trace();
  void scan() { parallel  = (size >= 275*275); } // for example
  Matrix zeros();
  void clear() { if (data) memset(data, 0, sizeof(T) * size); }
  Matrix identity();
  bool is_square() { return rows == cols; } 
  int nb_rows() const { return rows; }
  int nb_cols() const { return cols; }
  bool resize();
  Matrix inverse();
  T * diag();
  T * eigenvals();
  T * eigenvecs();
  T * get_data() const { return data; }
  void set_threads(int j) { nb_threads=j; }
  int get_threads() { return nb_threads; }
  Matrix minor(int r, int c);
  Vector<T> column(int c) const;
  void set_column(int c, const Vector<T>& vin) const;
  Matrix transpose();
  
  // matrix-matrix operations
  Matrix operator+(const Matrix& other);
  Matrix operator-(const Matrix& other);
  Matrix operator*(const Matrix& other);
  Matrix operator+=(const Matrix& other);
  Matrix operator-=(const Matrix& other);
  Matrix operator*=(const Matrix& other);
  Matrix operator*(const T& other);

  // matrix-vector operations
  Vector<T> operator*(const Vector<T>& other);

  // parallel operator overloads
  static void parallel_multiply(void * data);
  static void parallel_add(void * data);
  static void parallel_sub(void * data);

  // matrix vector operations

 private:
  int rows;
  int cols;
  size_t size;
  T * data;
  bool parallel;
  int nb_threads;
};

template<typename T>
struct parallel_data
{
  const Matrix<T> * thism;
  const Matrix<T> * other;
  T * thread_data;
  int start;
  int stop;
  size_t arr_size;
};

template<typename T>
void Matrix<T>::parallel_multiply(void * data)
{
  parallel_data<T> * md = (parallel_data<T>*) data;
  Matrix<T> * mr  = (Matrix<T>*) md->other;
  Matrix<T> * ml  = (Matrix<T>*) md->thism;

  assert(ml->nb_cols() == mr->nb_rows());
  
  int start_r = md->start;
  int stop_r = md->stop;
  int cols = ml->nb_cols();
  T * res = md->thread_data;
  
  for (int r = start_r; r < stop_r; ++r)
    {
      for (int c1 = 0; c1 < cols; ++c1)
	{	  
	  T tmp = 0;
	  for (int c=0; c < cols; ++c) tmp += (*ml)(r,c) * (*mr)(c,c1);	        	  
	  res[(r-start_r) * cols + c1] = tmp;
	}
    }  
}

template<typename T>
void Matrix<T>::parallel_add(void * data)
{
  parallel_data<T> * md = (parallel_data<T>*) data;
  Matrix<T> * mr  = (Matrix<T>*) md->other;
  Matrix<T> * ml  = (Matrix<T>*) md->thism;

  assert(ml->nb_cols() == mr->nb_cols() && ml->nb_rows() == mr->nb_rows());
  
  int start_r = md->start;
  int stop_r = md->stop;
  int cols = ml->nb_cols();
  T * res = md->thread_data;
  
  for (int r = start_r; r < stop_r; ++r)
    {
      for (int c1 = 0; c1 < cols; ++c1)
	  res[(r-start_r) * cols + c1] = (*ml)(r,c1) + (*mr)(r,c1);
    }  
}

template<typename T>
void Matrix<T>::parallel_sub(void * data)
{
  parallel_data<T> * md = (parallel_data<T>*) data;
  Matrix<T> * mr  = (Matrix<T>*) md->other;
  Matrix<T> * ml  = (Matrix<T>*) md->thism;

  assert(ml->nb_cols() == mr->nb_cols() && ml->nb_rows() == mr->nb_rows());
  
  int start_r = md->start;
  int stop_r = md->stop;
  int cols = ml->nb_cols();
  T * res = md->thread_data;
  
  for (int r = start_r; r < stop_r; ++r)
    {
      for (int c1 = 0; c1 < cols; ++c1)
	  res[(r-start_r) * cols + c1] = (*ml)(r,c1) - (*mr)(r,c1);
    }  
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& other)
{
  assert(this->cols == other.rows);
  Matrix res(this->rows, other.cols);    
  Timer clock;    
  if (parallel)
    {
      std::vector<Thread*> threads;
      std::vector<parallel_data<T>*> thread_data_v;

      int remainder = this->rows % nb_threads;
      int size_per_thread = (this->rows / nb_threads) * other.cols;
      //printf("..rows(%d), cols(%d), threads(%d), size(%d), remainder(%d)\n",this->rows, other.cols, nb_threads, size_per_thread, remainder);
      
      int start =0;
      for (int j=0; j<nb_threads; ++j)
	{
	  parallel_data<T> * md = new parallel_data<T>();
	  md->other = &other;
	  md->thism = this;

	  if (j == nb_threads - 1) size_per_thread += remainder * other.cols;
	  md->thread_data = new T[size_per_thread];
	  md->start = start;
	  md->stop = md->start + rows / nb_threads;
	  md->arr_size = size_per_thread;
	  if(j == nb_threads - 1) md->stop += remainder;

	  //printf("..summary thread(%d)), start(%d), stop(%d), size(%d)\n",j, md->start, md->stop, md->arr_size);

	  start = md->stop;
	  threads.push_back(new Thread(j, (thread_fnc)parallel_multiply, (void*) md) );	  
	  thread_data_v.push_back(md);
	}
      clock.start();
      for (int j=0; j<nb_threads; ++j) threads[j]->start();
      for (int j=0; j<nb_threads; ++j) threads[j]->join();
      
      // TODO : make parallel also
      int offset = 0;
      for (int j=0; j<nb_threads; ++j)
	{
	  for (int i=0; i<thread_data_v[j]->arr_size; ++i) 
	    {
	      res.set(offset+i, thread_data_v[j]->thread_data[i]);
	    }
	  offset += thread_data_v[j]->arr_size;
	}
      clock.stop();
      clock.print("..parallel matrix multiply");
      
    }
  else
    {
      clock.start();
      for (int r = 0; r < rows; ++r)
	{
	  for (int c1 = 0; c1 < cols; ++c1)
	    {	  
	      T tmp = T(0);
	      for (int c=0; c < cols; ++c) tmp += (*this)(r,c) * other(c,c1);	        	  
	      res.set(r,c1,tmp);
	    }
	}
      clock.stop();
      clock.print("..serial matrix multiply");
    }
  /*
      int err = 0;
      for (int r=0; r < rows; ++r)
	{
	  for (int c1 = 0; c1 < cols; ++c1)
	    {
	      T tmp = res(r,c1) - res2(r,c1);
	      if (tmp.real != 0 || tmp.imag != 0 ) 
		{
		  printf("[%d,%d] --> (%g,%g) != (%g,%g)\n", r,c1, res(r,c1).real, res(r,c1).imag, res2(r,c1).real, res(r,c1).imag);
		  err++;
		}
	    }
	}
      printf("..scan finished, %d errors!\n", err);
  */
  return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& other)
{
  assert(this->cols == other.cols && this->rows == other.rows);
  Matrix res(this->rows, this->cols);    
  Timer clock;    
  if (parallel)
    {
      std::vector<Thread*> threads;
      std::vector<parallel_data<T>*> thread_data_v;
      parallel_data<T> * md = new parallel_data<T>();

      int remainder = this->rows % nb_threads;
      int size_per_thread = (this->rows / nb_threads) * other.cols;
      //printf("..rows(%d), cols(%d), threads(%d), size(%d), remainder(%d)\n",this->rows, other.cols, nb_threads, size_per_thread, remainder);

      int start =0;
      for (int j=0; j<nb_threads; ++j)
	{
	  parallel_data<T> * md = new parallel_data<T>();
	  md->other = &other;
	  md->thism = this;
	  
	  if (j == nb_threads - 1) size_per_thread += remainder * other.cols;
	  md->thread_data = new T[size_per_thread];
	  md->start = start;
	  md->stop = md->start + rows / nb_threads;
	  md->arr_size = size_per_thread;
	  if(j == nb_threads - 1) md->stop += remainder;
	  
	  //printf("..summary thread(%d)), start(%d), stop(%d), size(%d)\n",j, md->start, md->stop, md->arr_size);
	  
	  start = md->stop;
	  threads.push_back(new Thread(j, (thread_fnc)parallel_add, (void*) md) );	  
	  thread_data_v.push_back(md);
	}
      clock.start();
      for (int j=0; j<nb_threads; ++j) threads[j]->start();
      for (int j=0; j<nb_threads; ++j) threads[j]->join();
      
      // TODO : make parallel also
      int offset = 0;
      for (int j=0; j<nb_threads; ++j)
	{
	  for (int i=0; i<thread_data_v[j]->arr_size; ++i) 
	    {
	      res.set(offset+i, thread_data_v[j]->thread_data[i]);
	    }
	  offset += thread_data_v[j]->arr_size;
	}
      clock.stop();
      clock.print("..parallel matrix multiply");
      
    }
  else
    {
      clock.start();
      for (int r = 0; r < rows; ++r)
	{
	  for (int c1 = 0; c1 < cols; ++c1)
	    {		      
	      T tmp = (*this)(r,c1) + other(r,c1);
	      res.set(r,c1,tmp);
	    }
	}
      clock.stop();
      clock.print("..serial matrix multiply");
    }
  return res;
}


template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix& other)
{
  assert(this->cols == other.cols && this->rows == other.rows);
  Matrix res(this->rows, this->cols);    
  Timer clock;    
  if (parallel)
    {
      std::vector<Thread*> threads;
      std::vector<parallel_data<T>*> thread_data_v;
      parallel_data<T> * md = new parallel_data<T>();

      int remainder = this->rows % nb_threads;
      int size_per_thread = (this->rows / nb_threads) * other.cols;
      //printf("..rows(%d), cols(%d), threads(%d), size(%d), remainder(%d)\n",this->rows, other.cols, nb_threads, size_per_thread, remainder);

      int start =0;
      for (int j=0; j<nb_threads; ++j)
	{
	  parallel_data<T> * md = new parallel_data<T>();
	  md->other = &other;
	  md->thism = this;
	  
	  if (j == nb_threads - 1) size_per_thread += remainder * other.cols;
	  md->thread_data = new T[size_per_thread];
	  md->start = start;
	  md->stop = md->start + rows / nb_threads;
	  md->arr_size = size_per_thread;
	  if(j == nb_threads - 1) md->stop += remainder;
	  
	  //printf("..summary thread(%d)), start(%d), stop(%d), size(%d)\n",j, md->start, md->stop, md->arr_size);
	  
	  start = md->stop;
	  threads.push_back(new Thread(j, (thread_fnc)parallel_sub, (void*) md) );	  
	  thread_data_v.push_back(md);
	}
      clock.start();
      for (int j=0; j<nb_threads; ++j) threads[j]->start();
      for (int j=0; j<nb_threads; ++j) threads[j]->join();
      
      // TODO : make parallel also
      int offset = 0;
      for (int j=0; j<nb_threads; ++j)
	{
	  for (int i=0; i<thread_data_v[j]->arr_size; ++i) 
	    {
	      res.set(offset+i, thread_data_v[j]->thread_data[i]);
	    }
	  offset += thread_data_v[j]->arr_size;
	}
      clock.stop();
      clock.print("..parallel matrix multiply");
      
    }
  else
    {
      clock.start();
      for (int r = 0; r < rows; ++r)
	{
	  for (int c1 = 0; c1 < cols; ++c1)
	    {		      
	      T tmp = (*this)(r,c1) - other(r,c1);
	      res.set(r,c1,tmp);
	    }
	}
      clock.stop();
      clock.print("..serial matrix multiply");
    }
  return res;
}

// matrix resize/minor/column/row/transpose operations
template<typename T>
Matrix<T> Matrix<T>::minor(int r, int c)
{
  assert(r <= rows && c <= cols);
  if (rows == 0) rows = 1;
  if (cols == 0) cols = 1;
  //if (rows == cols == 1) return Matrix(0,0);

  Matrix res(rows-1,cols-1);
  for(int row=0; row<rows; ++row)
    {
      if (row == r) continue;
      for (int col=0; col<cols; ++col)
	{
	  if (col == c) continue;	  	  
	  res.set(row>=rows-1?row-1:row,col>=cols-1?col-1:col, data[row*cols+col]);
	}
    }
  return res; 
}

template<typename T>
Vector<T> Matrix<T>::column(int c) const
{
  assert(c <= cols);
  Vector<T> res(rows);
  for(int r=0; r<rows; ++r) res.set(r,(*this)(r,c));
  return res;
}

template<typename T>
void Matrix<T>::set_column(int c, const Vector<T>& vin) const
{
  assert(vin.nb_rows() == rows && c <= cols);
  for(int r=0; r < rows; ++r)
    data[r*cols + c] = vin(r);
}

template<typename T>
Matrix<T> Matrix<T>::transpose()
{
  Matrix trans(cols, rows);
  for(int r=0; r<rows; ++r)
    {
      for(int c=0; c<cols; ++c) trans.set(c,r, (*this)(r,c));
    }
  return trans;
}

template<typename T>
Matrix<T> Matrix<T>::identity()
{
  Matrix id(rows,cols);
  for(int j=0; j<rows*cols; j += cols + 1)
    id.set(j,T(1));

  return id;
}

// matrix-vector multiplication
// todo-parallel implementation
template<typename T>
Vector<T> Matrix<T>::operator*(const Vector<T>& other)
{
  Vector<T> res(other);
  for (int r = 0; r < rows; ++r)
    {
      T tmp = T(0);
      for (int c=0; c < cols; ++c) tmp += (*this)(r,c) * other(c);	        	  
      res.set(r,tmp);      
    }  
  return res;
}
#endif

