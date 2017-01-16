#pragma once

#ifndef MATH_MATRIX_H
#define MATH_MATRIX_H

#include <cassert>
#include <string.h>
#include <vector>
#include <math.h>
#include <algorithm>    
#include <array>        
#include <random>       
#include <chrono>       

#include "../../system/threads.h"
#include "../../utils/timer.h"
#include "../../system/hardware.h"
#include "../../utils/array_utils.h"

template <typename T> class Matrix;

template<typename T>
class Vector
{
public:
	Vector(int rws) : rows(rws), size(rws), parallel(false), nb_threads(1)
	{
		data = new T[size];
		memset(data, 0, sizeof(T)*size);
		parallel = false;// (rows >= 275);
		if (parallel)
		{
			bool hyperthreading = false;
			nb_threads = Hardware::cpu_count(hyperthreading);
			nb_threads *= (hyperthreading ? 2 : 1);
		}
	}

	Vector(const Matrix<T>& other, int r, int c, bool is_col = true)
	{
		data = 0;
		rows = (is_col ? other.nb_rows() : other.nb_cols());
		size = rows;
		data = new T[size];
		for (int j = 0; j < size; ++j) data[j] = other(r + (is_col ? j : 0), c + (!is_col ? j : 0));
	}
	Vector<T>(const Vector& other)
	{
		data = 0;
		rows = other.rows;
		size = rows;
		data = new T[size];
		memcpy(this->data, other.get_data(), sizeof(T) * size);
	}
	~Vector() { if (data) { delete[] data; data = 0; } }
	T * get_data() const { return data; }
	T operator()(int r) const { return data[r]; }
	Vector<T>& operator=(const Vector& other) 
	{
		if (data) { delete[] data; data = 0; }
		memcpy(this->data, other.get_data(), sizeof(T) * size);  
		rows = other.nb_rows();
		cols = other.nb_cols();
		size = rows * cols;
		return *this;
	}
	void set(int j, T val) { data[j] = T(val); }
	int nb_rows() const { return rows; }
	int nb_cols() const { return rows; }
	void set_threads(int j) { nb_threads = j; }
	int get_threads() { return nb_threads; }
	T dot(const Vector& other);
	T norm();
	Vector normalize();
	Vector cross(const Vector& other);
	Vector conj() const;

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
	assert(this->rows == other.nb_rows());
	Vector<T> cc(other.conj());
	T res(0);
	for (int j = 0; j < this->rows; ++j) res += (data[j] * cc(j));

	return res;
}

template<typename T>
Vector<T> Vector<T>::operator/(const T& other)
{
	Vector res(*this);
	for (int j = 0; j < this->rows; ++j) res.set(j, data[j] / other);
	return res;
}

template<typename T>
Vector<T> Vector<T>::operator*(const T& other)
{
	Vector res(*this);
	for (int j = 0; j < rows; ++j) res.set(j, data[j] * other);
	return res;
}

template<typename T>
Vector<T> Vector<T>::operator-(const Vector& other)
{
	Vector res(*this);
	for (int j = 0; j < rows; ++j) res.set(j, data[j] - other(j));
	return res;
}

template<typename T>
inline Vector<T> operator*(const T& other, Vector<T>& vec)
{
	return vec * other;
}

// vector conjugation for complex types
template<typename T>
inline Vector<T> Vector<T>::conj() const
{
	Vector res(*this);
	for (int j = 0; j < rows; ++j) res.set(j, data[j].conj());
	return res;
}

template<>
inline Vector<int> Vector<int>::conj() const
{
	return *this;
}

template<>
inline Vector<float> Vector<float>::conj() const
{
	return *this;
}

template<>
inline Vector<double> Vector<double>::conj() const
{
	return *this;
}


template<typename T>
T Vector<T>::norm()
{
	T tmp = (*this).dot(*this);
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
	Vector res(*this);
	return res / res.norm();
}

template<typename T>
class Matrix
{
public:
	//Matrix() {};
	Matrix(const Vector<T>& other)
	{
		size = other.nb_rows();
		cols = 0;
		rows = other.nb_rows();
		data = new T[size];
		memset(data, 0, sizeof(T)*size);
		parallel = false;// (size >= 275);
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
		memset(data, 0, sizeof(T)*size);
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
		memcpy(this->data, other.get_data(), sizeof(T) * size);
	}

	~Matrix() 
	{ 
		if (data) 
		{ 
			size_t s = size;

			delete[] data; 
			data = 0; 
		} 
	}

	T operator()(int r, int c) const { return data[r*cols + c]; }

	Matrix<T>& operator=(const Matrix& other) 
	{
		if (rows != other.nb_rows() || cols != other.nb_cols())
		{
			free();
			data = new T[other.Size()];
		}
		memcpy(data, other.get_data(), sizeof(T) * other.Size());
		rows = other.nb_rows();
		cols = other.nb_cols();
		size = rows * cols;
		return *this;
	}
	void free()
	{
		if (data) { delete[] data; data = 0; }
	}
	void scan() { parallel = false; } // (size >= 275 * 275); } // for example
	void clear() { if (data) memset(data, 0, sizeof(T) * size); }
	Matrix identity();
	bool is_square() { return rows == cols; }
	void set_threads(int j) { nb_threads = j; }
	int get_threads() { return nb_threads; }
	Matrix minor(int r, int c);
	Matrix transpose();
	Matrix conj();
	void lower_triangle(Matrix& D);
	void upper_triangle(Matrix& D);

	// access
	T data_at(int j) { return T(data[j]); }
	T data_at(int r, int c) { return T(data[r * cols + c]); }
	void set(int r, int c, T val) { data[r*cols + c] = T(val); }
	void set(int j, T val) { data[j] = val; }
	void set(const Matrix<T>& other)
	{
		free();
		this->data = new T[other.Size()];
		memcpy(this->data, other.get_data(), sizeof(T) * size);
		rows = other.nb_rows();
		cols = other.nb_cols();
		size = rows * cols;
	}
	int nb_rows() const { return rows; }
	int nb_cols() const { return cols; }
	T * get_data() const { return data; }

	// utility
	Vector<T> get_col(int c) const;
	Matrix<T> * get_cols(int r, int n) const;
	Vector<T> get_row(int r) const;
	Matrix<T> * get_rows(int r, int n) const;
	void set_col(int c, const Vector<T>& vin) const;
	void set_row(int r, const Vector<T>& rin) const;
	//void reshape(int r, int c) const;
	void shuffle_rows(Matrix& storage);
	void shuffle_cols(Matrix& storage);
	void pad(int i); // inserts row/cols of identity into upper left portion of matrix (in place)
	void submatrix(int r, int c); // in place return of matrix starting at idx (r,c)
	size_t Size() const { return rows*cols; }

	// matrix-matrix operations
	Matrix operator+(const Matrix& other);
	Matrix operator-(const Matrix& other);
	Matrix operator*(const Matrix& other);
	Matrix operator+=(const Matrix& other);
	Matrix operator-=(const Matrix& other);
	Matrix operator*=(const Matrix& other);
	Matrix operator*(const T& other) const;
	//Matrix operator*(const T& other, const Matrix<T>& m) const;

	// matrix-vector operations
	Vector<T> operator*(const Vector<T>& other);

	// parallel operator overloads
	static void parallel_multiply(void * data);
	static void parallel_add(void * data);
	static void parallel_sub(void * data);

	// debug utilities
	void print(std::string Label = "")
	{
		if (Label != "") printf(" == %s ==\n", Label.c_str());
		for (int r = 0; r < rows; ++r)
		{
			for (int c = 0; c < cols; ++c)
			{
				printf("%.3f ", data[r*cols + c]);
			}
			printf("\n");
		}
	}

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
	Matrix<T> * mr = (Matrix<T>*) md->other;
	Matrix<T> * ml = (Matrix<T>*) md->thism;

	assert(ml->nb_cols() == mr->nb_rows());

	int start_r = md->start;
	int stop_r = md->stop;
	int cols = mr->nb_cols();
	T * res = md->thread_data;

	for (int r = start_r; r < stop_r; ++r)
	{
		for (int c1 = 0; c1 < cols; ++c1)
		{
			T tmp = 0;
			for (int c = 0; c < cols; ++c) tmp += (*ml)(r, c) * (*mr)(c, c1);
			res[(r - start_r) * cols + c1] = tmp;
		}
	}
}

template<typename T>
void Matrix<T>::parallel_add(void * data)
{
	parallel_data<T> * md = (parallel_data<T>*) data;
	Matrix<T> * mr = (Matrix<T>*) md->other;
	Matrix<T> * ml = (Matrix<T>*) md->thism;

	assert(ml->nb_cols() == mr->nb_cols() && ml->nb_rows() == mr->nb_rows());

	int start_r = md->start;
	int stop_r = md->stop;
	int cols = ml->nb_cols();
	T * res = md->thread_data;

	for (int r = start_r; r < stop_r; ++r)
	{
		for (int c1 = 0; c1 < cols; ++c1)
			res[(r - start_r) * cols + c1] = (*ml)(r, c1) + (*mr)(r, c1);
	}
}

template<typename T>
void Matrix<T>::parallel_sub(void * data)
{
	parallel_data<T> * md = (parallel_data<T>*) data;
	Matrix<T> * mr = (Matrix<T>*) md->other;
	Matrix<T> * ml = (Matrix<T>*) md->thism;

	assert(ml->nb_cols() == mr->nb_cols() && ml->nb_rows() == mr->nb_rows());

	int start_r = md->start;
	int stop_r = md->stop;
	int cols = ml->nb_cols();
	T * res = md->thread_data;

	for (int r = start_r; r < stop_r; ++r)
	{
		for (int c1 = 0; c1 < cols; ++c1)
			res[(r - start_r) * cols + c1] = (*ml)(r, c1) - (*mr)(r, c1);
	}
}

template<typename T>
void Matrix<T>::lower_triangle(Matrix& D)
{
	for (int r = 0; r < rows; ++r)
	{
		for (int c = 0; c < r; ++c)
		{
			D.set(r, c, data[r * cols + c]);
		}
	}
}

template<typename T>
void Matrix<T>::upper_triangle(Matrix& D)
{
	for (int r = 0; r < rows; ++r)
	{
		for (int c = r; c < cols; ++c)
		{
			D.set(r, c, data[r * cols + c]);
		}
	}
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& other)
{
	assert(this->cols == other.rows);
	Matrix res(this->rows, other.cols);
	Timer clock;
	if (false)
	{
		//std::vector<Thread*> threads;
		THREAD_HANDLE * threads = new THREAD_HANDLE[nb_threads];
		std::vector<parallel_data<T>*> thread_data;

		int remainder = this->rows % nb_threads;
		int size_per_thread = (this->rows / nb_threads) * other.cols;
		//printf("..rows(%d), cols(%d), threads(%d), size(%d), remainder(%d)\n",this->rows, other.cols, nb_threads, size_per_thread, remainder);

		int start = 0;
		for (int j = 0; j < nb_threads; ++j)
		{
			parallel_data<T> * md = new parallel_data<T>();
			md->other = &other;
			md->thism = this;

			if (j == nb_threads - 1) size_per_thread += remainder * other.cols;
			md->thread_data = new T[size_per_thread];
			md->start = start;
			md->stop = md->start + rows / nb_threads;
			md->arr_size = size_per_thread;
			if (j == nb_threads - 1) md->stop += remainder;

			//printf("..summary thread(%d)), start(%d), stop(%d), size(%d)\n",j, md->start, md->stop, md->arr_size);

			start = md->stop;
			thread_data.push_back(md);
		}
		clock.start();
		for (int j = 0; j < nb_threads; ++j) threads[j] = start_thread((thread_fnc)parallel_multiply, (void*)thread_data[j]);
		wait_threads_finish(threads, nb_threads);

		// TODO : make parallel also
		int offset = 0;
		for (int j = 0; j < nb_threads; ++j)
		{
			for (int i = 0; i < thread_data[j]->arr_size; ++i)
			{
				res.set(offset + i, thread_data[j]->thread_data[i]);
			}
			offset += thread_data[j]->arr_size;
		}
		clock.stop();
		//clock.print("..parallel matrix multiply");
	}
	else
	{
		clock.start();
		for (int r = 0; r < rows; ++r)
		{
			for (int c1 = 0; c1 < other.nb_cols(); ++c1)
			{
				T tmp = T(0);
				for (int c = 0; c < cols; ++c) tmp += (*this)(r, c) * other(c, c1);
				res.set(r, c1, tmp);
			}
		}
		clock.stop();
		//clock.print("..serial matrix multiply");
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
	if (false)//(parallel)
	{
		THREAD_HANDLE* threads = new THREAD_HANDLE[nb_threads];
		std::vector<parallel_data<T>*> thread_data;
		parallel_data<T> * md = new parallel_data<T>();

		int remainder = this->rows % nb_threads;
		int size_per_thread = (this->rows / nb_threads) * other.cols;
		//printf("..rows(%d), cols(%d), threads(%d), size(%d), remainder(%d)\n",this->rows, other.cols, nb_threads, size_per_thread, remainder);

		int start = 0;
		for (int j = 0; j < nb_threads; ++j)
		{
			parallel_data<T> * md = new parallel_data<T>();
			md->other = &other;
			md->thism = this;

			if (j == nb_threads - 1) size_per_thread += remainder * other.cols;
			md->thread_data = new T[size_per_thread];
			md->start = start;
			md->stop = md->start + rows / nb_threads;
			md->arr_size = size_per_thread;
			if (j == nb_threads - 1) md->stop += remainder;

			//printf("..summary thread(%d)), start(%d), stop(%d), size(%d)\n",j, md->start, md->stop, md->arr_size);

			start = md->stop;
			thread_data.push_back(md);
		}
		clock.start();
		for (int j = 0; j < nb_threads; ++j) threads[j] = start_thread((thread_fnc)parallel_add, (void*)thread_data[j]);
		wait_threads_finish(threads, nb_threads);

		// TODO : make parallel also
		int offset = 0;
		for (int j = 0; j < nb_threads; ++j)
		{
			for (int i = 0; i < thread_data[j]->arr_size; ++i)
			{
				res.set(offset + i, thread_data[j]->thread_data[i]);
			}
			offset += thread_data[j]->arr_size;
		}
		clock.stop();
		//clock.print("..parallel matrix multiply");

	}
	else
	{
		clock.start();
		for (int r = 0; r < rows; ++r)
		{
			for (int c1 = 0; c1 < cols; ++c1)
			{
				T tmp = (*this)(r, c1) + other(r, c1);
				res.set(r, c1, tmp);
			}
		}
		clock.stop();
		//clock.print("..serial matrix multiply");
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
		THREAD_HANDLE *  threads = new THREAD_HANDLE[nb_threads];
		std::vector<parallel_data<T>*> thread_data;
		parallel_data<T> * md = new parallel_data<T>();

		int remainder = this->rows % nb_threads;
		int size_per_thread = (this->rows / nb_threads) * other.cols;
		//printf("..rows(%d), cols(%d), threads(%d), size(%d), remainder(%d)\n",this->rows, other.cols, nb_threads, size_per_thread, remainder);

		int start = 0;
		for (int j = 0; j < nb_threads; ++j)
		{
			parallel_data<T> * md = new parallel_data<T>();
			md->other = &other;
			md->thism = this;

			if (j == nb_threads - 1) size_per_thread += remainder * other.cols;
			md->thread_data = new T[size_per_thread];
			md->start = start;
			md->stop = md->start + rows / nb_threads;
			md->arr_size = size_per_thread;
			if (j == nb_threads - 1) md->stop += remainder;

			//printf("..summary thread(%d)), start(%d), stop(%d), size(%d)\n",j, md->start, md->stop, md->arr_size);

			start = md->stop;
			thread_data.push_back(md);
		}
		clock.start();
		for (int j = 0; j < nb_threads; ++j) threads[j] = start_thread((thread_fnc)parallel_sub, (void*)thread_data[j]);
		wait_threads_finish(threads, nb_threads);

		// TODO : make parallel also
		int offset = 0;
		for (int j = 0; j < nb_threads; ++j)
		{
			for (int i = 0; i < thread_data[j]->arr_size; ++i)
			{
				res.set(offset + i, thread_data[j]->thread_data[i]);
			}
			offset += thread_data[j]->arr_size;
		}
		clock.stop();
		//clock.print("..parallel matrix multiply");

	}
	else
	{
		clock.start();
		for (int r = 0; r < rows; ++r)
		{
			for (int c1 = 0; c1 < cols; ++c1)
			{
				T tmp = (*this)(r, c1) - other(r, c1);
				res.set(r, c1, tmp);
			}
		}
		clock.stop();
		//clock.print("..serial matrix subtract");
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

	Matrix res(rows - 1, cols - 1);
	for (int row = 0; row < rows; ++row)
	{
		if (row == r) continue;
		for (int col = 0; col < cols; ++col)
		{
			if (col == c) continue;
			res.set(row >= rows - 1 ? row - 1 : row, col >= cols - 1 ? col - 1 : col, data[row*cols + col]);
		}
	}
	return res;
}

template<typename T>
Vector<T> Matrix<T>::get_row(int r) const
{
	assert(r <= rows);
	Vector<T> res(cols);
	for (int c = 0; c < cols; ++c) res.set(c, (*this)(r, c));
	return res;
}

// returns n x cols submatrix of n-rows, starting at row-index r.
template<typename T>
Matrix<T> * Matrix<T>::get_rows(int r, int n) const
{
	assert(r <= rows);
	Matrix<T> * res = new Matrix<T>(n, cols);
	for (int s = r, i = 0; s < r + n; ++s, ++i) res->set_row(i, get_row(s));
	return res;
}

template<typename T>
Vector<T> Matrix<T>::get_col(int c) const
{
	assert(c <= cols);
	Vector<T> res(rows);
	for (int r = 0; r < rows; ++r) res.set(r, (*this)(r, c));
	return res;
}

// returns n x cols submatrix of n-rows, starting at row-index r.
template<typename T>
Matrix<T> * Matrix<T>::get_cols(int c, int n) const
{
	assert(c <= cols);
	Matrix<T> * res = new Matrix<T>(rows, n);
	for (int s = c, i = 0; s < c + n; ++c, ++i) res->set_col(i, get_col(s));
	return res;
}

template<typename T>
void Matrix<T>::set_col(int c, const Vector<T>& vin) const
{
	assert(vin.nb_rows() == rows && c <= cols);
	for (int r = 0; r < rows; ++r)
		data[r*cols + c] = vin(r);
}

template<typename T>
void Matrix<T>::set_row(int r, const Vector<T>& rin) const
{
	assert(rin.nb_cols() == cols && r <= rows);
	for (int c = 0; c < cols; ++c)
		data[r * cols + c] = rin(r);
}

template<typename T>
Matrix<T> Matrix<T>::conj()
{
	Matrix cc(rows, cols);
	for (int r = 0; r < rows; ++r)
	{
		for (int c = 0; c < cols; ++c)
		{
			T tmp = data[r*cols + c]; tmp.imag *= -1;
			cc.set(r, c, tmp);
		}
	}
	return cc;
}

template<typename T>
Matrix<T> Matrix<T>::transpose()
{
	Matrix trans(cols, rows);
	for (int r = 0; r < rows; ++r)
	{
		for (int c = 0; c < cols; ++c) trans.set(c, r, (*this)(r, c));
	}
	return trans;
}

template<typename T>
Matrix<T> Matrix<T>::identity()
{
	Matrix id(rows, cols);
	for (int j = 0; j < rows*cols; j += cols + 1)
		id.set(j, T(1));

	return id;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& other) const
{
	Matrix res(*this);
	for (int j = 0; j < rows*cols; ++j) res.set(j, data[j] * other);
	return res;
}

//template<typename T>
Matrix<float> operator*(const float& other, const Matrix<float>& m)
{
	return m * other;
}

// inserts "i" identity rows (and corresponding cols)
// into the upper left portion of matrix, this is an inplace
// operation, which will resize the data array
template<typename T>
void Matrix<T>::pad(int i)
{
	int nrows = rows + i;
	int ncols = cols + i;
	T * new_data = new T[nrows * ncols];
	memset(new_data, T(0), nrows * ncols * sizeof(T));

	for (int r = 0; r < nrows; ++r)
	{
		for (int c = 0; c < ncols; ++c)
		{
			new_data[r * ncols + c] = (r == c && c < i) ? T(1) :
				(r >= i && c >= i) ? data[(r - i) * cols + (c - i)] : T(0);
		}
	}
	delete[] data; data = 0;
	rows = nrows;
	cols = ncols;
	size = nrows * ncols;
	data = new_data;
}

template <typename T>
void Matrix<T>::shuffle_cols(Matrix<T>& storage)
{
	int * indices = new int[cols];
	for (int j = 0; j < cols; ++j) indices[j] = j;

	shuffle<int>(indices, cols);

	for (int j = 0; j < cols; ++j) storage.set_col(indices[j], get_col(j));
}

template <typename T>
void Matrix<T>::shuffle_rows(Matrix<T>& storage)
{
	int * indices = new int[rows];
	for (int j = 0; j < rows; ++j) indices[j] = j;

	shuffle<int>(indices, rows);

	for (int j = 0; j < rows; ++j) storage.set_row(indices[j], get_row(j));
}

// transform (in place) current data to subset of matrix data
// which starts at index (r,c) .. deprecated (do not use)
template <typename T>
void Matrix<T>::submatrix(int r, int c)
{
	assert(r < rows && c < cols);
	int nrows = rows - r;
	int ncols = cols - c;
	T * new_data = new T[nrows * ncols];
	memset(new_data, T(0), nrows * ncols * sizeof(T));
	for (int row = 0; row < rows; ++row)
	{
		for (int col = 0; col < cols; ++col)
		{
			if (row >= r && col >= c) new_data[(row - r) * ncols + (col - c)] = data[row * cols + col];
		}
	}
	delete[] data; data = 0;
	rows = nrows;
	cols = ncols;
	size = nrows * ncols;
	data = new_data;
};

// matrix-vector multiplication
// todo-parallel implementation
template<typename T>
Vector<T> Matrix<T>::operator*(const Vector<T>& other)
{
	Vector<T> res(other);
	for (int r = 0; r < rows; ++r)
	{
		T tmp = T(0);
		for (int c = 0; c < cols; ++c) tmp += (*this)(r, c) * other(c);
		res.set(r, tmp);
	}
	return res;
}
#endif

