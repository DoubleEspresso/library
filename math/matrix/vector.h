#pragma once

#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H


template<typename T>
class Vector
{
public:
	Vector(int rws) : rows(rws), size(rws), parallel(false), nb_threads(1)
	{
		data = new T[size];
		memset(data, 0, sizeof(T)*size);
		parallel = (rows >= 275);
		if (parallel)
		{
			bool hyperthreading = false;
			nb_threads = Hardware::cpu_count(hyperthreading);
			nb_threads *= (hyperthreading ? 2 : 1);
		}
	}

	Vector<T>(const Vector& other)
	{
		data = 0;
		rows = other.rows;
		size = rows;
		data = new T[size];
		memcpy(this->data, other.data, sizeof(T) * size);
	}
	~Vector() { if (data) { delete[] data; data = 0; } }

	T operator()(int r) const { return data[r]; }
	void operator=(const Vector& other) { memcpy(this->data, other.data, sizeof(T) * size); }
	void set(int j, T val) { data[j] = T(val); }
	int nb_rows() { return rows; }
	int nb_cols() { return rows; }
	void set_threads(int j) { nb_threads = j; }
	int get_threads() { return nb_threads; }
	T dot(const Vector& other);
	Vector cross(const Vector& other);

	// vector-vector operations
	Vector operator+(const Vector& other);
	Vector operator-(const Vector& other);
	Vector operator+=(const Vector& other);
	Vector operator-=(const Vector& other);
	Vector operator*=(const T& other);
	Vector operator*(const T& other);

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


#endif
