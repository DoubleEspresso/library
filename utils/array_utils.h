#ifndef UTILS_ARRAY_H
#define UTILS_ARRAY_H

#include <cstdlib>
#include <iostream>
#include <ctime>

template<typename T>
inline void shuffle(T * data, size_t sz)
{
	std::srand(std::time(0));
	for (int i = sz - 1; i >= 0; --i)
	{
		int j = i * std::rand() / RAND_MAX;
		swap(data, i, j);
	}
}

template<typename T>
inline void swap(T * data, int i, int j)
{
	T v = data[i]; data[i] = data[j]; data[j] = v;
}

#endif

