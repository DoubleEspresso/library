#ifndef UTILS_BITOPS_H
#define UTILS_BITOPS_H

#include "../system/types.h"

// SWAR popcount algorithm for U64 datatypes
inline int count64(uint64 b)
{
	if (!b) return 0;
	b = b - ((b >> 1) & 0x5555555555555555ULL);
	b = (b & 0x3333333333333333ULL) + ((b >> 2) & 0x3333333333333333ULL);
	b = (b + (b >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
	return (int)((b * 0x0101010101010101ULL) >> 56);
}

// SWAR popcount algorithm for U32 datatypes
inline int count32(uint64 b)
{
	if (!b) return 0;
	uint32 bl = uint32(b ^ (b << 32)); // lower 32 bits
	uint32 bu = uint32(b ^ (b >> 32)); // upper 32 bits

	int low_count = 0;
	int high_count = 0;

	if (bl)
	{
		bl = bl - ((bl >> 1) & 0x55555555);
		bl = (bl & 0x33333333) + ((bl >> 2) & 0x33333333);
		low_count = (int)((((bl + (bl >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24);
	}
	if (bu)
	{
		bu = bu - ((bu >> 1) & 0x55555555);
		bu = (bu & 0x33333333) + ((bu >> 2) & 0x33333333);
		high_count = (int)((((bu + (bu >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24);
	}
	return (int)(low_count + high_count);
}


inline int count(uint64 b)
{
	return (_64BIT ? count64(b) : count32(b));
}

// use when most bits in "b" are 0.

int count64_max15(uint64 b)
{
int count;
for (count = 0; b; count++)
b &= b - 1;
return count;
}

int count32_max15(uint32 b)
{
int count;
for (count = 0; b; count++)
b &= b - 1;
return count;
}


const uint64 debruijn64 = const uint64(0x03f79d71b4cb0a89);
const uint32 debruijn32 = const uint32(0x077CB531U);

const int idx64[64] =
{
	0,  1, 48,  2, 57, 49, 28,  3,
	61, 58, 50, 42, 38, 29, 17,  4,
	62, 55, 59, 36, 53, 51, 43, 22,
	45, 39, 33, 30, 24, 18, 12,  5,
	63, 47, 56, 27, 60, 41, 37, 16,
	54, 35, 52, 21, 44, 32, 23, 11,
	46, 26, 40, 15, 34, 20, 31, 10,
	25, 14, 19,  9, 13,  8,  7,  6
};

const int idx32[64] =
{
	0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
	31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9,

	32, 33, 50, 34, 61, 46, 56, 35, 62, 54, 52, 47, 57, 49, 36, 40,
	63, 59, 45, 55, 53, 51, 48, 39, 58, 44, 50, 38, 43, 37, 42, 41
};

// returns the LSB of the bitmap using lookup tables
inline int lsb64(uint64 b)
{
	return idx64[(int)(((b & (~b + 1ULL)) * debruijn64) >> 58)];
}

// lsb for 32-bit systems (overly slow ?).
inline int lsb32(uint64 b)
{
	uint32 bl = uint32(b ^ (b << 32)); // lower 32 bits
	uint32 bu = uint32(b ^ (b >> 32)); // upper 32 bits

	uint32 bn = (bl ? bl : bu);
	int idx = (bl ? 0 : 32);

	return idx32[(int)((((bn & (~bn + 1ULL)) * debruijn32) >> 27) + idx)];
}

inline int lsb(uint64 b)
{
	return (_64BIT ? lsb64(b) : lsb32(b));
};

inline int pop_lsb(uint64& b)
{
	const int s = (_64BIT ? lsb64(b) : lsb32(b));
	b &= (b - 1);
	return s;
};

#endif
