// basic implementation of union-find
// quick union
// weighted quick union with path compression etc.

#ifndef ALGORITHM_UNION_H
#define ALGORITHM_UNION_H

#include <stdio.h>
#include <stdlib.h>
#include "../system/types.h"

class UnionFindBase
{
	uint N;
	uint * sites;
	uint * size;
public:
	UnionFindBase(uint n);
	virtual ~UnionFindBase();

	virtual void Union(uint p, uint q) = 0;
	virtual bool Connected(uint p, uint q) = 0;
	bool isOk(uint p, uint q) { return (p < N && q < N); }
	uint count() { return N; }
	uint id(uint i) { return sites[i]; }
	void set(uint i, uint val) { sites[i] = val; }
	uint sz(uint i) { return size[i]; }
	void add(uint i, uint j) { size[i] += j; }
};

class QuickFind : public UnionFindBase
{
	explicit QuickFind(uint n) : UnionFindBase(n) {}
	virtual void Union(uint p, uint q);
	virtual bool Connected(uint p, uint q);
};

class QuickUnion : public UnionFindBase
{
	explicit QuickUnion(uint n) : UnionFindBase(n) {}
	virtual void Union(uint p, uint q);
	virtual bool Connected(uint p, uint q);
	uint root(uint p);
};

#endif
