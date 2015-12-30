// basic implementation of union-find
// quick union
// weighted quick union with path compression etc.

#ifndef UNION_H_
#define UNION_H_

#include <stdio.h>
#include <stdlib.h>

typedef unsigned int UINT;

class UnionFindBase
{
 public:  
  UnionFindBase(UINT n);
  virtual ~UnionFindBase();

  // methods
  virtual void Union(UINT p, UINT q) =0;
  virtual bool Connected(UINT p, UINT q)=0;
  bool isOk(UINT p, UINT q) { return (p < N && q < N); }
  UINT count() { return N;}
  UINT id(UINT i) { return sites[i]; }
  void set(UINT i, UINT val) { sites[i] = val; }
  UINT sz(UINT i) { return size[i]; }
  void add(UINT i, UINT j) { size[i] += j; }

 private:
  UINT N; // nb of sites
  UINT * sites;
  UINT * size;
};

class QuickFind : public UnionFindBase
{
  explicit QuickFind(uint n) : UnionFindBase(n) {}
  virtual void Union(UINT p, UINT q);
  virtual bool Connected(UINT p, UINT q);
};

class QuickUnion : public UnionFindBase
{
 public: 
  explicit QuickUnion(uint n) : UnionFindBase(n) {}
  virtual void Union(UINT p, UINT q);
  virtual bool Connected(UINT p, UINT q);
  UINT root(UINT p);
};

#endif
