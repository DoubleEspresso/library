#include "union.h"

UnionFindBase::UnionFindBase(UINT n) : N(n), sites(0), size(0)
{
  if(!sites && (sites = new UINT[N]) == NULL)
    {
      printf("..failed to alloc %d sites\n",N);
      return;
    }
  if ( (size = new UINT[N]) == NULL )
    {
      printf("..failed to alloc %d sizes arr\n",N);
      return;
    }
  else
    {
      for (UINT j=0; j<N; ++j) sites[j] = j;
      for (UINT j=0; j<N; ++j) size[j] = 1;
    }
}

UnionFindBase::~UnionFindBase()
{
  if (sites) {delete [] sites; sites = 0;}
  if (size) {delete [] size; size = 0;}
}


// -------------quick find------------------------
bool QuickFind::Connected(UINT p, UINT q)
{
  if (!isOk(p,q)) return false;
  else return id(p) == id(q);
}

void QuickFind::Union(UINT p, UINT q)
{
  if (!isOk(p,q)) return;
  else if (id(p) == id(q)) return;
  UINT sp = id(p);
  UINT sq = id(q);
  for (UINT j=0; j<count(); ++j)
    {
      if (id(j) == sp) set(j, sq); 
    }
}


// -------------union find------------------------
UINT QuickUnion::root(UINT p)
{
  while (p != id(p)) 
    {
      set(p,  id( id(p) ) ); // path compression
      p = id(p);
    }
  return p;
}

bool QuickUnion::Connected(UINT p, UINT q)
{
  //printf("CONNECTED: root of %d = %d, root of %d = %d\n",p,root(p),q,root(q));
  if (!isOk(p,q)) return false;
  return root(p) == root(q);
}

void QuickUnion::Union(UINT p, UINT q)
{
  if (!isOk(p,q)) return;
  UINT rp = root(p);
  UINT rq = root(q);
  //printf("UNION: root of %d = %d, root of %d = %d\n",p,root(p),q,root(q));
  if (rp == rq) return;

  // keep trees flat-ish
  if (sz(rp) < sz(rq)) { set(rp, rq); add(rq, sz(rp)); }
  else { set(rq,rp); add(rp, sz(rq)); } 
  //printf("after union: root of %d = %d, root of %d = %d\n",p,root(p),q,root(q));
}

