#include "union.h"

UnionFindBase::UnionFindBase(uint n) : N(n), sites(0), size(0)
{
  if(!sites && (sites = new uint[N]) == NULL)
    {
      printf("..failed to alloc %d sites\n",N);
      return;
    }
  if ( (size = new uint[N]) == NULL )
    {
      printf("..failed to alloc %d sizes arr\n",N);
      return;
    }
  else
    {
      for (uint j=0; j<N; ++j) sites[j] = j;
      for (uint j=0; j<N; ++j) size[j] = 1;
    }
}
UnionFindBase::~UnionFindBase()
{
  if (sites) {delete [] sites; sites = 0;}
  if (size) {delete [] size; size = 0;}
}

// quick find
bool QuickFind::Connected(uint p, uint q)
{
  if (!isOk(p,q)) return false;
  else return id(p) == id(q);
}
void QuickFind::Union(uint p, uint q)
{
  if (!isOk(p,q)) return;
  else if (id(p) == id(q)) return;
  uint sp = id(p);
  uint sq = id(q);
  for (uint j=0; j<count(); ++j)
    {
      if (id(j) == sp) set(j, sq); 
    }
}

// quick union
uint QuickUnion::root(uint p)
{
  while (p != id(p)) 
    {
      set(p,  id( id(p) ) ); // path compression
      p = id(p);
    }
  return p;
}
bool QuickUnion::Connected(uint p, uint q)
{
  //printf("CONNECTED: root of %d = %d, root of %d = %d\n",p,root(p),q,root(q));
  if (!isOk(p,q)) return false;
  return root(p) == root(q);
}
void QuickUnion::Union(uint p, uint q)
{
  if (!isOk(p,q)) return;
  uint rp = root(p);
  uint rq = root(q);
  //printf("UNION: root of %d = %d, root of %d = %d\n",p,root(p),q,root(q));
  if (rp == rq) return;

  // keep trees flat-ish
  if (sz(rp) < sz(rq)) { set(rp, rq); add(rq, sz(rp)); }
  else { set(rq,rp); add(rp, sz(rq)); } 
  //printf("after union: root of %d = %d, root of %d = %d\n",p,root(p),q,root(q));
}

