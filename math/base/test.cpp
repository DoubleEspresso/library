#include <stdio.h>
#include <cmath>
#include <vector>

#include "complex.h"
#include "../matrix/matrix.h"
#include "../matrix/linalg.h"
#include "../../utils/timer.h"


struct int_struct
{  
  int id;
  int * arr;
};

void add(void * data) 
{
  if (data == NULL) 
    {
      printf("work-fnc null\n");
      return;
    }
  int_struct * is = (int_struct*) data;
  printf("%d begins work\n",is->id);
  int start = is->id * 200000000/4;
  int end = start + 200000000/4;
  for(int j=start; j<end; ++j) is->arr[j] = exp(sin(is->arr[j] * is->arr[j]));
}

int main(int argc, char ** argv)
{
  // initialize
  /*
  int sz = 200000000;
  int * data = new int[200000000];
  for (int j =0; j<sz; ++j) data[j] = j;
 
  Timer clock;
  std::vector<Thread*> threads;
  for (int j=0; j<4; ++j)
    {
      int_struct * is = new int_struct();
      is->id = j;
      is->arr = data;
      threads.push_back(new Thread(j, (thread_fnc)add, (void*) is));
    }
  clock.start();

  // start
  for (int j=0; j<4; ++j) threads[j]->start();

  // join
  for (int j=0; j<4; ++j) threads[j]->join();

  clock.stop();
  clock.print("Parallel Addition");  


  for (int j =0; j<sz; ++j) data[j] = j;

  clock.start();
  for(int j=0 ;j<sz; ++j) data[j] = exp(sin(data[j] * data[j])); 
  clock.stop();
  clock.print("Serial Addition");
  */

  //bool hyperthreading = false;
  //int nb_cpus = Hardware::cpu_count(hyperthreading);
  //printf("..dbg physical cpus = %d, hyperthreading enabled = %s\n",nb_cpus, hyperthreading ? "true" : "false");


  Vector<Complex_f> v(3);

  Matrix<Complex_f> m1(3,3); // should compute in parallel automatically
  
  for (int j=0; j<3*3; ++j) m1.set(j, Complex_f(j+1,1)); 
  for (int j=0; j<3; ++j) v.set(j, Complex_f(j+1,1)); 

// what about real data types .. norm is different..
  Matrix<Complex_f> m2 = m1.minor(2,2);

  printf("..norm = %g,%g\n", v.norm().real, v.norm().imag);
  Vector<Complex_f> m3 = v.normalize();
  
  for(int r = 0; r < 3; ++r)
    {
      for (int c=0; c < 3; ++c)
	{
	    printf(" (%g,%g) ", m1(r,c).real, m1(r,c).imag);
	}
      printf("\n");
    }  

  Matrix<Complex_f> m1t = m1.transpose();
  printf("\n");
  for(int r = 0; r < 3; ++r)
    {
      for (int c=0; c < 3; ++c)
	{
	  printf(" (%g,%g) ", m1t(r,c).real, m1t(r,c).imag);
	}
      printf("\n");
    }    
  printf("\n");
  Matrix<Complex_f> asdf(13,13);
  Matrix<Complex_f> id = asdf.identity();
  for(int r = 0; r < 13; ++r)
    {
      for (int c=0; c < 13; ++c)
	{
	  printf(" (%g,%g) ", id(r,c).real, id(r,c).imag);
	}
      printf("\n");
    }    

 
 
  printf("\n");
  for(int r = 0; r < 2; ++r)
    {
      for (int c=0; c < 2; ++c)
	{
	    printf(" (%g,%g) ", m2(r,c).real, m2(r,c).imag);
	}
      printf("\n");
    }
  
  // LinearAlgebra::gram_schmidt<Complex<float> >(Matrix<Complex<float> > const&)
  Matrix<Complex_f> mx = LinearAlgebra::gram_schmidt<Complex_f>(m1);
  printf("\n");
  for(int r = 0; r < 3; ++r)
    {
      for (int c=0; c < 3; ++c)
	{
	  printf(" (%g,%g) ", mx(r,c).real, mx(r,c).imag);
	}
      printf("\n");
    }

  // test orthogonality
  Vector<Complex_f> v1 = mx.column(0);
  Vector<Complex_f> v2 = mx.column(1);
  Vector<Complex_f> v3 = mx.column(2);

  Complex_f v12 = (v1.conj()).dot(v2);
  Complex_f v13 = (v2.conj()).dot(v3);
  Complex_f v23 = (v3.conj()).dot(v1);
  
  printf("\n\n");
  printf(" (%g,%g)\n", v12.real, v12.imag);
  printf("\n\n");

  printf(" (%g,%g)\n", v13.real, v13.imag);
  printf("\n\n");

  printf(" (%g,%g)\n", v23.real, v23.imag);
  printf("\n\n");
    
//for(int r=0; r<3; ++r) printf(" (%g,%g)\n", m2(r).real, m2(r).imag);

  //printf("\n");
  //printf(" (%g,%g)\n", m3.norm().real, m3.norm().imag);

return 0;
}
