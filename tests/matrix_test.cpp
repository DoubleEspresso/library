#include <stdio.h>
#include <cmath>
#include <vector>

#include "../math/base/complex.h"
#include "../math/matrix/matrix.h"
#include "../math/matrix/linalg.h"
#include "../utils/timer.h"

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




  Matrix<float> m1(4,4); // should compute in parallel automatically
  Matrix<float> m2(4,4);
  /*
  m1.set(0,Complex_f(1,0));
  m1.set(1,Complex_f(-1,0));
  m1.set(2,Complex_f(0,0));
  m1.set(3,Complex_f(0,0));
  m1.set(4,Complex_f(0,1));
  m1.set(5,Complex_f(-1,0));
  m1.set(6,Complex_f(0,1));
  m1.set(7,Complex_f(1,0));
  m1.set(8,Complex_f(1,1));
  */

  m1.set(0,4); m1.set(1,1); m1.set(2,-2); m1.set(3,2);
  m1.set(4,1); m1.set(5,2); m1.set(6,0); m1.set(7,1);
  m1.set(8,-2); m1.set(9,0); m1.set(10,3); m1.set(11,-2);
  m1.set(12,2); m1.set(13,1); m1.set(14,-2); m1.set(15,-1);

  m2.set(0,1); m2.set(1,-1); m2.set(2,2); m2.set(3,2);
  m2.set(4,-1); m2.set(5,2); m2.set(6,1); m2.set(7,-1);
  m2.set(8,2); m2.set(9,1); m2.set(10,3); m2.set(11,2);
  m2.set(12,2); m2.set(13,-1); m2.set(14,2); m2.set(15,1);
  
  for(int r = 0; r < 4; ++r)
    {
      for (int c=0; c < 4; ++c)
	{
	  printf(" (%g) ", m2(r,c));//.real, m1(r,c).imag);
	}
      printf("\n");
    }  

  /*
  Matrix<float> res = LinearAlgebra::hessenberg_form<float>(m1, 0);
  
  printf("\n");
  for(int r = 0; r < 4; ++r)
    {
      for (int c=0; c < 4; ++c)
	{
	  printf(" (%g) ", res(r,c));//, res(r,c).imag);
	}
      printf("\n");
    }
  */
  Matrix<float> res2 = LinearAlgebra::tridiagonal_householder<float>(m2);
  printf("\n");
  for(int r = 0; r < 4; ++r)
    {
      for (int c=0; c < 4; ++c)
	{
	  printf(" (%g) ", res2(r,c));//, res(r,c).imag);
	}
      printf("\n");
    }
  
  /*
  //Matrix<Complex_f> mx = LinearAlgebra::gram_schmidt<Complex_f>(m1);
  Matrix<Complex_f> Q(3,3);
  Matrix<Complex_f> R(3,3);
  LinearAlgebra::QR<Complex_f>(m1, Q, R);

  printf("\n");
  for(int r = 0; r < 3; ++r)
    {
      for (int c=0; c < 3; ++c)
	{
	  printf(" (%g,%g) ", Q(r,c).real, Q(r,c).imag);
	}
      printf("\n");
    }
      printf("\n");      printf("\n");
  for(int r = 0; r < 3; ++r)
    {
      for (int c=0; c < 3; ++c)
	{
	  printf(" (%g,%g) ", R(r,c).real, R(r,c).imag);
	}
      printf("\n");
    }
  */

return 0;
}
