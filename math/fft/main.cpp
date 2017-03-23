#include <iostream>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <cmath> 
#include <stdlib.h> // rand
#include <time.h>

#include "fft.h"


// To compile this code using gcc:
// g++ -o fastexp -O3 -fomit-frame-pointer -msse2 -mfpmath=sse -ffast-math

void parse_args(int argc, char * argv[]);

int main(int argc, char* argv[])
{

  parse_args(argc, argv);
  fft f;

  // fft data
  size_t N = 2048;
  double * data = new double[N];
  double freq = 3;
  
  for (unsigned int j=0; j<N; ++j)
    {
      data[j] = sin(PI2 * freq * ((double)((double)(j+1) / (double)N)));
    }
  
  // dbg
  //double * data = new double[N];
  //data[0] = 1; data[1] = 2; data[2] = 3; data[3] = 4;
  
  complex * coeffs = new complex[N];
  double start = clock();
  f.forward_1d(data, coeffs, N);
  printf("..fft time=%3.3fms\n", (clock() - start) / CLOCKS_PER_SEC * 1000.0 );
  /*
  for (unsigned int j=0; j<N; ++j)
    {
      if (coeffs[j].re != 0 || coeffs[j].im != 0) printf("%d (%3.3f, %3.3f)\n", j, coeffs[j].re, coeffs[j].im);
    }
  */
  /*dbg - naive method comparison*/
  start = clock();
  f.forward_naive(data, coeffs, N);  
  printf("..naive time=%3.3fms\n", (clock() - start) / CLOCKS_PER_SEC * 1000.0 );
  /*for (unsigned int j=0; j<N; ++j)
    {
    if (coeffs[j].re != 0 || coeffs[j].im != 0) printf("%d (%3.3f, %3.3f)\n", j, coeffs[j].re, coeffs[j].im);
    }*/

  //complex tmp = QEXP::root(5, 1024);
  //printf("(%3.6f, %3.6f)\n", tmp.re, tmp.im);
  return 0;
}


void parse_args(int argc, char * argv[])
{
  if (argc < 2)
    {
      printf("..usage: ./fft [options] <input data>\n");
      printf("  where options include:\n");
      printf("  -f <input filename> : input data file name to be transformed\n");
      printf("  -o <output filename> : input data file name to be transformed\n");
    }
  for(int j=1; j<argc; ++j) // skip first option
    {
      if (strcmp(argv[j], "-f") == 0)
	{ 
	  printf("..parse input filename (todo)\n");
	}
      else
	{
	  printf("..unknown option: %s\n", argv[j]);
	}
    }
}
