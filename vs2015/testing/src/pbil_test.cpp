#include <stdio.h>
#include <iostream>

#include "../../../algorithm/mlearning/PBIL.h"

#define ABS(x) x < 0 ? -x : x

struct Param
{
	int bit_length;
	int answer;
};

struct Param2
{
	size_t rows, cols;
	int bits_per_int;
	int total_length;
	int answer;
};

int integer(int * sample, int sidx, int eidx)
{
	int result = 0;
	for (int j = sidx, s = eidx - 1 - sidx; j < eidx; ++j, --s)
	{
		result |= (sample[j] << s); // bit ordering is e.g. 4 = {1,0,0}
	}
	return result;
}

// compute error for toy problem : find binary representation of a given "answer" integer.
float residual_function1(int * sample, void * params)
{
	Param * p = ((Param*)(params));	
	float r = integer(sample, 0, p->bit_length) - p->answer;
	return ABS(r);
}

// compute error for toy problem : find set of numbers summing to 
// "answer" along the rows/columns of an MxN matrix.
float residual_function2(int * sample, void * params)
{
	Param2 * p = ((Param2*)(params));
	int * rowsum = new int[p->rows]; int * colsum = new int[p->cols];
	float rowerr = 0, colerr = 0;
	memset(rowsum, 0, p->rows * sizeof(int)); memset(colsum, 0, p->cols * sizeof(int));

	for (int r = 0, idx = 0; r < p->rows; ++r)
	{
		for (int c = 0; c < p->cols; ++c, idx += p->bits_per_int)
		{
			rowsum[r] += integer(sample, idx, idx + p->bits_per_int);
			colsum[c] += integer(sample, idx, idx + p->bits_per_int); 
		}
		
	}
	for (int r = 0; r < p->rows; ++r) rowerr += (rowsum[r] - p->answer) * (rowsum[r] - p->answer);
	for (int c = 0; c < p->cols; ++c) colerr += (colsum[c] - p->answer) * (colsum[c] - p->answer);
	//printf("rowerr=%.3f, colerr=%.3f\n", rowerr, colerr);
	// free memory
	if (rowsum) { delete[] rowsum; rowsum = 0; }
	if (colsum) { delete[] colsum; colsum = 0; }
	return rowerr + colerr;
}

int main()
{
	/*toy problem 1 : find binary representation of a given input integer*/
	/*
	Param p;
	p.answer = 671265;
	p.bit_length = (int) log2(p.answer)+1;

	PBIL * plearn = new PBIL(300, p.bit_length);
	uint iterations = 0;
	float * probabilities = plearn->optimize((residual_func)&residual_function1, (void*)&p, 0.15, 0.015, 0.3, 0.05, iterations, 1e-6);
	printf("..pbil converged in %d iterations\n", iterations);
	*/

	/*toy problem 2 : find digits representing rows/cols of MxN matrix whos rows/cols sum to answer*/
	Param2 p2;
	p2.answer = 94;
	p2.rows = 4; 
	p2.cols = 4;
	p2.bits_per_int = ((int)log2(p2.answer) + 1);
	p2.total_length = p2.rows * p2.cols * p2.bits_per_int; // total bits for each member of the population
	PBIL * plearn = new PBIL(8000, p2.total_length);
	uint iterations = 0;
	plearn->optimize((residual_func)&residual_function2, (void*)&p2, 0.15, 0.015, 0.3, 0.05, iterations, 1e-6);
	//plearn->optimize_parallel((residual_func)&residual_function2, (void*)&p2, 0.15, 0.015, 0.3, 0.05, iterations, 4, 1e-6);
	printf("..pbil converged in %d iterations\n", iterations);
	for (int r = 0, idx = 0; r < p2.rows; ++r)
	{
		for (int c = 0; c < p2.cols; ++c, idx += p2.bits_per_int)
		{
			printf("%d ", integer(plearn->best(), idx, idx + p2.bits_per_int)); 
		}
		printf("\n");
	}

	std::cin.get(); // wait here.
}