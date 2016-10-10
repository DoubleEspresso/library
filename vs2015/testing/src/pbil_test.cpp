#include <stdio.h>
#include <iostream>

#include "../../../algorithm/mlearning/PBIL.h"

#define ABS(x) x < 0 ? -x : x

struct Param
{
	int bit_length;
	int answer;
};

float residual_function(int * sample, void * params)
{
	Param * p = ((Param*)(params));	
	int result = 0;
	for (int j = 0, s = p->bit_length-1; j < p->bit_length; ++j, --s)
	{
		result |= (sample[j] << s); // bit ordering is e.g. 4 = {1,0,0}
	}
	float r = result - p->answer;
	//printf("r=%d\n", r);
	return ABS(r);
}

int main()
{
	Param p;
	p.answer = 671265;
	p.bit_length = (int) log2(p.answer)+1;

	PBIL * plearn = new PBIL(300, p.bit_length);
	uint iterations = 0;
	float * probabilities = plearn->optimize((residual_func)&residual_function, (void*)&p, 0.15, 0.015, 0.3, 0.05, iterations, 1e-6);
	printf("..pbil converged in %d iterations\n", iterations);

	std::cin.get(); // wait here.
}