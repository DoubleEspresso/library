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
	int r = result - p->answer;
	float e = (float)((float)ABS(r));// / (float)p->answer);
	//printf("..e = %.6f\n", e);
	return e;
}

int main()
{
	Param p;
	p.answer = 671265;// 3235791;
	p.bit_length = (int) log2(p.answer)+1;
	//printf("..bit_length = %d", p.bit_length);
	//int * sample = new int[p.bit_length] { 1,1,0,1,1,1,1,0,1,0,1,1,0 };
	//float err = residual_function(sample, (void*)&p);
	//printf("...sanity check=%.3f\n", err);

	PBIL * plearn = new PBIL(800, p.bit_length);
	plearn->optimize((residual_func)&residual_function, (void*)&p, 0.7, 0.075, 0.02, 0.05, 500);
	//residual_func rf, void * params, float learn_rate, float neg_learn_rate, float mutation_probabilty, float mutation_shift, uint iterations
	//printf("best=");
	//for (int j = 0; j < p.bit_length; ++j) printf("%d", best[j]);
	//printf("\n");
	std::cin.get(); // wait here.
}