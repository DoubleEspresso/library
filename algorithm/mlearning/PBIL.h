#ifndef ALGORITHM_PBIL_H
#define ALGORITHM_PBIL_H

/*
Population based incremental learning -- 
*/
#include <vector>
#include <random>

#include "../../utils/bitops.h"

typedef float (*cost_func)(void*);

class PBIL
{
	size_t population, num_bits;
	float best_cost;	
	std::mt19937 rng; // mersenne Twister: Good quality random number generator
public:
	PBIL(size_t pop_size, size_t nb_bits) : population(pop_size), num_bits(nb_bits), best_cost(-1) 
	{		
		rng.seed(std::random_device{}());
	}
	~PBIL() {}

	void optimize(cost_func cf, void*params, uint iterations);
	int * load_genes(float * probabilities);
};

int * PBIL::load_genes(float * probabilities)
{
	std::uniform_real_distribution<double> dist(0, 1);
	int * genes = new int[population * num_bits];
	for (int i = 0, idx = 0; i < population; ++i)
	{
		for (int j = 0; j < num_bits; ++j, ++idx)
		{
			if (dist(rng) < probabilities[j]) genes[idx] = 1;
			else genes[idx] = 0;
		}
	}
	return genes;
}

void PBIL::optimize(cost_func cf, void * params, uint iterations)
{
	if (population <= 0) population = 512;
	float * probabilities = new float[num_bits];
	
	for (int i = 0; i < num_bits; ++i) probabilities[i] = 0.5f;

	// ...or until convergence ?
	for (int i = 0; i < iterations; i++)
	{
		// initialize population
		int * genes = load_genes(probabilities);

		// compute costs
		float * costs = new float[population];
		for (int j = 0; j < population; ++j) costs[j] = cf(params);

		// get max/min bounds
		int * min_gene = 0; int * max_gene = 0;
		double min_cost = 1e10; double max_cost = -1e10;
		for (int j = 0; j < population; ++j)
		{
			float c = costs[j];
			if (min_cost > c)
			{
				min_cost = c;
				//min_gene = genes[j];
			}
			if (max_cost < c)
			{
				max_cost = c;
				//max_gene = genes[j];
			}
		}

		if (best_cost > min_cost)
		{
			best_cost = min_cost;
			//best_gene = min_gene;
		}
	}


}

#endif