#ifndef ALGORITHM_PBIL_H
#define ALGORITHM_PBIL_H

/*
Population based incremental learning -- 
*/
#include <vector>
#include <random>

#include "../../system/types.h"

typedef float (*cost_func)(void*);

class PBIL
{
	size_t population, num_bits;
	float best_cost;	
	std::mt19937 rng; // mersenne twister
public:
	PBIL(size_t pop_size, size_t nb_bits, float lrate) : population(pop_size), num_bits(nb_bits), best_cost(1e10)
	{		
		rng.seed(std::random_device{}());
	}
	~PBIL() {}

	void optimize(cost_func cf, void * params, float learn_rate, float neg_learn_rate, float mutation_probabilty, float mutation_shift, uint iterations);
	int * load_genes(float * probabilities);
	void update_probabilities(float * probabilities, int * min_gene, int * max_gene, float learn_rate, float neg_learn_rate);
	void mutate(float * probabilities, float mutation_probabilty, float mutation_shift);
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

void PBIL::update_probabilities(float * probabilities, int * min_gene, int * max_gene, float learn_rate, float neg_learn_rate)
{
	for (int j = 0; j < num_bits; j++)
	{
		if (min_gene[j] == max_gene[j])
		{
			probabilities[j] = probabilities[j] * (1.0 - learn_rate) + min_gene[j] * learn_rate;
		}
		else
		{
			float learn_rate2 = learn_rate + neg_learn_rate;
			probabilities[j] = probabilities[j] * (1.0 - learn_rate2) + min_gene[j] * learn_rate2;
		}
	}
}

void PBIL::mutate(float * probabilities, float mutation_probabilty, float mutation_shift)
{
	std::uniform_real_distribution<double> dist(0, 1);
	for (int j = 0; j < num_bits; j++)
	{
		if (dist(rng) < mutation_probabilty)
		{
			probabilities[j] = probabilities[j] * (1.0 - mutation_shift) + (dist(rng) < 0.5 ? 1.0 : 0) * mutation_shift;
		}
	}
}

void PBIL::optimize(cost_func cf, void * params, float learn_rate, float neg_learn_rate, float mutation_probabilty, float mutation_shift, uint iterations)
{
	if (population <= 0) population = 512;
	float * probabilities = new float[num_bits];	
	for (int i = 0; i < num_bits; ++i) probabilities[i] = 0.5f;
	int * best_gene = new int[num_bits];

	for (int i = 0; i < iterations; i++)
	{
		// initialize population
		int * genes = load_genes(probabilities);

		// compute costs
		float * costs = new float[population];
		for (int j = 0; j < population; ++j) costs[j] = cf(params);

		// get max/min bounds
		int * min_gene = new int[num_bits]; int * max_gene = new int[num_bits]; 
		double min_cost = 1e10; double max_cost = -1e10;
		for (int j = 0; j < population; ++j)
		{
			float c = costs[j];
			if (min_cost > c)
			{
				min_cost = c;
				for (int k = population*num_bits, idx = 0; k < (population + 1)*num_bits; ++k, ++idx) min_gene[idx] = genes[k];
			}
			if (max_cost < c)
			{
				max_cost = c;
				for (int k = population*num_bits, idx = 0; k < (population + 1)*num_bits; ++k, ++idx) max_gene[idx] = genes[k];
			}
		}

		if (best_cost > min_cost)
		{
			best_cost = min_cost;
			for (int k = population*num_bits, idx = 0; k < (population + 1)*num_bits; ++k, ++idx) best_gene[idx] = genes[k];
		}
		
		update_probabilities(probabilities, min_gene, max_gene, learn_rate, neg_learn_rate);
		mutate(probabilities, mutation_probabilty, mutation_shift);

		// free memory
		if (min_gene) { delete min_gene; min_gene = 0; }
		if (max_gene) { delete max_gene; max_gene = 0; }
		if (costs) { delete costs; costs = 0; }
	}
}

#endif