#ifndef ALGORITHM_PBIL_H
#define ALGORITHM_PBIL_H

/*
Population based incremental learning -- 
*/
#include <vector>
#include <random>

#include "../../system/types.h"

typedef float (*residual_func)(int*, void*);

class PBIL
{
	size_t population, num_bits;
	float best_err;	
	std::mt19937 rng; // mersenne twister
public:
	PBIL(size_t pop_size, size_t nb_bits, float lrate) : population(pop_size), num_bits(nb_bits), best_err(1e10)
	{		
		rng.seed(std::random_device{}());
	}
	~PBIL() {}

	void optimize(residual_func cf, void * params, float learn_rate, float neg_learn_rate, float mutation_probabilty, float mutation_shift, uint iterations);
	int ** create_samples(float * probabilities);
	void update_probabilities(float * probabilities, int * min_gene, int * max_gene, float learn_rate, float neg_learn_rate);
	void mutate(float * probabilities, float mutation_probabilty, float mutation_shift);
};



int ** PBIL::create_samples(float * probabilities)
{
	std::uniform_real_distribution<double> dist(0, 1);
	int ** genes = new int*[population]; 
	for (int k = 0; k < population;++k)  genes[k] = new int[num_bits];
	for (int p = 0; p < population; ++p)
	{
		for (int j = 0; j < num_bits; ++j)
		{
			if (dist(rng) < probabilities[j]) genes[p][j] = 1;
			else genes[p][j] = 0;
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

void PBIL::optimize(residual_func rf, void * params, float learn_rate, float neg_learn_rate, float mutation_probabilty, float mutation_shift, uint iterations)
{
	if (population <= 0) population = 512;
	float * probabilities = new float[num_bits];	
	for (int i = 0; i < num_bits; ++i) probabilities[i] = 0.5f;
	int * best_sample = 0;

	for (int i = 0; i < iterations; i++)
	{
		// initialize population
		int ** samples = create_samples(probabilities);

		// compute costs
		float * errors = new float[population];
		for (int j = 0; j < population; ++j) errors[j] = rf(samples[j], params);

		// get max/min bounds
		int * min_sample = 0; int * max_sample = 0; 
		double min_err = 1e10; double max_err = -1e10;
		for (int j = 0; j < population; ++j)
		{
			float e = errors[j];
			if (min_err > e)
			{
				min_err = e;
				min_sample = samples[j];
			}
			if (max_err < e)
			{
				max_err = e;
				max_sample = samples[j];
			}
		}

		if (best_err > min_err)
		{
			best_err = min_err;
			best_sample = min_sample;
		}
		
		update_probabilities(probabilities, min_sample, max_sample, learn_rate, neg_learn_rate);
		mutate(probabilities, mutation_probabilty, mutation_shift);

		// free memory		
		if (samples) { for (int j = 0; j < population; ++j) { delete[] samples[j]; samples[j] = 0; } delete [] samples; samples = 0; }
		if (errors) { delete [] errors; errors = 0; }
	}
	if (probabilities) { delete[] probabilities; probabilities = 0; }
}

#endif