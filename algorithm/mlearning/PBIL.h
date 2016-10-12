#ifndef ALGORITHM_PBIL_H
#define ALGORITHM_PBIL_H

/*
Population based incremental learning -- 
*/
#include <vector>
#include <random>

#include "../../system/types.h"
#include "../../system/threads.h"
#include "../../utils/timer.h"

typedef float (*residual_func)(int *, void*);

struct pbil_thread_data
{
	int start, stop;
	int ** samples;
	void * work_params;
	residual_func * work_task;
	float * error;
};

class PBIL
{
	size_t population, num_bits;
	float best_err;	float * probabilities;
	int * best_sample;
	std::mt19937 rng; // mersenne twister
public:
	PBIL(size_t pop_size, size_t nb_bits) : population(pop_size), num_bits(nb_bits), best_err(1e10), probabilities(0), best_sample(0)
	{		
		rng.seed(std::random_device{}());
	}
	~PBIL() 
	{
		if (probabilities) { delete[] probabilities; probabilities = 0; }
		if (best_sample) { delete[] best_sample; best_sample = 0; }
	}
	int * best() { return best_sample; }
	void optimize(residual_func cf, void * params, float learn_rate, float neg_learn_rate, float mutation_probabilty, float mutation_shift, uint& iterations, double tolerance);
	void optimize_parallel(residual_func rf, void * params, float learn_rate, float neg_learn_rate, float mutation_probabilty, float mutation_shift, uint& iterations, uint nb_threads, double tolerance);
	static void residual_parallel(void * params);
	void educate(int ** samples, float * probabilities);
	void update_probabilities(float * probabilities, int * min_gene, int * max_gene, float learn_rate, float neg_learn_rate);
	void mutate(float * probabilities, float mutation_probabilty, float mutation_shift);
};


void PBIL::educate(int ** samples, float * probabilities)
{
	std::uniform_real_distribution<double> dist(0, 1);
	for (int p = 0; p < population; ++p)
	{
		for (int j = 0; j < num_bits; ++j)
		{
			if (dist(rng) < probabilities[j]) samples[p][j] = 1;
			else samples[p][j] = 0;
		}
	}

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

void PBIL::optimize(residual_func rf, void * params, float learn_rate, float neg_learn_rate, float mutation_probabilty, float mutation_shift, uint& iterations, double tolerance)
{
	if (population <= 0) population = 512;
	probabilities = new float[num_bits];	
	for (int i = 0; i < num_bits; ++i) probabilities[i] = 0.5f;
	int * best = 0;
	Timer clock;
	std::vector<double> times;
	// initialize population
	int ** samples = new int*[population];
	for (int k = 0; k < population; ++k)
	{
		samples[k] = new int[num_bits];
		for (int i = 0; i < num_bits; ++i) samples[k][i] = 0;
	}

	
	while(best_err >= tolerance)
	{
		educate(samples, probabilities);

		// compute costs
		clock.start();
		float * errors = new float[population];
		for (int j = 0; j < population; ++j) errors[j] = rf(samples[j], params);
		clock.stop();

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
			best = min_sample;
		}
		
		update_probabilities(probabilities, min_sample, max_sample, learn_rate, neg_learn_rate);
		mutate(probabilities, mutation_probabilty, mutation_shift);
		if (++iterations % 200 == 0)
		{
			times.push_back(clock.elapsed_ms());
			printf("iter(%d) ", iterations);
			for (int j = 0; j < num_bits; ++j) printf("%d", best[j]);
			printf("\n");
		}

		// free memory		
		if (errors) { delete [] errors; errors = 0; }
	}	
	if (best != 0)
	{
		best_sample = new int[num_bits];
		memcpy(best_sample, best, num_bits * sizeof(int));
		printf("best: ");
		for (int j = 0; j < num_bits; ++j) printf("%.3f ", probabilities[j]);
		printf("\n");
		for (int j = 0; j < num_bits; ++j) printf("%d", best[j]);
		printf("\n");
		printf("error = %.3f\n", best_err);
		float avgtime = 0;
		for (int j = 0; j < times.size(); ++j) avgtime += times[j] / times.size();
		printf("avgerage update time = %.3fms\n", avgtime);
	}
	if (samples) { for (int j = 0; j < population; ++j) { delete[] samples[j]; samples[j] = 0; } delete[] samples; samples = 0; }

}

void PBIL::optimize_parallel(residual_func rf, void * params, float learn_rate, float neg_learn_rate, float mutation_probabilty, float mutation_shift, uint& iterations, uint nb_threads, double tolerance)
{
	if (population <= 0) population = 512;
	probabilities = new float[num_bits];
	for (int i = 0; i < num_bits; ++i) probabilities[i] = 0.5f;
	int * best = 0;
	Timer clock;

	// initialize population
	int ** samples = new int*[population];
	for (int k = 0; k < population; ++k)
	{
		samples[k] = new int[num_bits];
		for (int i = 0; i < num_bits; ++i) samples[k][i] = 0;
	}

	std::vector<pbil_thread_data*> thread_data;
	std::vector<double> times;
	int remainder = population % nb_threads;
	int batch_size = (population / nb_threads) + remainder;
	THREAD_HANDLE * threads = new THREAD_HANDLE[nb_threads];
	int start = 0;

	// thread initialization
	for (int j = 0; j < nb_threads; ++j)
	{
		pbil_thread_data * td = new pbil_thread_data();
		td->work_task = &rf; 
		td->work_params = params; // thread independent parameters
		td->start = start;
		td->stop = td->start + population / nb_threads;
		start = td->stop;
		// allocate error arrays
		uint error_size = td->stop - td->start;
		td->samples = new int*[error_size];
		for (int k = 0; k < error_size; ++k)
		{
			td->samples[k] = new int[num_bits];
			for (int i = 0; i < num_bits; ++i) td->samples[k][i] = 0;
		}
		td->error = new float[error_size]; memset(td->error, 0, error_size * sizeof(float));
		//threads.push_back(new Thread(j, (thread_fnc)residual_parallel, (void*)td));
		
		thread_data.push_back(td);
	}

	while (best_err >= tolerance)
	{
		educate(samples, probabilities);

		// update samples & compute costs (in parallel)
		for (int tid = 0, idx = 0, start = 0; tid < nb_threads; ++tid, start = idx)
		{
			for (int j = start, sidx = 0; j < start + population / nb_threads; ++j, ++idx, ++sidx)
			{
				thread_data[tid]->samples[sidx] = samples[j];
			}
		}
		// catch remaining members..

		// launch threads
		clock.start();
		for (int j = 0; j < nb_threads; ++j)
		{
			threads[j] = start_thread((thread_fnc)residual_parallel, (void*)thread_data[j]);;
		}
		wait_threads_finish(threads, nb_threads);
		clock.stop();

		// copy errors from worker threads to main thread
		float * errors = new float[population];
		for (int tid = 0, idx = 0, start = 0; tid < nb_threads; ++tid, start = idx)
		{
			for (int j = start, eidx = 0; j < start + population / nb_threads; ++j, ++idx, ++eidx)
			{
				errors[j]= thread_data[tid]->error[eidx]; 
			}
		}
		// catch remaining errors

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
			best = min_sample;
		}

		update_probabilities(probabilities, min_sample, max_sample, learn_rate, neg_learn_rate);
		mutate(probabilities, mutation_probabilty, mutation_shift);
		if (++iterations % 200 == 0)
		{
			printf("iter(%d) ", iterations);
			for (int j = 0; j < num_bits; ++j) printf("%d", best[j]);
			times.push_back(clock.elapsed_ms());
			printf("\n");
		}

		// free memory		
		if (errors) { delete[] errors; errors = 0; }
	}
	if (best != 0)
	{
		best_sample = new int[num_bits];
		memcpy(best_sample, best, num_bits * sizeof(int));
		printf("best: ");
		for (int j = 0; j < num_bits; ++j) printf("%.3f ", probabilities[j]);
		printf("\n");
		for (int j = 0; j < num_bits; ++j) printf("%d", best[j]);
		printf("\n");
		printf("error = %.3f\n", best_err);
		float avgtime = 0;
		for (int j = 0; j < times.size(); ++j) avgtime += times[j] / times.size();
		printf("avgerage update time = %.3fms\n", avgtime);

	}
	if (samples) { for (int j = 0; j < population; ++j) { delete[] samples[j]; samples[j] = 0; } delete[] samples; samples = 0; }
}

void PBIL::residual_parallel(void * params)
{
	pbil_thread_data * p = (pbil_thread_data*)params;
	for (int j = p->start, idx = 0; j < p->stop; ++j, ++idx)
	{
		// note : samples is passed as a sub-array (indexed from 0 to (stop-start)) of total population
		// same as error array.
		 float ret = (*p->work_task)(p->samples[idx], p->work_params);
		 p->error[idx] = ret;
	}
}

#endif