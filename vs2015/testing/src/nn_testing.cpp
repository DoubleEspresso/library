#include <stdio.h>
#include <iostream>
#include <random>

#include "../../../math/matrix/matrix.h"
#include "../../../algorithm/mlearning/NN.h"


void parse_args(int argc, char*argv[]);
float _sig(void *p)
{
	float * pf = (float*)(p);
	return 1.0f / (1.0f + exp(-1.0f * (*pf)));
}
float _dsig(void *p)
{
	return _sig(p) * (1.0 - _sig(p));
}
float _cost(void *p)
{
	return 0; // not needed
}
float _dcost(void *p)
{
	float * pf = (float*)p;
	return *pf;
}

int main(int argc, char * argv[])
{
	parse_args(argc, argv);

	std::mt19937 _rng; // mersenne twister
	std::uniform_real_distribution<double> dist(0, 1);


	/*load data here*/
	int training_size = 50;
	Matrix<float> * training_data = new Matrix<float>(50, 2);

	for (int j = 0; j < training_size; ++j)
	{
		float x = j; float y = sin(x * x) + dist(_rng);
		training_data->set(j, 0, x);
		training_data->set(j, 1, y);
	}
	/*network params*/
	float lrate = 0.1; // learning rate
	size_t sample_size = 2; // batch size for sgd 
	size_t tepochs = 4; // nb of training epochs

	int * nn_dims = new int[3];
	nn_dims[0] = training_data->nb_rows();
	nn_dims[1] = 4;
	nn_dims[2] = 1;

	//Network(Matrix<float> * indata, size_t hidden_layers, int * layer_dims, size_t num_layers, size_t n_runs)
	Network * net = new Network(training_data, 1, nn_dims, 3, 2);
	net->init_funcs((net_func)&_sig, (net_func)&_dsig, (net_func)&_cost, (net_func)&_dcost);
	if (!net->sgd(training_data, lrate, sample_size, tepochs))
	{
		printf("..NN stochastic gradient error, failed to train network\n");
		return -1;
	}
	printf("..NN finished training OK\n");
	std::cin.get();
	return 0;
}

/* -- methods implementations -- */
void parse_args(int argc, char * argv[])
{

}
