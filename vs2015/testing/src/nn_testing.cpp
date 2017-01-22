#include <stdio.h>
#include <iostream>
#include <random>

#include "../../../math/matrix/matrix.h"
#include "../../../algorithm/mlearning/NN.h"


void parse_args(int argc, char*argv[]);

float _sig(void *p)
{
	float * pf = (float*)(p);
	return 1.0 / (1.0 + exp(-1.0 * (*pf)));
}
float _sig_inv(float x)
{
	return -1.0 * log(1.0 / x - 1.0);
}
float train_func(void *x)
{
	float * xf = (float*)(x);
	// (*xf < 0.0 ? 0.0 : 1.0); 
	// sin(2.0*(*xf))*sin(2.0*(*xf));
	// ((*xf > -0.5 && *xf  < 0.5) ? 0.0 : 1.0);
	return sin(4.0*(*xf))*sin(4.0*(*xf));
}
float _dsig(void *p)
{
	return  1.0 * _sig(p) * (1.0 - _sig(p));
}
float _cost(void *p)
{
	return 0; // not needed
}
float _dcost(void *p)
{
	float * pf = (float*)p;
	return (float)1.0 * (*pf);
}

int main(int argc, char * argv[])
{
	parse_args(argc, argv);

	//std::mt19937 _rng; // mersenne twister
	//std::uniform_real_distribution<double> dist(0, 1);


	/*load data here*/
	int training_size = 40000;
	Matrix<float> * training_data = new Matrix<float>(training_size, 2);


	std::uniform_real_distribution<double> dist(-1, 1);
	std::mt19937 rng;
	for (int j = 0; j < training_size; ++j)
	{	
		//(x < 0.0 ? 0.0 : 1.0);
		//sin(2.0*x)*sin(2.0*x);
		// ((x > -0.5 && x < 0.5) ? 0.0 : 1.0);
		double x = dist(rng); float y = sin(4.0*x)*sin(4.0*x);
		training_data->set(j, 0, x);
		training_data->set(j, 1, y);
	}

	/*network params*/
	float lrate = 15.14; // learning rate
	size_t sample_size = 600; // batch size for sgd 
	size_t tepochs = 16000; // nb of training epochs

	int * nn_dims = new int[3];
	nn_dims[0] = training_data->nb_rows();
	nn_dims[1] = 5;
	nn_dims[2] = 1;

	//Network(Matrix<float> * indata, size_t hidden_layers, int * layer_dims, size_t num_layers)
	Network * net = new Network(training_data, 1, nn_dims, 3);
	net->init_funcs((net_func)&_sig, (net_func)&_dsig, (net_func)&_cost, (net_func)&_dcost);
	if (!net->sgd(training_data, lrate, sample_size, tepochs))
	{
		printf("..NN stochastic gradient error, failed to train network\n");
		return -1;
	}
	printf("..NN finished training OK\n");

	// test trained network
	int tpts = 20;
	Matrix<float> xdata(20, 1);
	for (int j = 0; j < 20; ++j)
	{
		xdata.set(j, 0, dist(rng));
	}
	net->verify(&xdata, (net_func)&train_func);

	std::cin.get();
	return 0;
}

/* -- methods implementations -- */
void parse_args(int argc, char * argv[])
{

}
