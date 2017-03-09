#include <stdio.h>
#include <iostream>
#include <random>

#include "../../../math/matrix/matrix.h"
#include "../../../algorithm/mlearning/NN.h"
#include "../../../utils/fileio.h"
#include "../../../algorithm/mlearning/NaiveBayesClassifier.h"

void parse_args(int argc, char*argv[]);

double _sig(void *p)
{
	float * pf = (float*)(p);
	return  1.0 / (float)(1.0 + expf(-1.0 * (*pf)));
}
double train_func(void *x)
{
	float * xf = (float*)(x);
	// (*xf < 0.0 ? 0.0 : 1.0); 
	// sin(2.0*(*xf))*sin(2.0*(*xf));
	// ((*xf > -0.5 && *xf  < 0.5) ? 0.0 : 1.0);
	// sin(4.0*(*xf))*sin(4.0*(*xf))
	return sin(2.0*(*xf))*sin(2.0*(*xf));
}
double _dsig(void *p)
{
	return 1.0 * _sig(p) * (1.0 - _sig(p));
}
double _cost(void *p) // cross entropy
{
	Matrix<float> * ay = (Matrix<float>*)(p);
	float a = ay->data_at(0, 0);
	float y = ay->data_at(0, 1);
	return -1.0 * ( y * log(a) + (1.0 - y) * log(1 - a));
}
double _dcost(void *p)
{
	float * pf = (float*)p;
	return 1.0 * (*pf);
}

int main(int argc, char * argv[])
{
	parse_args(argc, argv);

	NaiveBayes nb("A:\\code\\test-data\\prima-indians-data.csv", 9, "%f%c");
	nb.train(2);
	std::vector<float> d = nb.data();
	float acc = nb.validate(d);
	printf("..%.3f percent accuracy\n", acc * 100.0);
	
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
		// 0.25 * (2.0*x*x + cos(10.0*x) * cos(10.0*x) + sin(2.0*x)*sin(2.0*x))
		double x = dist(rng); float y = sin(2.0*x)*sin(2.0*x);
		training_data->set(j, 0, x);
		training_data->set(j, 1, y);
	}

	/*network params*/
	float lrate = 0.14; // learning rate
	size_t sample_size = 10000; // batch size for sgd 
	size_t tepochs = 16000; // nb of training epochs

	int * nn_dims = new int[3];
	nn_dims[0] = training_data->nb_rows();
	nn_dims[1] = 4;
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
