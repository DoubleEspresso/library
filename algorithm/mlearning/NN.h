#pragma once
#ifndef ALGORITHM_NN_H
#define ALGORITHM_NN_H

#include <random>

#include "..\..\math\matrix\matrix.h"
#include "..\..\math\matrix\vector.h"
#include "..\..\math\matrix\linalg.h"

// sandbox code for simple artificial neural network
// on the CPU.

typedef double(*net_func)(void*); // sigmoid, cost func & derivatives. 

// abstraction of a network layer
// to start: a layer is defined by a single transfer matrix W_ij, and bias b_i
class Layer
{
	// transfer matrices definitions
	// basic transfer operation on input vector 
	// Vi' = W_{ij}.Vi_{j} + B_i (similar for each layer)
	bool _input; // specialized operations needed
	bool _output; // specialized operations needed

	// relevant layer data
	Matrix<float>* _W; // transfer matrix [rows(L+1:dim) x cols(L:dim)], e.g. transfer from L-1->L layer
	Matrix<float>* _b; // bias vector [rows(L+1:dim)]
	Matrix<float>* _z; // storage -- neuron value vector (result of the transfer operation from layer L-1)
	Matrix<float>* _a; // storage -- activation vector for this layer
	Matrix<float>* _d; // dCost/dBias derivative (used to compute the gradients on the fly)

public:
	Layer() : _input(false), _output(false), _W(0), _b(0), _z(0), _a(0), _d(0) { }
	Layer(bool input_layer, bool output_layer, size_t r, size_t c) : _input(input_layer), _output(output_layer)
	{
		if (input_layer)
		{
			_W = new Matrix<float>(r, 1);
			_a = new Matrix<float>(1, 1);
		}
		else if (output_layer)
		{
			_W = new Matrix<float>(r, 1);
			_a = new Matrix<float>(1, 1);
		}
		else
		{
			_W = new Matrix<float>(r, c); // c = #rows of this layer, r = #rows of next layer
			_a = new Matrix<float>(r, 1);
		}
		_b = new Matrix<float>(r, 1);
		_z = new Matrix<float>(r, 1);
		_d = new Matrix<float>(r, 1);
	}
	~Layer()
	{
		if (_W) { delete _W; _W = 0; }
		if (_b) { delete _b; _b = 0; }
		if (_z) { delete _z; _z = 0; }
		if (_a) { delete _a; _a = 0; }
		if (_d) { delete _d; _d = 0; }
	}
	bool init(std::mt19937& rng)
	{
		std::uniform_real_distribution<double> dist(0, 1);
		for (int r = 0; r < _W->nb_rows(); ++r)
		{
			for (int c = 0; c < _W->nb_cols(); ++c)
			{
				_W->set(r, c, (float)dist(rng));
				if (r == 0) _b->set(0, c, (float)dist(rng));
			}
		}
		return true;
	}
	void trace()
	{
		_W->print("w-vals");
		_b->print("b-vals");
	}
	Matrix<float> * W() { return _W; }
	Matrix<float> * b() { return _b; }
	Matrix<float> * a() { return _a; }
	Matrix<float> * d() { return _d; }
	Matrix<float> * z() { return _z; }
	// @params
	//	_data : training data for layer L-1 .. nx2 matrix of input/output tuples
	//	_f : user-defined activation function
	// @output
	//	stores result in activation "_a", _data is updated (by reference) 
	//  to values for this layer L (to be used for layer L+1) 
	void forward_step(float xin, net_func _f, Layer * l = 0) // e.g. moves input layer to first hidden layer 
	{
		if (_input)
		{
			// the input activation is assumed to be a scalar for now [1x1] matrix
			_a->set(0, 0, (float)_f((void*)&xin));
			return;
		}
		// transfer operation from (L-1)-layer to L layer
		// note : W_ij(L-1) * data_i + b_i --> produces vector(size ix1)
		Matrix<float> w = (*l->W()); 
		Matrix<float> a = (*l->a());
		_z->set((*l->W()) * (*l->a()) + (*_b));
		for (int r = 0; r < _z->nb_rows(); ++r)
		{
			float v = _z->data_at(r, 0);
			_a->set(r, 0, (float)_f((void*)&v)); // apply activation function, and store result
		}
	}
	// @params
	//	 dy : dCost vector = dC/dy_i (where C = cost_func(y_i-t_i), t=target vector (answer-data)
	//		  and y_i is the result of a feed-forward operation.  dy is assumed to be a column vector of size (output-dim,1)
	// sigp : function pointer to sigma'(y_i) @ layer L (derivative of cost function)
	//   l  : pointer to L+1 layer (used only if not output layer)
	// @output
	//   _d : stores the elements of an error vector used to compute the gradients for 
	//        gradient descent (backpropagation step of training)
	//		  note, formula is : _d_i(L) = d/dy (Cost(y_i)) * sig'(z_i)
	void compute_delta(net_func sigp, Matrix<float>* dy, Layer * l = 0)
	{
		if (sigp == 0) return;
		Matrix<float> ds(_d->nb_rows(), 1);
		for (int j = 0; j < ds.nb_rows(); ++j)
		{
			float v = _z->data_at(j, 0);
			ds.set(j, 0, (float)(sigp)((void*)&v));
		}
		if (_output)
		{
			for (int j = 0; j < _d->nb_rows(); ++j)
			{
				_d->set(j, 0, dy->data_at(j, 0) * ds.data_at(j, 0));
			}
		}
		else // layer L delta is computed from layer L+1 delta, and layer L W_ij
		{
			// backpropagate _d@(L+1) by applying the transpose W_ij
			if (l == 0) return;
			Matrix<float> v = _W->transpose() * (*l->d()); // W_{L+1, L}_T * d_{L+1,1}
			for (int j = 0; j < _d->nb_rows(); ++j)
			{
				_d->set(j, 0, v.data_at(j, 0) * ds.data_at(j, 0));
			}
		}
	}
	// note : in stochastic gradient descent, we average over the batch of gradients
	// to update (fixing L), these inputs are given by dW = 1/m*sum_m (dC/dW) etc.
	void update_params(float rate, const Matrix<float>& dW, const Matrix<float>& dB)
	{
		// occasionally fails .. (different sizes in - operator).
		(*_W) = (*_W) - rate * dW; // gradient descent for weights
		(*_b) = (*_b) - rate * dB; // gradient descent for biases
	}
};

class Network // to be made a base class in the future ?
{
	Matrix<float> * _data; // input training data (nx2 matrix of (x,y) tuples)
	Layer ** _L; // params for transfer from input to first hidden layer layout: [input->hidden, hidden->hidden, hidden->end]
	size_t _links; // number of hidden layers
	int * _dims; // dimensions for each layer (including in/output layers)
	size_t _runs; // number of training iterations
	net_func _sig; // activation function for the network
	net_func _dsig; // derivative of sigmoid
	net_func _cost; // error function for the network
	net_func _dcost; // derivative error function for the network
	bool _valid; // status variable of network
	std::mt19937 _rng; // mersenne twister
public:
	Network() : _data(0), _L(0), _links(0), _dims(0), _runs(0),
		_sig(0), _dsig(0), _cost(0), _dcost(0), _valid(false)
	{
		_rng.seed(std::random_device{}());
	}
	Network(Matrix<float> * indata, size_t hidden_layers, int * layer_dims, size_t n_runs)
		: _data(indata), _L(0), _links(hidden_layers + 1), _dims(layer_dims), _runs(n_runs),
		_sig(0), _dsig(0), _cost(0), _dcost(0), _valid(false)
	{
		if (_links < 2)
		{
			printf("..ERROR : need at least 1-hidden layer of neurons\n");
			return;
		}
		_rng.seed(std::random_device{}());

		_L = new Layer*[_links]; // includes input/output layer		

		for (int j = 0, k = _dims[j]; j < _links; k = _dims[++j])
		{
			if (k > 0)
			{
				_L[j] = new Layer((j == 0), (j == _links - 1), (j>1 ? _dims[j - 1] : _dims[0]), _dims[j]);

				if (!_L[j] || !_L[j]->init(_rng))
				{
					printf("..ERROR: failed to initialize parameters of network, aborting.\n");
					return;
				}
			}
			else
			{
				printf("..ERROR: invalid layer dim(0) for network-layer(%d)\n", j);
				return;
			}
		}
	}
	~Network()
	{
		if (_data) { delete _data; _data = 0; }
		if (_L) { for (int j = 0; j < _links; ++j) delete _L[j]; delete[] _L; _L = 0; }
		if (_dims) { delete _dims; _dims = 0; }
		//if (_sig) { delete _sig; _sig = 0; }
		//if (_dsig) { delete _dsig; _dsig = 0; }
		//if (_cost) { delete _cost; _cost = 0; }
		//if (_dcost) { delete _dcost; _dcost = 0; }
	}
	void init_funcs(net_func s, net_func sp, net_func ef, net_func ep)
	{
		_sig = (net_func)s; _dsig = (net_func)sp; _cost = (net_func)ef; _dcost = (net_func)ep;
		_valid = (_sig != 0 && _dsig != 0 && _cost != 0 && _dcost != 0);
	}
	bool isValid() { return _valid; }

	float compute_loss(Matrix<float>& y)
	{
		// check dims[y] == _dims[_layers-1] .. 
		float r = 0;
		for (int j = 0; j < _dims[_links - 1]; ++j)
		{
			float delta = y.data_at(j, 0) - y.data_at(j, 1); // y *assumed be* a nx2 matrix
			r += (float)(*_cost)((void*)&delta);
		}
		return r;
	}
	// @input
	//		xin : an (possible) nx2 matrix of training data to be fed to the network .. 
	//			forward pass will compute each node value, and activation, given a hidden layer
	//			of weights and biases, from the input layer to the output layer
	void forward_pass(float xin)
	{
		for (int j = 0; j < _links; ++j) _L[j]->forward_step(xin, *_sig, (j == 0 ? 0 : _L[j-1]));
	}
	// @input
	//		cost : an (possible) nx2 matrix of computed costs .. 
	//			backward pass computes the dW, and dB deltas used to construct the gradients
	//			for gradient descent updates of the weights and biases.
	void backward_pass(Matrix<float> * cost_vec)
	{
		for (int j = _links - 1; j >= 0; --j)
		{
			_L[j]->compute_delta(_dsig,
				(j == _links - 1 ? cost_vec : 0),
				(j < _links - 1 ? _L[j + 1] : 0));
		}
	}

	// stochastic (sampled) gradient descent
	// @input
	//		_data: matrix of nx2 tuples (x,y) input training data
	//		lrate: learning rate for the network
	//		size: sample size of the training data to be used (batch size)
	//		samples: number of training epochs
	bool sgd(Matrix<float> * _data, float lrate, size_t size, size_t samples)
	{
		// TODO: sanity checking: enforce size <= _data->nb_rows()/2, else set size = _data->nb_rows(); 
		// if !valid return false; etc.
		for (int sample = 0; sample < samples; ++sample) // loop over training epochs (number of samplings)
		{
			// 1. randomize input data for each training epoch
			Matrix<float> * shuffled = new Matrix<float>(*_data);
			_data->shuffle_rows(*shuffled);

			// allocate storage for network parameters needed in gradient descent step
			int out_size = _dims[_links - 1]; // for 1d data, this should be a 1x1 scalar
			Matrix<float> * costs = new Matrix<float>(out_size, 1);
			std::vector<Matrix<float>> dW;
			std::vector<Matrix<float>> dB;

			// 2. loop over training sample - of n-rows (defines a sampled batch from the training data)
			Matrix<float> * batch = shuffled->get_rows(0, size); // batch = [size]x[2] matrix

			// note: we assume 1d data, so here a batch = [nx2] array .. loop over each x,y pair 
			// todo: generalize the concept of "x" it could also be an array/image for example, as could the y-value.
			for (int i = 0; i < batch->nb_rows(); ++i)
			{
				// 3. compute network parameters relevant for gradient descent (back-propagation)
				forward_pass(batch->data_at(i, 0)); // updates all node values, and activation vectors

				// 4. load cost-vector for backward pass
				Matrix<float> * y = _L[_links - 1]->z(); // exception (?)
				for (int r = 0; r < costs->nb_rows(); ++r)
				{
					float v = y->data_at(r, 0) - batch->data_at(r, 1); // todo : _dcost should accept (y_i, t_i) as input + void-params, and return a floating point result 
					costs->set(r, 0, (*_dcost)((void*)&v));
				}

				// 5. backward pass with cost vector - will store all delta params at each layer of the network 
				// to be used for gradient descent update of weights and biases for each layer of the network.
				backward_pass(costs); 

				// 6. store the computed delta's for each batch (the deltas are overwritten for each training epoch)
				if (i == 0) // first data point
				{
					for (int j = 0; j < _links; ++j)
					{
						dW.push_back((*_L[j]->d()) * _L[j]->a()->transpose());
						dB.push_back((*_L[j]->d()));
					}
				}
				else
				{
					for (int j = 0; j < _links; ++j)
					{
						dW[j] = dW[j] + (*_L[j]->d()) * _L[j]->a()->transpose();
						dB[j] = dB[j] + (*_L[j]->d());
					}
				}
			}
		
			// 7. gradient descent update of the weights/biases using stored dW & dB arrays
			for (int j = 1; j < _links; ++j)
			{
				_L[j]->update_params(lrate, 1.0f / dW.size() * dW[j], 1.0f / dB.size() * dB[j]);
			}

			// 8. trace progress
			printf("..finished training epoch %d/%d\n", sample, samples);
		
		}
		return true;
	}
};

#endif