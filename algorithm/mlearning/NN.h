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


								  // abstraction of a network transfer parameters from one layer to the next
								  // to start: a link-layer is defined by a transfer matrix W_ij, and bias b_i
class Link
{
	// transfer matrices definitions
	// basic transfer operation on input vector 
	// Vi' = W_{ij}.Vi_{j} + B_i (similar for each layer)

	// relevant layer data
	Matrix<float>* _W; // transfer matrix [rows(L+1:dim) x cols(L:dim)], e.g. transfer from L-1->L layer
	Matrix<float>* _b; // bias vector [rows(L+1:dim)]

public:
	Link() : _W(0), _b(0) { }
	Link(size_t r, size_t c)
	{
		_W = new Matrix<float>(r, c);
		_b = new Matrix<float>(r, 1);
	}
	~Link()
	{
		if (_W) { delete _W; _W = 0; }
		if (_b) { delete _b; _b = 0; }
	}
	bool init(std::mt19937& rng)
	{
		std::uniform_real_distribution<double> dist(-0.5, 0.5);
		for (int r = 0; r < _W->nb_rows(); ++r)
		{
			for (int c = 0; c < _W->nb_cols(); ++c)
			{
				_W->set(r, c, (float)dist(rng));
				if (c == 0) _b->set(r, c, (float)dist(rng));
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
	// note : in stochastic gradient descent, we average over the batch of gradients
	// to update (fixing L), these inputs are given by dW = 1/m*sum_m (dC/dW) etc.
	void update_params(const Matrix<float>& dW, const Matrix<float>& dB)
	{
		// occasionally fails .. (different sizes in - operator).
		(*_W) = (*_W) - dW; // gradient descent for weights
		(*_b) = (*_b) - dB; // gradient descent for biases
	}
};

// represents an array of activations and raw-values for all "neurons"
// in a "hidden-layer" .. the hidden portion of the neural network
class Nodes
{
	// note : these are all *hidden* nodes by construction
	Matrix<float>* _z; // storage -- neuron value vector (result of the transfer operation from layer L-1)
	Matrix<float>* _a; // storage -- activation vector for this layer
	Matrix<float>* _d; // dCost/dBias derivative (used to compute the gradients on the fly)
	size_t _size;
public:
	Nodes() : _z(0), _a(0), _d(0), _size(0) {};
	Nodes(size_t num)
		: _z(0), _a(0), _d(0), _size(num)
	{
		_z = new Matrix<float>(num, 1);
		_a = new Matrix<float>(num, 1);
		_d = new Matrix<float>(num, 1);
		assert(_z->nb_rows() == _size && _a->nb_rows() == _size && _d->nb_rows() == _size);
	}
	~Nodes()
	{
		if (_z) { delete _z; _z = 0; }
		if (_a) { delete _a; _a = 0; }
		if (_d) { delete _d; _d = 0; }
	}
	Matrix<float> * a() { return _a; }
	Matrix<float> * d() { return _d; }
	Matrix<float> * z() { return _z; }
	size_t dim() { return _size; }
	void clear()
	{
		_z->clear();
		_a->clear();
		_d->clear();
	}
	void compute_activations(net_func _f)
	{
		assert(_f != 0 && _z != 0);
		for (int r = 0; r < _size; ++r)
		{
			float v = _z->data_at(r, 0); // assumed to be columnwise-data
			_a->set(r, 0, (float)(_f((void*)&v)));
		}
	}

	void compute_zvals(Matrix<float> * a, Link * l)
	{
		assert(l != 0);
		assert(l->W()->nb_cols() == a->nb_rows());
		assert(l->b()->nb_rows() == l->W()->nb_rows());
		assert(_z->nb_rows() == l->W()->nb_rows() && _z->nb_cols() == 1);
		(*_z) = (*l->W()) * (*a) + (*l->b());
	}

	// note: these are hidden layer deltas only, the input/output layer data
	// is stored in the main network class
	void compute_deltas(Matrix<float> * ds, Link * l, Matrix<float> * dsig)
	{
		// W_ij^^ * di --> W_ji di --> d'j 
		assert(l != 0);
		assert(l->W()->nb_rows() == ds->nb_rows());
		assert(dsig->nb_rows() == l->W()->nb_cols());

		//assert(ds->nb_rows() == dsig->nb_rows());
		Matrix<float> wd(l->W()->nb_cols(), 1);
		wd = l->W()->transpose() * (*ds);

		// set the delta-vector by multiply previous result by dsigma
		for (int r = 0; r < dsig->nb_rows(); ++r)
		{
			_d->set(r, 0, wd.data_at(r, 0) * dsig->data_at(r, 0)); //schur product
		}
	}
};


class Network // to be made a base class in the future ?
{
	Matrix<float> * _data; // input training data (nx2 matrix of (x,y) - input, output tuples)
	Link ** _L; // params for transfer from input to first hidden layer layout: [input->hidden, hidden->hidden, hidden->end]
	Nodes ** _nodes; // stores all activations and z-values for each "neuron" in a given "hidden" layer.
	// means if there are N-hidden layers, there are N+1-link layers
	size_t _links; // number of hidden layers
	int * _dims; // dimensions for each layer (including in/output layers)
	net_func _sig; // activation function for the network
	net_func _dsig; // derivative of sigmoid
	net_func _cost; // error function for the network
	net_func _dcost; // derivative error function for the network
	bool _valid; // status variable of network
	std::mt19937 _rng; // mersenne twister
public:
	Network() : _data(0), _L(0), _links(0), _dims(0),
		_sig(0), _dsig(0), _cost(0), _dcost(0), _valid(false)
	{
		_rng.seed(std::random_device{}());
	}
	Network(Matrix<float> * indata, size_t hidden_layers, int * layer_dims, size_t num_layers)
		: _data(indata), _L(0), _links(hidden_layers + 1), _dims(layer_dims),
		_sig(0), _dsig(0), _cost(0), _dcost(0), _valid(false)
	{
		if (num_layers != hidden_layers + 2)
		{
			printf("..ERROR : mismatch between layer dimension array size and requested number of network layers, abort.\n");
			return;
		}
		if (_links < 2)
		{
			printf("..ERROR : need at least 1-hidden layer of neurons\n");
			return;
		}
		_rng.seed(std::random_device{}());

		_L = new Link*[_links]; 		
		for (int j = 0; j < _links; ++j)
		{
			if (_dims[j] <= 0 || _dims[j + 1] <= 0)
			{
				printf("..ERROR: invalid layer dim(0) for network-layer(%d)\n", (_dims[j] <= 0 ? j : j + 1));
				return;
			}

			if (j == 0) _L[j] = new Link(_dims[j + 1], 1); // input->hidden link
			else if (j == _links - 1) _L[j] = new Link(1, _dims[j]); // hidden->output link
			else _L[j] = new Link(_dims[j + 1], _dims[j]); // hidden->hidden link
			
			if (!_L[j] || !_L[j]->init(_rng))
			{
				printf("..ERROR: failed to initialize parameters of network, aborting.\n");
				return;
			}
		}
		// note: there are n=hidden_layers of nodes* and n+1 Link data pointers!!
		_nodes = new Nodes*[hidden_layers];
		for (int j = 0; j < hidden_layers; ++j)
		{
			// note: this allocates a set of vectors to store activations/z-vals etc. (column vector)
			_nodes[j] = new Nodes(_dims[j + 1]); // (j+1) since _dims starts @ input layer
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
	// @params
	//	_data : training data for layer L-1 .. nx2 matrix of input/output tuples
	//  ia : input activation (for first layer only)
	//	_f : user-defined activation function
	// @output
	//	stores result in activation "_a", _data is updated (by reference) 
	//  to values for this layer L (to be used for layer L+1) 
	void forward_step(int idx, Matrix<float> * xin) //moves input layer to first hidden layer 
	{
		// transfer operation from layer @ layer=idx, to layer=idx+1
		// note : W_ij(L-1) * data_i + b_i --> produces vector(size ix1)
		// compute z-vals === Link(idx) * activations(idx-1) + b(idx)
		// note that nodes start idx @ hidden layer, and Links start idx @ input data (off by 1 each).
		Nodes * curr_nodes = _nodes[idx]; // @ idx

		// compute the activations for the input data on the first step
		if (idx == 0)
		{
			assert(xin->nb_rows() == _L[idx]->W()->nb_cols() && xin->nb_cols() == 1);
			for (int r = 0; r < xin->nb_rows(); ++r)
			{
				float v = xin->data_at(r, 0); // assumed to be columnwise-data
				xin->set(r, 0, (float)(_sig((void*)&v))); // note: this is updated by reference
			}
			curr_nodes->compute_zvals(xin, _L[idx]);
			curr_nodes->compute_activations(*_sig);
			return;
		}
		Nodes * prev_nodes = _nodes[idx - 1];
		curr_nodes->compute_zvals(prev_nodes->a(), _L[idx]);
		curr_nodes->compute_activations(*_sig);
	}

	void backward_step(int idx, Nodes * output,  Matrix<float> * dcost)
	{
		if (idx == _links)
		{
			// load the output deltas
			// dcost represents the derivative of the cost vector
			assert(output != 0 && dcost != 0);
			assert(dcost->nb_rows() == output->z()->nb_rows());
			for (int r = 0; r < dcost->nb_rows(); ++r)
			{
				float v = output->z()->data_at(r, 0); // assumed to be columnwise-data
				// note : this udpates the output-deltas by references (we need to store those)
				float tmp = (float)(_dsig((void*)&v) * dcost->data_at(r, 0));
				output->d()->set(r, 0, tmp);
			}
			return;
		}
		Nodes * curr_nodes = _nodes[idx - 1]; // pointer to the hidden nodes 
		Nodes * prev_nodes = 0;
		prev_nodes = (idx == _links - 1 ? output : _nodes[idx]);
		assert(curr_nodes != 0);
		assert(prev_nodes != 0);
		curr_nodes->compute_deltas(prev_nodes->d(), _L[idx], compute_sigp(curr_nodes->z()));
	}

	// @input
	//		xin : an (possible) nx1 matrix of training data to be fed to the network .. 
	//			forward pass will compute each node value, and activation, given a hidden layer
	//			of weights and biases, from the input layer to the output layer
	Nodes * forward_pass(Matrix<float> * xin)
	{
		{
			// accounts for input layer and hidden layers 
			for (int j = 0; j < _links - 1; ++j) forward_step(j, xin);
		}
		Nodes * output = new Nodes(_dims[_links]); // the output layer dimensions
		assert((int)output->dim() == xin->nb_rows()); // these dimensions should match! (x,y) tuples for data is assumed!
		{
			// output layer here ...
			output->compute_zvals(_nodes[_links - 2]->a(), _L[_links - 1]);
			output->compute_activations(*_sig);
		}
		assert(output != 0);
		return output;
	}
	// @input
	//		cost : an (possible) nx2 matrix of computed costs .. 
	//			backward pass computes the dW, and dB deltas used to construct the gradients
	//			for gradient descent updates of the weights and biases.
	void backward_pass(Nodes * output, Matrix<float> * dcosts)
	{
		assert(output != 0 && dcosts != 0);
		{
			// backward step -- iterate over the layer-idx (from links + 1 to 0) 
			for (int j = _links; j >= 1; --j) backward_step(j, output, dcosts);
		}
	}

	// given an rx1 column vector, return an rx1 column
	// vector of sig' values
	Matrix<float> * compute_sigp(Matrix<float> * Z)
	{
		assert(_dsig != 0);
		assert(Z->nb_rows() > 0 && Z->nb_cols() == 1);
		Matrix<float> * res = new Matrix<float>(*Z);
		for (int r = 0; r < Z->nb_rows(); ++r)
		{
			float v = Z->data_at(r, 0);
			res->set(r, 0, (float)_dsig((void*)&v));
		}
		return res;
	}

	// stochastic (sampled) gradient descent
	// @input
	//		_data: matrix of nx2 tuples (x,y) input training data
	//		lrate: learning rate for the network
	//		size: sample size of the training data to be used (batch size)
	//		samples: number of training epochs
	bool sgd(Matrix<float> * _data, float lrate, size_t size, size_t samples)
	{
		for (int sample = 0; sample < samples; ++sample) // loop over training epochs (number of samplings)
		{
			// 1. randomize input data for each training epoch
			Matrix<float> * shuffled = new Matrix<float>(*_data);
			_data->shuffle_rows(*shuffled);

			// allocate storage for network parameters needed in gradient descent step
			int out_size = _dims[_links]; // for 1d data, this should be a 1x1 scalar
			Matrix<float> * dcosts = new Matrix<float>(out_size, 1);
			std::vector<Matrix<float>> dW;
			std::vector<Matrix<float>> dB;

			// 2. loop over training sample - of n-rows (defines a sampled batch from the training data)
			Matrix<float> * batch = shuffled->get_rows(0, size); // batch = [size]x[2] matrix

			// note: we assume 1d data, so here a batch = [nx2] array .. loop over each x,y pair 
			// todo: generalize the concept of "x" it could also be an array/image for example, as could the y-value.
			for (int i = 0; i < batch->nb_rows(); ++i)
			{
				// for testing, the input data is a collection of (xy) tuples
				Matrix<float> singleton(1, 1);
				singleton.set(0, 0, batch->data_at(i, 0));
				
				// 3. compute network parameters relevant for gradient descent (back-propagation)
				Nodes * output = forward_pass(&singleton); // note: singleton is updated by reference

				// 4. load cost-vector for backward pass
				Matrix<float> * y = output->a();
				assert(y->nb_rows() == dcosts->nb_rows());
				float v = y->data_at(0, 0) - batch->data_at(i, 1);
				// float tmp = (float)_dcost((void*)&v); returning 0 always !?
				dcosts->set(0, v);
				//dcosts->print("..dcost value first pass..");
				// 5. backward pass with dcost vector - will store all delta params at each layer of the network 
				// to be used for gradient descent update of weights and biases for each layer of the network.
				backward_pass(output, dcosts); 

				// 6. store the computed delta's for each batch (the deltas are overwritten for each training epoch)
				// for each link-layer, adjust the stored weights with current values of deltas and activations
				for (int j = 0; j < _links; ++j)
				{
					Nodes * n = (j == _links - 1 ? output : _nodes[j]);
					if (j == 0)
					{
						if (i == 0) dW.push_back((*n->d()) * singleton.transpose());
						// note : singleton = activation (by ref)
						// note : nodes[j=0] is the first layer of hidden nodes!
						else dW[j] = dW[j] + ((*n->d()) * singleton.transpose());
					}
					else
					{
						if (i == 0) dW.push_back((*n->d()) * _nodes[j - 1]->a()->transpose());
						else dW[j] = dW[j] + ((*n->d()) * _nodes[j - 1]->a()->transpose());
					}
					if (i == 0) dB.push_back((*n->d()));
					else dB[j] = dB[j] + (*n->d());
				}

				for (int j = 0; j < _links - 1; ++j) _nodes[j]->clear();
			}
		
			// 7. gradient descent update of the weights/biases using stored dW & dB arrays

			float scale = lrate / batch->nb_rows();
			for (int j = 0; j < _links; ++j)
			{
				_L[j]->update_params(scale * dW[j] , scale * dB[j]);
			}

			// 8. trace progress
			if (sample%100==0) printf("..finished training epoch %d/%d\n", sample, samples);
		
		}

		// 9. dbg - readout final network bias/weight params 
		for (int j = 0; j < _links; ++j)
		{
			_L[j]->trace();
		}

		return true;
	}

	void verify(Matrix<float> * testdata, net_func _f)
	{
		printf("..validation\n");
		float rms = 0.0f;
		for (int r = 0; r < testdata->nb_rows(); ++r)
		{
			Matrix<float> S(1, 1);
			S.clear();
			S.set(0, 0, testdata->data_at(r, 0));
			// forward pass on input data
			for (int j = 0; j < _links; ++j)
			{
				Matrix<float> * wij = _L[j]->W(); Matrix<float> * bi = _L[j]->b();
				
				S = (*wij) * S + (*bi);

				Matrix<float> T(S);
				//if (j != _links - 1)
				{
					for (int i = 0; i < S.nb_rows(); ++i)
					{
						float sval = S.data_at(i, 0);
						T.set(i, 0, _sig((void*)&sval));
					}
				}
				S = T;
			}
			//net_results.set(r, 0, singleton.data_at(0, 0));
			float tmpv = testdata->data_at(r, 0);
			float y = (float)_f((void*)&tmpv);
			printf("%d \t%1.6f \t%1.6f \t%1.6f\n", r, S.data_at(0, 0), y, S.data_at(0, 0) - y);
			rms += (S.data_at(0, 0) - y) * (S.data_at(0, 0) - y);
		}
		rms = sqrt(rms / testdata->nb_rows());
		printf("..rms-deviation=%.3f\n", rms);
	}
};

#endif