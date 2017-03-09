#pragma once

#ifndef ALGORITHM_MLEARNING_NBAYES_H
#define ALGORITHM_MLEARNING_NBAYES_H

#include <random>

#include "..\..\math\matrix\matrix.h"
#include "..\..\math\matrix\vector.h"
#include "..\..\math\matrix\linalg.h"
#include "..\..\utils\fileio.h"

class NaiveBayes
{
	std::vector<float> csv_data;
	std::vector<std::vector<float>> means; // [sets][features]
	std::vector<std::vector<float>> sigmas; // [sets][features]
	std::vector<float> counts; // count of members in each set
	size_t ncols;
public:
	NaiveBayes(std::string fname, int cols, std::string format_string) : ncols(cols)
	{
		FileIO::read_csv<float>(fname, cols, format_string , csv_data);
	}
	~NaiveBayes() {}

	std::vector<float> data() { return csv_data; }
	void train(size_t sets)
	{
		if (csv_data.size() <= 0) return;
		
		// init an empty feature data array .. stores the means/sigmas of all data
		// for each class.
		means.clear();
		sigmas.clear();
		counts.clear();

		std::vector<float> tmp;
		for (int n = 0; n < ncols - 1; ++n) tmp.push_back(0);
		for (int j = 0; j < sets; ++j)
		{
			means.push_back(tmp);
			sigmas.push_back(tmp);
			counts.push_back(0);
		}

		// purpose: compute the means for each feature in each class (based on the following assumptions)..
		// note: last element is assumed to be the "classifier" .. data that tells us which class the feature set belongs to
		// further assumption: data is formated such that the last index per row will range from [0..k] for k-classes, i.e. it directly
		// indexes the class we are interested in.
		// even further assumption: that the data is not "jagged" e.g. rows * ncols = data.size(), so we can get nrows = data.size() / ncols;
		float nrow = (float) (ncols > 0 ? (float) csv_data.size() / (float) ncols : 1);
		for (int j = ncols - 1, row = 0; j < csv_data.size(); j += ncols, ++row)
		{
			++counts[csv_data[j]];
			// do not include the last col in this loop
			for (int k = 0; k < ncols - 1; ++k)
			{
				means[csv_data[j]][k] += csv_data[row * ncols + k];
			}
		}

		for (int c = 0; c < means.size(); ++c)
		{
			for (int j = 0; j < ncols - 1; ++j)
			{
				means[c][j] /= (float)counts[c];
			}
		}

		// purpose: compute the stdev for each feature given the means
		for (int j = ncols - 1, row = 0; j < csv_data.size(); j += ncols, ++row)
		{
			// do not include the last col in this loop
			int set = csv_data[j];
			for (int k = 0; k < ncols - 1; ++k)
			{
				int idx = row*ncols + k;
				sigmas[set][k] += (csv_data[idx] - means[set][k]) * (csv_data[idx] - means[set][k]);
			}
		}

		// take sqrt of each element. / n-1
		for (int c = 0; c < means.size(); ++c)
		{
			for (int k = 0; k < ncols - 1; ++k)
			{
				sigmas[c][k] = sqrt(sigmas[c][k] / ((float)counts[c] - 1));
			}
		}

		// dbg.
		/*for (int c = 0; c < means.size(); ++c)
		{
			for (int j = 0; j < ncols - 1; ++j)
			{
				printf("mean[%d][%d] = %.3f\tsigma[%d][%d] = %.3f\n", c, j, means[c][j], c,j, sigmas[c][j]);
			}
		}*/
	}
	int classify(std::vector<float> in)
	{
		if (in.size() != ncols - 1)
		{
			printf("..error: invalid input data for this training set, abort classification\n");
			return -1;
		}

		std::vector<float> probabilities;
		for (int c = 0; c < means.size(); ++c)
		{
			float log_sum = 0;
			float prob_class = counts[c] / (counts[0] + counts[1]);
			for (int j = 0; j < ncols - 1; ++j)
			{
				log_sum += log(prob_class) -
					(log(sigmas[c][j]) 
						+ (in[j] - means[c][j]) * (in[j] - means[c][j]) / (2.0 * sigmas[c][j] * sigmas[c][j]));
			}
			probabilities.push_back(log_sum);
		}
		
		float max = -99999999999;
		int class_idx = -1;
		for (int j = 0; j < probabilities.size(); ++j)
		{
			if (probabilities[j] > max)
			{
				max = probabilities[j];
				class_idx = j;
			}
		}
		return class_idx;
	}

	float validate(std::vector<float> validation_data)
	{
		int npass = 0; int row = 0;
		for (int j = ncols - 1; j < validation_data.size(); j += ncols, ++row)
		{
			int result = (int)validation_data[j];
			std::vector<float> test;
			for (int k = 0; k < ncols - 1; ++k) test.push_back(validation_data[row * ncols + k]);

			if (classify(test) == result)
			{
				++npass;
			}
		}
		return (float)npass / (float)(row);
	}
};


#endif