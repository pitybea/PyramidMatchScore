#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
using namespace std;
#include <stdio.h>
#include <cstring>
#include <string>

struct dircFisher
{
	double thresh;
	int direction;
	int dimensn;
	double alpha;
	dircFisher()
	{

	}
};

class adaboost
{
public:
	int dimension;
	int num_wk_cl;
	int num_smps;
	int num_P_smps;
	int num_N_smps;

	vector<vector<double>> train_smps;
	vector<vector<int>> orders;
	vector<int> train_label;
	vector<double> smp_weight;
	//vector<double> alphas;
	vector<int> needTo;

	vector<dircFisher> str_cl;
	vector<dircFisher> candida_str_cl;
	vector<double> wk_cl_weight;
	vector<double> wk_cl_errors;

	void intitWts();
	void ordering();
	double train_wk_cl(int);
	void training();
	double classfy(vector<double>);
	int singleClassfy(dircFisher dF,vector<double> smp);
	void outputTofile(FILE* fp);
	void loadFromfile(FILE* fp);

};