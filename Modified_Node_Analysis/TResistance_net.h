#pragma once
#include <vector>
#include<fstream>
#include<cmath>
#include "TCNT_net.h"
using namespace std;
int maxval(int a, int b);
struct Sresistance;
struct Snode
{
	vector<int> _connected_node;
	vector<double> _connected_resistance;

	double Calculate_diagonal_element();
	double Calculate_non_diagonal_element(int j);
};
class TResisitance_net
{
public:
	void input_resistance_net(ifstream &input);
	void input_resistance_net(vector<Sresistance>);
	//Input the resistance net
	void Generate_Sparse_Matrix();
	//Generate the sprase matrix of the resistance net
	double Calculate_effective_resistance();
	//Calculate the effective resistance of the net
private:
	vector<Snode> _node_list;  //The resistance net
	int _nv1, _nv2;            //The positive and negative polar of the volgate
	double _v;                 //The value of the volgate
	int _num_freedom;          //Number of unknown variables
	int _num_node;             //Number of node

	int* _iK;
	int* _jK;
	double* _K;                   //Sparse matrix of the resistance system

	//void Generate_Sparse_Matrix();
	////Generate the sprase matrix of the resistance net
};

