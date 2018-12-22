#include <fstream>
#include <iostream>
#include <vector>
#include "TResistance_net.h"
using namespace std;
int main()
{
	ifstream input;
	input.open("Resistance_net.k");
	TResisitance_net Rnet;
	Rnet.input_resistance_net(input);
	Rnet.Generate_Sparse_Matrix();
	double r = Rnet.Calculate_effective_resistance();
	cout << abs(r) << endl;
	system("Pause");
	return 0;
}