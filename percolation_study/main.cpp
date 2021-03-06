#include <iostream>
#include <fstream>
#include "TCNT_net.h"
#include "CNT_net_generation.h"
#include <vector>
#include <string>
using namespace std;
int main()
{
	ifstream input;
	input.open("RVE_Info.k");
	//Create the CNT network
	TCNT_net_generator CNT_net_generator;
	CNT_net_generator.input_RVE_info(input);
	double dE = 1 * 1.60219e-19;    //The height of the barrier 1eV
	int M = 2;                    //The total number of conduction channels
	double sigma_CNT = 1e4;         //The sigma of the CNT
	double lx, ly, lz;
	CNT_net_generator.G_RVE_size(lx, ly, lz);

	int num_test=100;
	ofstream output;
	output.open("R_list.k");
	for (int j = 0; j < 50; j++)
	{
		double fraction = (j + 1)*0.0005;
		CNT_net_generator.S_fraciton(fraction);
		double S_total = 0, S;
		for (int i = 0; i < num_test; i++)
		{
			vector<vector<vec3D>> CNT_net;
			vector<vector<int>> bc_Info;
			CNT_net_generator.Create_CNT_net(CNT_net, bc_Info,1);
			TCNT_net effective_resistance_net(dE, M, sigma_CNT);
			effective_resistance_net.input_CNT_net(CNT_net, bc_Info);
			double R = effective_resistance_net.Calculate_effective_resistance("x_min", "x_max",1);
			S = lx / (ly*lz*abs(R));
			//cout << i<<" "<<R << endl;
			//output << i << " " << R << endl;
			S_total = S_total + S;
		}
		double S_average = S_total / num_test;
		cout << fraction << " " << S_average << endl;
		output << fraction << " " << S_average << endl;
	}

	system("Pause");
	return 0;
}