#include <iostream>
#include <fstream>
#include <string>
#include "TCNT_net.h"
#include <vector>
using namespace std;
int main()
{
	double dE = 1 * 1.60219e-19;    //The height of the barrier 1eV
	int M = 460;                    //The total number of conduction channels
	double sigma_CNT = 1e4;         //The sigma of the CNT
	TCNT_net R_net(dE,M,sigma_CNT);

	vector<vector<vec3D>> CNT_net(2);
	vector<vector<int>> bc_Info(2);
	vec3D p1(0, 5, 0), p2(10, 5, 0), p3(0, 0, 0.1), p4(10, 10, 0.1);
	int num_divided=1;
	for (int i = 0; i < num_divided+1; i++)
	{
		vec3D pos1 = p1 + (p2 - p1)*i / num_divided;
		CNT_net[0].push_back(pos1);
		vec3D pos2 = p3 + (p4 - p3)*i / num_divided;
		CNT_net[1].push_back(pos2);
	}
	//CNT_net[0].push_back(p1); CNT_net[0].push_back(p2);
	//CNT_net[1].push_back(p3); CNT_net[1].push_back(p4);
	vector<int> index(12);
	for (int i = 0; i < 12; i++)
	{
		index[i] = 0;
	}
	bc_Info[0] = bc_Info[1] = index;
	bc_Info[0][0] = 1; bc_Info[1][9] = 1;
	R_net.input_CNT_net(CNT_net, bc_Info);
	R_net.Calculate_effective_resistance("x_min", "x_max");
	system("Pause");
	return 0;
	R_net.input_CNT_net();
	R_net.Calculate_effective_resistance("x_min", "x_max");
	system("Pause");
	return 0;
}