#include <vector>
#include "data_structure.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "TCNT_grid.h"
#include "functions.h"
#include <Windows.h>
#include "CNT_net_generation.h"
using namespace std;
void create_sample_CNT();
int main()
{
	ifstream input;
	input.open("RVE_Info.k");
	//Create the CNT network
	TCNT_net_generator CNT_net_generator;
	CNT_net_generator.input_RVE_info(input);
	vector<vector<vec3D>> CNT_net;
	vector<vector<int>> bc_Info;
	//vector<vec3D> test_CNT;
	//vec3D p1(10, 10, -10), p2(10, 10, 77), p3(30, 30, 77), p4(30, 30, -10), p5(50, 50, -10), p6(60, 60, 20);
	//test_CNT.push_back(p1), test_CNT.push_back(p2);
	//test_CNT.push_back(p3), test_CNT.push_back(p4);
	//test_CNT.push_back(p5), test_CNT.push_back(p6);
	//CNT_net_generator.Clipping_curved_CNT(test_CNT, CNT_net, bc_Info);
	//return 0;
	ofstream output_k, output_tec, output_bc;
	output_k.open("CNT_grid.k"); output_tec.open("CNT_grid.dat"); output_bc.open("bc_Info.k");
	//CNT_net.resize(2);
	//bc_Info.resize(2);
	//vec3D p1(0, 5, 10), p2(10, 5, 10), p3(0, 0, 11), p4(15, 15, 11);
	//CNT_net[0].push_back(p1); CNT_net[0].push_back(p2);
	//CNT_net[1].push_back(p3); CNT_net[1].push_back(p4);
	//bc_Info.resize(2);
	//bc_Info[0].resize(12); bc_Info[1].resize(12);
	//for (int i = 0; i < 12; i++)
	//{
	//	bc_Info[0][i] = bc_Info[1][i] = 0;
	//}
	//bc_Info[0][0] = 1; bc_Info[1][9] = 1;
	//CNT_net_generator.Generate_CNT_net_grid(CNT_net, bc_Info, output_k, output_tec, output_bc);
	//return 0;
	//if (CNT_net_generator.G_type()==0)
	//{
	//	//Create stragiht CNT net
	//	CNT_net_generator.Create_straight_CNT_net(CNT_net, bc_Info);
	//	CNT_net_generator.Generate_CNT_net_grid(CNT_net, bc_Info, output_k, output_tec, output_bc);
	//	//Create_straight_CNT_net(num_CNT, 67, 67, 67, 67);
	//}
	//else if (CNT_net_generator.G_type()==1)
	//{
	//	//Create wavy CNT net
	//	CNT_net_generator.Create_wavy_CNT_net(CNT_net, bc_Info);
	//	CNT_net_generator.Generate_CNT_net_grid(CNT_net, bc_Info, output_k, output_tec, output_bc);
	//}
	CNT_net_generator.Create_CNT_net(CNT_net, bc_Info);
	CNT_net_generator.Generate_CNT_net_grid(CNT_net, bc_Info, output_k, output_tec, output_bc);
	system("Pause");
	return 0;
}
void create_sample_CNT()
{
	ofstream output_k, output_tec;
	output_k.open("CNT_grid.k");
	output_tec.open("CNT_grid.dat");
	ifstream input;
	input.open("CNT_info.k");
	vector<vec3D> p1, p2;
	int num_CNT;
	input >> num_CNT;
	p1.resize(num_CNT); p2.resize(num_CNT);
	double dl = 2;
	for (int i = 0; i < num_CNT; i++)
	{
		input >> p1[i].x >> p1[i].y >> p1[i].z >> p2[i].x >> p2[i].y >> p2[i].z;
		double l = (p2[i] - p1[i]).get_length();
		cout << l << endl;
		int n_divided = maxval(int(l / dl),2);
		Generate_straight_CNT(p1[i], p2[i], n_divided, 0.335, 0.335*0.6, 0, 1, output_k, output_tec, 860, 0.17);
	}
	output_k << "*END" << endl;
	system("Pause");
	//vec3D p1[11], p2[11];
	//p1[0].value(31.38, 61.37, 0), p2[0].value(12.33907, 40.00896, 20.10098);
	//p1[1].value(0.301, 59.98483, 39.38792); p2[1].value(30.5128, 6.702012, 66.53935);
	//p1[2].value(3.67, 59.24, 27.72); p2[2].value(35.7607, 7.308873, 55.3306);
	//p1[3].value(35.25034, 48.88831, 47.38919); p2[3].value(45.01, 0, 24.04);
	//p1[4].value(0, 7.26, 15.37); p2[4].value(55.73623, 26.83868, 8.440469);
	//p1[5].value(60.51, 67, 50.82); p2[5].value(55.38022, 7.732504, 58.23969);
	//p1[6].value(43.32, 67, 13.87); p2[6].value(48.96965, 34.332, 59.70719);
	//p1[7].value(17.66929, 10.66876, 17.1103); p2[7].value(49.61, 67, 3.41);
	//p1[8].value(67, 22.44, 66.61); p2[8].value(16.23993, 15.41999, 48.79998);
	//p1[9].value(48.56213, 21.3542, 56.4203); p2[9].value(67, 57.76, 59);
	//p1[10].value(36.4, 0, 61.7); p2[10].value(48.82045, 13.94051, 48.57952);
	//double dl = 0.67;
	//for (int i = 0; i < 11; i++)
	//{
	//	double l = (p2[i] - p1[i]).get_length();
	//	int n_divided = int(l / dl);
	//	Generate_straight_CNT(p1[i], p2[i], n_divided, 0.335, 0.335*0.6, 1, 3, output_k,output_tec, 860, 0.17);
	//}
	//output_k << "*END" << endl;
	return;
}

