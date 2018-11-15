#include <vector>
#include "data_structure.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "TCNT_grid.h"
#include "functions.h"
#include <Windows.h>

using namespace std;
int main()
{
	vec3D p_begin(0, 0, 0);
	int n_divided = 50;
	ofstream output_k,output_tec;
	output_k.open("CNT_grid.k");
	output_tec.open("CNT_grid.dat");
	//vec3D p1(10, 4, 1), p2(10, 3, 1), p3(10, 2, 1), p4(10, 1, 1), p5(10, -1, 1), p6(10, -2, 1), p7(10, -3, 1);
	//vector<vec3D> curve;
	//curve.push_back(p1); curve.push_back(p2); curve.push_back(p3); curve.push_back(p4); curve.push_back(p5);
	//curve.push_back(p6); curve.push_back(p7);
	//vector<vector<vec3D>> temp;
	//vec3D mini(0, 0, 0), maxi(100, 100, 100);
	//effective_curve(curve, temp, mini, maxi);
	Create_wavy_CNT_net(200, 67, 45,3.0, 67, 67, 67,860,0.17);
	//Generate_wavy_CNT(p_begin, 4.0, n_divided, 67);
	//Create_straight_CNT_net(300, 67, 67, 67, 67);
	//ofstream output_k,output_tec;
	//output_k.open("CNT_grid.k");
	//output_tec.open("CNT_grid.dat");
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
	//int i = 1, j = 2;

	system("Pause");
	return 0;
}

