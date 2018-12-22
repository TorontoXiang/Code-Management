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
<<<<<<< HEAD
	//double r1, r2;
	//vec3D p1(78.314218338546411, 69.130193159469130, 80.760595164779161);
	//vec3D p2(36.222923029098986, 34.726814916524965, 119.92350134511798);
	//vec3D p3(6.2790739307357146, 122.72091788022473, 64.981533767145748);
	//vec3D p4(19.481940098969059, 117.65977187026066, 130.47254331420984);
	//vector<vec3D> curve1, curve2;
	//for (int i = 0; i < 36; i++)
	//{
	//	curve1.push_back(p1 + (p2 - p1)*i / 35.0);
	//	curve2.push_back(p3 + (p4 - p3)*i / 35.0);
	//}
	//Selectrical_node node1, node2;
	//cout << min_curve_dis(curve1, curve2, node1, node2) << endl;
	//cout << min_dis(p1, p2, p3, p4) << endl;
	//cout << min_dis(p1, p2, p3, p4, r1, r2)<<" ";
	//cout << r1 << " " << r2 << endl;
=======
	//create_sample_CNT();
	//return 0;
>>>>>>> 9845fca6fed3347720f97efaed443beeec089468
	ifstream input;
	input.open("CNT_grid_generation.k");
	int type, num_CNT, n_divided;
	double ratio;
	input >> type >> num_CNT >> n_divided >> ratio;
	if (type==0)
	{
		//Create stragiht CNT net
		Create_straight_CNT_net(num_CNT, 67, 67, 67, 67);
	}
	else if (type==1)
	{
		//Create wavy CNT net
		Create_wavy_CNT_net(num_CNT, 67, n_divided, ratio, 67, 67, 67);
	}
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

