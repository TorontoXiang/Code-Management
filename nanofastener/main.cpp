#include "Tsubstrate.h"
#include <iostream>
#include <fstream>
#include <omp.h>
using namespace std;
double straight_straight_coupled(double density, double r_NW, double r_P, double l_NW, double theta1, double theta2);
double wavy_wavy_coupled(double density, double r_NW, double r_P, double l_NW, double curvature1, double curvature2);
double straight_wavy_coupled(double density, double r_NW, double r_P, double l_NW, double curvature, double theta);
int main()
{
	omp_set_num_threads(16);
	srand((int)time(0));
	ofstream output_tec,output_p;
	output_tec.open("subtrate.dat");
	output_p.open("result_p.k");
	Tsubstrate coupled_substrate;
	bool is_success = false;
	double density,p;
	//for (int i = 1; i < 41; i++)
	//{
	//	cout << i << endl;
	//	density = 0.01*i;
	//	//p = straight_straight_coupled(density, 15, 0, 150, 0, 0);
	//	//p = straight_straight_coupled(density, 15, 0, 150, 0.1, 0.1);
	//	//p = wavy_wavy_coupled(density, 15, 0, 150, 0.2, 0.2);
	//	//p = straight_wavy_coupled(density, 15, 0, 150, 0.2, 0.1);
	//	output_p << density << " " << p<<endl;
	//	if (p==0)
	//	{
	//		return 0;
	//	}
	//}
	//Tsubstrate straight_substrate(300, 300, 15, 0, 300, 40, 0);
	//straight_substrate.output_substrate(output_tec,0);
	//straight_substrate.generate_coupled_straight_substrate(0, coupled_substrate, is_success);
	Tsubstrate wavy_substrate(500,500,15,0,150,0.23,0.2,20);
	wavy_substrate.output_substrate(output_tec,0);
	wavy_substrate.generate_coupled_wavy_substrate(300, 300, 15, 0, 150, 0.23, 0.2, 20, coupled_substrate, is_success);
	coupled_substrate.output_substrate(output_tec, 1);
	//cout << straight_straight_coupled(0.29, 15, 0, 150, 0, 0) << endl;
	//cout << wavy_wavy_coupled(0.2, 15, 0, 150, 0.2, 0.2) << endl;
	system("Pause");
	return 0;
}
double straight_straight_coupled(double density, double r_NW, double r_P, double l_NW, double theta1, double theta2)
{
	int N_test=50;
	int N_success=0;
	#pragma omp parallel for reduction(+: N_success) schedule(dynamic)
	for (int i = 0; i < N_test; i++)
	{
		Tsubstrate straight_substrate(400,400,r_NW,r_P,l_NW,density,theta1);
		Tsubstrate coupled_substrate;
		bool is_successed = false;
		straight_substrate.generate_coupled_straight_substrate(300, 300, r_NW, r_P, l_NW, density, theta2, coupled_substrate, is_successed);
		if (is_successed)
		{
			N_success = N_success + 1;
		}
	}
	return 1.0*N_success/N_test;
}
double wavy_wavy_coupled(double density, double r_NW, double r_P, double l_NW, double curvature1, double curvature2)
{
	int n_divided = 20;
	int N_test = 50;
	int N_success = 0;
	#pragma omp parallel for reduction(+: N_success) schedule(dynamic)
	for (int i = 0; i < N_test; i++)
	{
		
		Tsubstrate straight_substrate(500, 500, r_NW, r_P, l_NW, density, curvature1,n_divided);
		Tsubstrate coupled_substrate;
		bool is_successed = false;
		straight_substrate.generate_coupled_wavy_substrate(300, 300, r_NW, r_P, l_NW, density, curvature2, n_divided,coupled_substrate, is_successed);
		if (is_successed)
		{
			N_success = N_success + 1;
		}
	}
	return 1.0*N_success/N_test;
}
double straight_wavy_coupled(double density, double r_NW, double r_P, double l_NW, double curvature, double theta)
{
	int n_divided = 20;
	int N_test = 50;
	int N_success = 0;
	#pragma omp parallel for reduction(+: N_success) schedule(dynamic)
	for (int i = 0; i < N_test; i++)
	{
		Tsubstrate straight_substrate(500, 500, r_NW, r_P, l_NW, density, curvature, n_divided);
		Tsubstrate coupled_substrate;
		bool is_successed = false;
		straight_substrate.generate_coupled_straight_substrate(300, 300, r_NW, r_P, l_NW, density, theta, coupled_substrate, is_successed);
		if (is_successed)
		{
			N_success = N_success + 1;
		}
	}
	return 1.0*N_success/N_test;
}