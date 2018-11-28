#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include "data_structure.h"
#include "functions.h"
#include "TCNT_grid.h"
using namespace std;
class Tsubstrate
{
public:
	Tsubstrate(){};
	Tsubstrate(double lx, double ly, double r_NW, double r_P, double l_NW,int num_NW, double curvature, int n_divide);
	//Generate a substrate with wavy nanowires
	//lx,ly - the size of the substrate
	//r_NW,r_P - the raduis of nanowire and the pattern
	//curvature - the curvature of the nanowire
	//n_divide - the number of segment in each nanowire
	Tsubstrate(double lx, double ly, double r_NW, double r_P, double l_NW, int num_NW, double theta_max);
	//Generate a substrate with straight nanowires
	//lx,ly - the size of the substrate
	//r_NW,r_P - the raduis of nanowire and the pattern
	//theta_max - the orientation of each nanowire

	void generate_coupled_straight_substrate(double theta_max,Tsubstrate &coupled_substrate,bool &is_successed);
	//Generate a coupled substrate with straight NW
	void generate_coupled_wavy_substrate(double curvature, Tsubstrate &coupled_substrate, bool &is_successed);
	//Generate a coupled substrate with wavy NW
	//Generate a coupled substrate with straight NW
	void output_substrate(ofstream &output_tec,int value);
	//Output the substrate

private:
	vector<vector<vec3D>> _NW_list;
	vector<vec3D> _root_list;
	double _lx, _ly;
	double _r_NW,_l_NW;
	double _r_P;
	int _num_NW, _n_divide;

	bool is_root_valid(vec3D x0);
	//Determine whether the root of a new NW is valid (not inside the pattern region)
	double max_z();
	//Calculate the maximal z coordinate of all NWs
	vector<vec3D> generate_valid_straight_NW(double z, int direction, double theta_max);
	//Generate a straight nanowire which is validated for this substrate
	vector<vec3D> generate_valid_wavy_NW(double z, int direction, double curvature,int n_divide);
	//Generate a straight nanowire which is validated for this substrate
};
