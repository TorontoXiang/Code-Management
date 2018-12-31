#pragma once
#include "TCNT_grid.h"
class TCNT_net_generator
{
public:
	TCNT_net_generator() {};

	void input_RVE_info(ifstream &input);
	//Input the RVE and CNT information

	void Create_straight_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info);
	//Create a straight CNT network (CNT_net) and the boundary information of each CNT (bc_Info)

	void Create_wavy_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info);
	//Create a wavy CNT network (CNT_net) and the boundary information of each CNT (bc_Info)

	void Generate_CNT_net_grid(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info, ofstream &output_k, ofstream &output_tec, ofstream &output_bc);
	//Output the k file, Tecplot file and boundary information of a CNT network

	int _m, _n;                 //The division on the truncated face of CNT
	int _n_divided;             //The division along the axis of CNT
	int _type;                  //0-straight CNT network; 1-curved CNT network
	double _curvature_ratio;    //Describe the curvature of a curved CNT network
	double _lx, _ly, _lz;       //The size of the RVE
	double _d_CNT, _l_CNT;      //The size of each CNT
	double _d_VDW;              //The van der Waals distance of CNT
	double _E, _mu;             //The Young's moduli and Possion's ratio of CNT
	double _volume_fraction;    //The volume fraction of CNT in the composite
};
void Create_straight_CNT_net(int num_CNT, double l_CNT, double lx, double ly, double lz);
//Create a random straight CNT net
//lx,ly,lz - The length of the box
//num_CNT - The number of the CNT put in the box
//l_CNT - The length of the CNT
void Create_wavy_CNT_net(int num_CNT, double l_CNT, int n_divided, double ratio, double lx, double ly, double lz);
//Create a random wavy CNT net
//lx,ly,lz - The length of the box
//num_CNT - The number of the CNT put in the box
//l_CNT - The length of the CNT
//n_divided - The division of each CNT
//ratio - control the curvature of the CNT
//E,mu - The material property of the CNT
void Generate_straight_CNT(vec3D p1, vec3D p2, int n_divided, double r, double l, int m, int n, ofstream& output_k, ofstream& output_tec, double E, double mu);
//Generate a straight CNT grid and input file if the two endpoints are given
//p1,p2 - The two endpoints
//n_divided - The number of the division along the length
//r,l - The radius and length of the CNT
//m,n - The parameters for generating the truncation mesh
//output_k,output_tec - The output of the CNT grid
//E,mu - The material property of the CNT
void Generate_CNT_net_grid(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info, int m, int n, double E, double mu, ofstream &output_k, ofstream &output_tec, ofstream &output_bc);
//Output the k file and Tecplot file of a CNT net