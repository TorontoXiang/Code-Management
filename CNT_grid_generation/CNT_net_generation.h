#pragma once
#include "TCNT_grid.h"
struct Slink_node
{
	Slink_node(vec3D coor,int State=-3);
	//Create a new node

	vec3D pos;
	double sign_dis;
	int state;
	Slink_node* next_node;
	Slink_node* pre_node;
};
class TCNT_net_generator
{
public:
	TCNT_net_generator() { /*srand((int)time(0));*/ };

	void input_RVE_info(ifstream &input);
	//Input the RVE and CNT information

	void Create_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info,int output_control=0);
	//Create a CNT network

	void Generate_CNT_net_grid(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info, ofstream &output_k, ofstream &output_tec, ofstream &output_bc);
	//Output the k file, Tecplot file and boundary information of a CNT network

	void Clipping_curved_CNT(vector<vec3D> &curved_CNT, vector<vector<vec3D>> &sub_CNT, vector<vector<int>> &bc_info);
	//Generate the subCNTs inside the RVE form a curved CNT

	int G_type() { return _type; };
	void G_RVE_size(double &lx, double &ly, double &lz) { lx = 1e-9*_lx, ly = 1e-9*_ly, lz = 1e-9*_lz; };
	void S_fraciton(double fraction) { _volume_fraction = fraction; };

private:
	int _m, _n;                 //The division on the truncated face of CNT
	int _n_divided;             //The division along the axis of CNT
	int _type;                  //0-straight CNT network; 1-curved CNT network
	double _curvature_ratio;    //Describe the curvature of a curved CNT network
	double _lx, _ly, _lz;       //The size of the RVE
	double _d_CNT, _l_CNT;      //The size of each CNT
	double _d_VDW;              //The van der Waals distance of CNT
	double _E, _mu;             //The Young's moduli and Possion's ratio of CNT
	double _volume_fraction;    //The volume fraction of CNT in the composite

	void Create_straight_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info,int output_control=0);
	//Create a straight CNT network (CNT_net) and the boundary information of each CNT (bc_Info)

	void Create_wavy_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info,int output_control=0);
	//Create a wavy CNT network (CNT_net) and the boundary information of each CNT (bc_Info)

	int segment_state(vec3D p1, vec3D p2);
	//Determine the state of a segment

	bool is_inside_RVE(vec3D p);
	//Whether p is inside the RVE
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