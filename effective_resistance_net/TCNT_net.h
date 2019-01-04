#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "data_structure.h"
#include "functions.h"
#include "TResistance_net.h"
#include "readin.h"
using namespace std;
struct SCNT
{
	SCNT(int mm, int nn) { m = mm; n = nn; };
	vector<vec3D> node_list;
	vector<vector<int>> IEN_list;
	int nump, nume;
	int m, n;

	vector<vec3D> calculate_centroid_line();
};
struct Sresistance_T
{
	Sresistance_T(Selectrical_node n1, Selectrical_node n2, double RR)
	{
		node1 = n1; node2 = n2; R = RR;
	}
	Selectrical_node node1;
	Selectrical_node node2;
	double R;
};
struct Sresistance
{
	Sresistance(int ii, int jj, double RR)
	{
		i = ii; j = jj; R = RR;
	}
	int i, j;
	double R;
};
void calculate_moment(vec3D &p1, vec3D &p2, vec3D &p3, vec3D &p4, vec3D &moment, double &area);
//Calculate the area and moment of a quadrilateral
double triagnle_area(vec3D &p1, vec3D &p2, vec3D &p3);
//Calculate the area of a triangle
int identify_RVE_face(string face_name);
//Identify the face number of a face_name
class TCNT_net
{
public:
	TCNT_net(double dE, int M, double sigma_CNT):_dE(dE), _M(M), _sigma_CNT(sigma_CNT) {};

	void input_CNT_net();
	void input_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> bc_Info);
	//Input the CNT net and the bc information
	double Calculate_effective_resistance(string positive_face,string negative_face,int output_control=0);
	//positive_face,negative_face="x(y,z)_min" and "x(y,z)_max"
	//positive_face: The RVE face connected to positive polar of the source
	//negative_face: The RVE face connected to negative polar of the source
	//is_open: Whether no positive or negative node in the RVE
private:
	vector<vector<vec3D>> _CNT_net;
	vector<vector<int>> _bc_info;
	vector<Sresistance_T> _RT_list;
	vector<Sresistance> _resistance_net;
	int _num_CNT;
	int _num_tunneling_R;
	double _dE,_sigma_CNT;
	int _M;

	void Create_tunneling_resistance();
	void Create_intrinsic_resistance();
	void Intrinsic_resistance_on_CNT(int CNT_ID);
};