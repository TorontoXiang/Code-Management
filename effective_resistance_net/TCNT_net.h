#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "data_structure.h"
#include "functions.h"
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
void new_line(ifstream& input);     
//Move to the beginning of the next linee
void calculate_moment(vec3D &p1, vec3D &p2, vec3D &p3, vec3D &p4, vec3D &moment, double &area);
//Calculate the area and moment of a quadrilateral
double triagnle_area(vec3D &p1, vec3D &p2, vec3D &p3);
//Calculate the area of a triangle
class TCNT_net
{
public:
	TCNT_net(double d_tunneling, int M, double r_CNT):_d_tunneling(d_tunneling), _M(M), _r_CNT(r_CNT) {};

	void input_CNT_net();
	//Input the CNT net and the bc information
	//void Create_resistance_net();
	void Create_tunneling_resistance();
	void Create_intrinsic_resistance();
	void Intrinsic_resistance_on_CNT(int CNT_ID);
private:
	vector<vector<vec3D>> _CNT_net;
	vector<vector<int>> _bc_info;
	vector<Sresistance_T> _RT_list;
	vector<Sresistance> _resistance_net;
	int _num_CNT;
	int _num_tunneling_R;
	double _d_tunneling,_r_CNT;
	int _M;


	//void Create_tunneling_resistance();
	//void Create_intrinsic_resistance();
	//void Intrinsic_resistance_on_CNT(int CNT_ID);
};