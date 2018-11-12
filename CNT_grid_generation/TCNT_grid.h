#ifndef TCNT_GRID
#define TCNT_GRID
#include "data_structure.h"
#include <vector>
#include <fstream>
using namespace std;
class Ttruncation
{
public:
	Ttruncation() {};

	void generate_truncation_grid(double r, double l,int m, int n);
private:
	vec3D _r0;
	double _A[3][3];
	double _c_t;
	int _nump;
	int _nume;
	vec3D* _node_list;
	int** _IEN_list;

	friend class TCNT_grid;
};

class TCNT_grid
{
public:
	TCNT_grid(vector<vec3D> node_list);
	void generate_CNT_grid(double r, double l, double m, double n);
	//r:radiu of CNT, l:inner square length
	//m:rudial direction division, n:ring direction division
	void output_CNT_tecplot(ofstream& output);
	void output_CNT_k_file(ofstream& output,double E,double mu);
private:
	int _num_node;
	vector<Ttruncation> _truncation_list;
};

void T_matrix(vec3D z, double(&T)[3][3]);
//Calculate the Transformation matrix between two coordinate system 1 and 2
//z: the coordinate of z-axis of system 2 in system 1
void matrix_multiply(double(&A1)[3][3], double(&A2)[3][3], double(&A)[3][3]);
int calculate_index(int i, int j, int i_max2, int i_max1, int j_max1, int nump1);
#endif