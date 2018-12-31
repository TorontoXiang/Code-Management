#ifndef DATA_STRUCTURE
#define DATA_STRUCTURE
#include <vector>
#include <fstream>
using namespace std;
double const pi = 3.141592653;
struct Saxis_load
{
	int direction;
	double value;
};
struct Slocation_info
{
	int cell_id;        //The cell_id where a node locates
	double iso_coor[3]; //The isoparametric coordinate of the node in the cell
};
struct Sstress
{
	double sxx, syy, szz, sxy, syz, sxz;
};
struct Snode_stress
{
	Snode_stress() { num = 0; sxx = syy = szz = sxy = syz = sxz = 0; };
	int num;
	double sxx, syy, szz, sxy, syz, sxz;

	void Average() 
	{
		if (num>0)
		{
			sxx = sxx / num; syy = syy / num; szz = szz / num; 
			sxy = sxy / num; sxz = sxz / num; syz = syz / num;
		}
	};
	void add(Sstress const &other)
	{
		sxx = sxx + other.sxx; syy = syy + other.syy; szz = szz + other.szz;
		sxy = sxy + other.sxy; sxz = sxz + other.sxz; syz = syz + other.syz;
	};
	double equivalent_stress()
	{
		double J2 = 0.5*(sxx*sxx + syy * syy + szz * szz) + syz * syz + sxz * sxz + sxy * sxy;
		return sqrt(J2 * 3);
	};
};
struct Smat
{
	double _E, _mu;            //Young's modulus and Poission's ratio
};
struct vec3D
{
	double x, y, z;

	vec3D();
	vec3D(double xx, double yy, double zz) { x = xx, y = yy, z = zz; };
	//Normalize a vec3D
	void normalize();
	//Value a vec3D
	void value(double const xx, double const yy, double const zz);
	//Value the ith variable 
	void value(int i, double num);
	//Self multiply
	double self_multuply() const;
	//Get length
	double get_length() const;
	//Multiply by a matrix
	vec3D multiply_by_matrix(double(&A)[3][3]);
	//Transfer the point into a new coordinate system
	vec3D multiply_by_matrix_transpose(double(&A)[3][3]);
	vec3D transfer_into_new_system(vec3D &n1, vec3D &n2, vec3D &n3, vec3D &x0);
	//Get ith variable
	double access(int i);

	//Opreations
	bool operator ==(vec3D const &other);
	vec3D operator +(vec3D const &other);
	vec3D operator -(vec3D const &other);
	vec3D operator -();
	vec3D operator *(double const a);
	vec3D operator /(double const a);
	vec3D operator %(vec3D const &other);    //vector product
	double operator *(vec3D const &other);  //dot product
	bool operator <(vec3D const &other);
	bool operator >(vec3D const &other);
};
struct link_node
{
	link_node(vec3D coor);
	//Create a new node

	vec3D pos;
	int state;
	link_node* next_node;
	link_node* pre_node;
};
void inverse(double a[][3], double &det, bool &degeneratre);
//The inverse of a 3*3 matrix
int local_node_id(int i, int j, int k);
template<class T> T maxval(T a, T b)
{
	if (a >= b)
	{
		return a;
	}
	return b;
}

template<class T> T minval(T a, T b)
{
	if (a >= b)
	{
		return b;
	}
	return a;
}
#endif