#ifndef DATA_STRUCTURE
#define DATA_STRUCTURE
#include <vector>
#include <fstream>
using namespace std;
struct Slocation_info
{
	int cell_id;        //The cell_id where a node locates
	double iso_coor[3]; //The isoparametric coordinate of the node in the cell
};
struct Snode
{
	Snode() {};
	Snode(int id, double(&coor)[3]) 
	{ 
		_id = id; _pos[0] = coor[0]; _pos[1] = coor[1]; _pos[2] = coor[2]; 
		_bc[0] = _bc[1] = _bc[2] = 0; _bc_value[0] = _bc_value[1] = _bc_value[2] = 0;
		_is_surface = false;
	};

	int _id;
	double _pos[3];
	int _bc[3];                //Boundary condition  of a node
	double _bc_value[3];       //The value of boundary condition
	bool _is_surface;          //Whether a node is a surface node
	Slocation_info* _location; //The location information of a CNT node (Only allocate memory for the surface node in CNT grid)

	void calculate_location(double(&x_min)[3], double(&x_max)[3], double(&interval)[3],int nx_max,int ny_max,int nz_max);
	//Calculate the location of a surface node in CNT grid
	
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
