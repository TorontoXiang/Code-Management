#include "data_structure.h"
#include <iostream>
using namespace std;
void Snode::calculate_location(double(&x_min)[3], double(&x_max)[3], double(&interval)[3], int nx_max, int ny_max, int nz_max)
{
	if (_is_surface)
	{
		if (_pos[0]<x_min[0] || _pos[1] < x_min[1] || _pos[2] < x_min[2] || _pos[0] > x_max[0] || _pos[1] > x_max[1] || _pos[2] > x_max[2])
		{
			cout << _pos[0] << " " << _pos[1] << " " << _pos[2] << endl;
			cout << "CNT node is out of the polymer" << endl;
			system("Pause");
			//return;
			//exit(0);
		}
		_location = new Slocation_info;
		double delta[3];
		delta[0] = _pos[0] - x_min[0]; delta[1] = _pos[1] - x_min[1]; delta[2] = _pos[2] - x_min[2];
		int nx = (int)floor(delta[0] / interval[0]);
		int ny = (int)floor(delta[1] / interval[1]);
		int nz = (int)floor(delta[2] / interval[2]);
		nx = minval(nx, nx_max - 1); ny = minval(ny, ny_max - 1); nz = minval(nz, nz_max - 1);
		_location->cell_id= ny_max * nz_max*nx + nz_max * ny + nz;
		double dx = _pos[0] - nx * interval[0];
		_location->iso_coor[0] = -1 + 2 * dx / interval[0];
		double dy = _pos[1] - ny * interval[1];
		_location->iso_coor[1] = -1 + 2 * dy / interval[1];
		double dz = _pos[2] - nz * interval[2];
		_location->iso_coor[2] = -1 + 2 * dz / interval[2];
	}
	return;
}
vec3D::vec3D()
{
	x = y = z = 0;
	return;
}
void vec3D::normalize()
{
	double l;
	l = get_length();
	if (l < 1e-20)
	{
		cout << "Warning: zero vector can not be normalized!!" << endl;
	}
	x = x / l; y = y / l; z = z / l;
	return;
}
void vec3D::value(double const xx, double const yy, double const zz)
{
	x = xx; y = yy; z = zz;
	return;
}
void vec3D::value(int i, double num)
{
	if (i == 1)
	{
		x = num;
	}
	else if (i == 2)
	{
		y = num;
	}
	else if (i == 3)
	{
		z = num;
	}
	else
	{
		cout << "Error: the subscript of vec3D is " << i << " in calling vec3D::value()" << endl;
		system("Pause");
		exit(0);
	}
	return;
}
double vec3D::self_multuply() const
{
	return x * x + y * y + z * z;
}
double vec3D::get_length() const
{
	return sqrt(self_multuply());
}
vec3D vec3D::multiply_by_matrix(double(&A)[3][3])
{
	vec3D temp;
	temp.x = A[0][0] * x + A[0][1] * y + A[0][2] * z;
	temp.y = A[1][0] * x + A[1][1] * y + A[1][2] * z;
	temp.z = A[2][0] * x + A[2][1] * y + A[2][2] * z;
	return temp;
}
vec3D vec3D::multiply_by_matrix_transpose(double(&A)[3][3])
{
	vec3D temp;
	temp.x = A[0][0] * x + A[1][0] * y + A[2][0] * z;
	temp.y = A[0][1] * x + A[1][1] * y + A[2][1] * z;
	temp.z = A[0][2] * x + A[1][2] * y + A[2][2] * z;
	return temp;
}
vec3D vec3D::transfer_into_new_system(vec3D &n1, vec3D &n2, vec3D &n3, vec3D &x0)
{
	vec3D diff(x - x0.x, y - x0.y, z - x0.z);
	vec3D temp;
	temp.x = diff * n1; temp.y = diff * n2; temp.z = diff * n3;
	return temp;
}
double vec3D::access(int i)
{
	if (i == 1)
	{
		return x;
	}
	else if (i == 2)
	{
		return y;
	}
	else if (i == 3)
	{
		return z;
	}
	else
	{
		cout << "Error:The subscript of vec3D is " << i << " in calling vec3D::access()" << endl;
		system("Pause");
		exit(0);
	}
}
bool vec3D::operator ==(vec3D const &other)
{
	double epsilon = 1e-10;
	if (abs(x - other.x) + abs(y - other.y) + abs(z - other.z) < epsilon)
	{
		return true;
	}
	else
	{
		return false;
	}
}
vec3D vec3D::operator +(vec3D const &other)
{
	vec3D temp;
	temp.x = x + other.x; temp.y = y + other.y; temp.z = z + other.z;
	return temp;
}
vec3D vec3D::operator -(vec3D const &other)
{
	vec3D temp;
	temp.x = x - other.x; temp.y = y - other.y; temp.z = z - other.z;
	return temp;
}
vec3D vec3D::operator -()
{
	vec3D temp(-x, -y, -z);
	return temp;
}
vec3D vec3D::operator *(double const a)
{
	vec3D temp;
	temp.x = a * x; temp.y = a * y; temp.z = a * z;
	return temp;
}
vec3D vec3D::operator /(double const a)
{
	vec3D temp;
	if (abs(a) == 0)
	{
		cout << "Warning: The denominator is zero in calculating / for vec3D" << endl;
		system("Pause");
	}
	temp.x = x / a; temp.y = y / a; temp.z = z / a;
	return temp;
}
vec3D vec3D::operator %(vec3D const &other)
{
	vec3D temp;
	temp.x = y * other.z - z * other.y; temp.y = -(x*other.z - z * other.x); temp.z = x * other.y - y * other.x;
	return temp;
}
double vec3D::operator *(vec3D const &other)
{
	return x * other.x + y * other.y + z * other.z;
}
bool vec3D::operator <(vec3D const &other)
{
	if (x < other.x && y < other.y && z < other.z)
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool vec3D::operator >(vec3D const &other)
{
	if (x > other.x && y > other.y && z > other.z)
	{
		return true;
	}
	else
	{
		return false;
	}
}
void inverse(double a[][3], double& det, bool& degenerate)
{
	double inva[3][3];
	degenerate = false;
	det = a[0][0] * a[1][1] * a[2][2] - a[0][0] * a[1][2] * a[2][1] - a[0][1] * a[1][0] * a[2][2] + a[0][1] * a[1][2] * a[2][0] + a[0][2] * a[1][0] * a[2][1] - a[0][2] * a[1][1] * a[2][0];
	if (abs(det) < 1e-20)
	{
		degenerate = true;
		return;
	}
	inva[0][0] = a[1][1] * a[2][2] - a[1][2] * a[2][1]; inva[0][1] = a[0][2] * a[2][1] - a[0][1] * a[2][2]; inva[0][2] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
	inva[1][0] = a[1][2] * a[2][0] - a[1][0] * a[2][2]; inva[1][1] = a[0][0] * a[2][2] - a[0][2] * a[2][0]; inva[1][2] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
	inva[2][0] = a[1][0] * a[2][1] - a[1][1] * a[2][0]; inva[2][1] = a[0][1] * a[2][0] - a[0][0] * a[2][1]; inva[2][2] = a[0][0] * a[1][1] - a[0][1] * a[1][0];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			a[i][j] = inva[i][j] / det;
		}
	}
	return;
}
int local_node_id(int i, int j, int k)
{
	if (i==0 && j==0)
	{
		return 4 * k;
	}
	else if (i==1 && j==0)
	{
		return 1 + 4 * k;
	}
	else if (i==1 && j==1)
	{
		return 2 + 4 * k;
	}
	else 
	{
		return 3 + 4 * k;
	}
}