#include "data_structure.h"
#include <iostream>
#include <cmath>
using namespace std;
//Operations for vec3D
vec3D::vec3D()
{
	x=y=z=0;
	return;
}
void vec3D::normalize()
{
	double l;
	l=get_length();
	if (l<1e-20)
	{
		cout<<"Warning: zero vector can not be normalized!!"<<endl;
	}
	x=x/l;y=y/l;z=z/l;
	return;
}
void vec3D::value(double const xx,double const yy,double const zz)
{
	x=xx;y=yy;z=zz;
	return;
}
void vec3D::value(int i,double num)
{
	if (i==1)
	{
		x=num;
	}
	else if (i==2)
	{
		y=num;
	}
	else if (i==3)
	{
		z=num;
	}
	else
	{
		cout<<"Error: the subscript of vec3D is "<<i<<" in calling vec3D::value()"<<endl;
		system("Pause");
		exit(0);
	}
	return;
}
double vec3D::self_multuply() const
{
	return x*x+y*y+z*z;
}
double vec3D::get_length() const
{
	return sqrt(self_multuply());
}
vec3D vec3D::multiply_by_matrix(double (&A)[3][3])
{
	vec3D temp;
	temp.x=A[0][0]*x+A[0][1]*y+A[0][2]*z;
	temp.y=A[1][0]*x+A[1][1]*y+A[1][2]*z;
	temp.z=A[2][0]*x+A[2][1]*y+A[2][2]*z;
	return temp;
}
vec3D vec3D::transfer_into_new_system(vec3D &n1,vec3D &n2,vec3D &n3,vec3D &x0)
{
	vec3D diff(x-x0.x,y-x0.y,z-x0.z);
	vec3D temp;
	temp.x=diff*n1;temp.y=diff*n2;temp.z=diff*n3;
	return temp;
}
double vec3D::access(int i)
{
	if (i==1)
	{
		return x;
	}
	else if (i==2)
	{
		return y;
	}
	else if (i==3)
	{
		return z;
	}
	else
	{
		cout<<"Error:The subscript of vec3D is "<<i<<" in calling vec3D::access()"<<endl;
		system("Pause");
		exit(0);
	}
}
bool vec3D::operator ==(vec3D const &other)
{
	double epsilon=1e-10;
	if (abs(x-other.x)+abs(y-other.y)+abs(z-other.z)<epsilon)
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
	temp.x=x+other.x;temp.y=y+other.y;temp.z=z+other.z;
	return temp;
}
vec3D vec3D::operator -(vec3D const &other)
{
	vec3D temp;
	temp.x=x-other.x;temp.y=y-other.y;temp.z=z-other.z;
	return temp;
}
vec3D vec3D::operator -()
{
	vec3D temp(-x,-y,-z);
	return temp;
}
vec3D vec3D::operator *(double const a)
{
	vec3D temp;
	temp.x=a*x;temp.y=a*y;temp.z=a*z;
	return temp;
}
vec3D vec3D::operator /(double const a)
{
	vec3D temp;
	if (abs(a)==0)
	{
		cout<<"Warning: The denominator is zero in calculating / for vec3D"<<endl;
		system("Pause");
	}
	temp.x=x/a;temp.y=y/a;temp.z=z/a;
	return temp;
}
vec3D vec3D::operator %(vec3D const &other)
{
	vec3D temp;
	temp.x=y*other.z-z*other.y;temp.y=-(x*other.z-z*other.x);temp.z=x*other.y-y*other.x;
	return temp;
}
double vec3D::operator *(vec3D const &other)
{
	return x*other.x+y*other.y+z*other.z;
}
bool vec3D::operator <(vec3D const &other)
{
	if (x<other.x && y<other.y && z<other.z)
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
	if (x>other.x && y>other.y && z>other.z)
	{
		return true;
	}
	else
	{
		return false;
	}
}
//Struct matrix2D:
matrix2D matrix2D::operator +(matrix2D other)
{
	matrix2D temp;
	temp.a11=a11+other.a11;temp.a12=a12+other.a12;
	temp.a21=a21+other.a21;temp.a22=a22+other.a22;
	return temp;
}
matrix2D matrix2D::operator *(matrix2D other)
{
	matrix2D temp;
	temp.a11=a11*other.a11+a12*other.a21;
	temp.a12=a11*other.a12+a12*other.a22;
	temp.a21=a21*other.a11+a22*other.a21;
	temp.a22=a21*other.a12+a22*other.a22;
	return temp;
}
double matrix2D::operator %(vec2D v)
{
	return a11*v.x*v.x+(a12+a21)*v.x*v.y+a22*v.y*v.y;
}
matrix2D matrix2D::operator/(double a)
{
	matrix2D temp;
	temp.a11=a11/a;temp.a12=a12/a;temp.a21=a21/a;temp.a22=a22/a;
	return temp;
}
matrix2D matrix2D::transfer()
{
	matrix2D temp;
	temp.a11=a11;temp.a22=a22;
	temp.a12=a21;temp.a21=a12;
	return temp;
}
vec2D matrix2D::operator *(vec2D v)
{
	vec2D temp;
	temp.x=a11*v.x+a12*v.y;
	temp.y=a21*v.x+a22*v.y;
	return temp;
}
matrix2D matrix2D::operator*( double a )
{
	matrix2D temp;
	temp.a11=a11*a;temp.a12=a12*a;temp.a21=a21*a;temp.a22=a22*a;
	return temp;
}
matrix2D matrix2D::operator-( matrix2D other )
{
	matrix2D temp;
	temp.a11=a11-other.a11;temp.a12=a12-other.a12;
	temp.a21=a21-other.a21;temp.a22=a22-other.a22;
	return temp;
}
//Struct vec2D
double vec2D::operator *(vec2D other)
{
	return x*other.x+y*other.y;
}
vec2D vec2D::operator *(double a)
{
	vec2D temp;
	temp.x=a*x;
	temp.y=a*y;
	return temp;
}
vec2D vec2D::operator -(vec2D other)
{
	vec2D temp;
	temp.x=x-other.x;temp.y=y-other.y;
	return temp;
}
vec2D vec2D::operator/(double a)
{
	vec2D temp;
	temp.x=x/a;temp.y=y/a;
	return temp;
}
double vec2D::get_length()
{
	return sqrt(x*x+y*y);
}
