#ifndef PUBLIC_FUNCTION
#define PUBLIC_FUNCTION
#include "data_structure.h"
#include <string>
#include <vector>
#include <iostream>
using namespace std;
vec3D gauss_for_centroid(vec3D &p1,vec3D &p2,vec3D &p3);        //Gauss integration for centroid calculation
vec3D gauss_for_MoF_derivative(vec3D &p1,vec3D &p2,vec3D &p3);  //Gauss integration for the derivative in MoF method
double gauss(vec3D &p1,vec3D &p2,vec3D &p3,string type="z");    //Gauss integration on a triangle piece
void area(vec3D &p1,vec3D &p2,vec3D &p3,double &s,double &nz);  //Calculate the area of a triangle and the third component of the unit outward normal
                                                                //The direction of the triangle is p1->p2->p3
void area(vec3D &p1,vec3D &p2,vec3D &p3,double &s,vec3D &n);    //Calculate the area of a triangle and its unit outward normal
                                                                //The direction of the triangle is p1->p2->p3
template<class T> T maxval(T a,T b);
//Return max of a and b
template<class T> T minval(T a,T b);
//Return min of a and b
double random(double a,double b);
//Return a random number between a and b
void inverse(double a[][3],double &det,bool &degeneratre);
//The inverse of a 3*3 matrix
double tetredron_volume(vec3D &p1,vec3D &p2,vec3D &p3,vec3D &p4);
//Calculate the volume of a tetrahedron
void absmax_component(vec3D &n,int &id,double &absmax);
//Calculate the abs maximum component of n
int sign(double a);
//Calculate the sign of a
bool comparison_cell_coordinate(vec3D const &coor_min_old,vec3D const &coor_max_old,vec3D const &coor_min_new,vec3D const &coor_max_new);
//Campare the coordinate to determine whether intersected
vec3D interpolate(vec3D_double &p1,vec3D_double &p2);
//Interpolate the middle point between two points with different sigh distance
bool is_exist(int i,vector<int> &list);
//Whether i is in list
vec3D convert_to_sphere_coor(vec3D coor);
int find_int(double number);
double calcualte_alpha(double d0,int N,double r_in,double r_out);
//Calculate the index for implosion node moving
double f_calculate_alpha(int N,double k,double x);
//-------------------------------------------------------------------------------
//Functions for optimization
//-------------------------------------------------------------------------------
void bisect(double& a,double& b,double& f);
//Bisect method
//a:variable with objective<0
//b:variable with objective>0
//f:objective at (a+b)/2
//Return a,b in next iteration
void secant(double& a,double& b,double& fa,double& fb,double& f);
//Secant method
//a:variable with fa<0
//b:variable with fb>0
//f:objective at the next iteration point
//Update a or b by
double newton(double& x,double& f,double& df);
//Newton method
//Input the former information and output the next point
//Functions used in BFGS method
double quad(double a1,double a2,double f1,double df1,double f2);
//Quadratic interpolation
//Input the former information and output the next iteration point
double cubic1(double a1,double a2,double f1,double df1,double f2,double df2);
//Cubic interpolation if a2 exists
//Input the former information and output the next iteration point
double cubic2(double a1,double a,double f1,double df1,double f,double df);
//Cubic interpolation if a2 is infinite
//Input the former information and output the next iteration point
double root_of_approximate_volume_function(vec3D &min_info,vec3D &max_info,double &refence_fraction);
//Calcualte the root of a cubic function in [d_min,d_max]
//min_info:the information at d_min including d_min,f_min,df_min
//max_info:the information at d_max including d_max,f_max,df_max
//--------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
//Template functions' realize
//---------------------------------------------------------------------------------------
template<class T> T maxval(T a,T b)
{
	if (a>=b)
	{
		return a;
	}
	return b;
}

template<class T> T minval(T a,T b )
{
	if (a>=b)
	{
		return b;
	}
	return a;
}
//------------------------------------------------------------------------------------------
#endif