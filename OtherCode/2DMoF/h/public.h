#ifndef PUBLIC
#define PUBLIC
#include <cmath>
#include <vector>
using namespace std;
double const pi=3.141592653589793;
struct vec2D
{
	vec2D():x(0),y(0){};
	vec2D(double xx,double yy):x(xx),y(yy){};
	double x,y;

	void value(double xx,double yy){x=xx;y=yy;};
	void normalize(){double length=sqrt(x*x+y*y);x=x/length;y=y/length;};
	double operator *(vec2D const &other){return x*other.x+y*other.y;};         
	vec2D operator *(double const &num){vec2D temp(num*x,num*y);return temp;};
	vec2D operator +(vec2D const &other){vec2D temp(other.x+x,other.y+y);return temp;};
	vec2D operator -(vec2D const &other){vec2D temp(x-other.x,y-other.y);return temp;};
	vec2D operator /(double const &num){vec2D temp(x/num,y/num);return temp;};
};
struct discontinuous_point
{
	double _theta;
	double _length;
	double _xc,_yc;
	double _d1;
};
struct iteration_point
{
	iteration_point():theta(0),d1(0),d2(0),d3(0){};
	double theta;
	double value,d1,d2,d3;
	
	void calculate_d2(discontinuous_point& dp,double s1,double s2,double c1,double c2,vec2D& centroid_ref,double area_below);
};

double minval(double a,double b);
double maxval(double a,double b);

double area(vec2D& p1,vec2D& p2,vec2D& p3);
//Calculate the area of a triagnle with p1,p2,p3
double calculate_inclination(vec2D& dp);
//Calculte the inclination of dp
double root_of_cubic_interpolation(iteration_point& iteration_min,iteration_point& iteration_max);
//Calcualte the root of a cubic function in the feasible region
double root_of_fifth_interpolation(iteration_point& iteration_min,iteration_point& iteration_max);
bool is_intersect(vec2D x_range_old,vec2D y_range_old,vec2D x_range_new,vec2D y_range_new);

vec2D calculate_moment(double xl,double yl,double xr,double yr);
#endif 
