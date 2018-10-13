#include "public.h"
#include <cmath>
#include <iostream>
using namespace std;
void iteration_point::calculate_d2(discontinuous_point& dp,double s1,double s2,double c1,double c2,vec2D& centroid_ref,double area_below)
{
	double length=dp._length;
	double temp=length*length*length/(12*area_below);
	theta=dp._theta;
	d1=dp._d1;
	double s=sin(theta),c=cos(theta);
	double ss1=s1*c-c1*s,ss2=s2*c-c2*s;
	double cc1=c1*c+s1*s,cc2=c2*c+s2*s;
	double a=ss1/cc1+ss2/cc2;
	vec2D x2(-3*a*s+2*c,3*a*c+2*s);
	d2=2*temp*temp+((dp._xc-centroid_ref.x)*(-3*a*s+2*c)+(dp._yc-centroid_ref.y)*(3*a*c+2*s))*temp;
	//double b=-(1/(cc1*cc1)+1/(cc2*cc2));
	//d3=-3*temp*temp*(-s*(-3*a*s+2*c)+c*(3*a*c+2*s))+((dp._xc-centroid_ref.x)*(-1.5*a*(-3*a*s+2*c)-3*(b*s+a*c)-2*s)+(dp._yc-centroid_ref.y)*(-1.5*a*(3*a*c+2*s)+3*(b*c-a*s)+2*c))*temp;
	//d2=2*temp*temp+((dp._centroid_below-centroid_ref)*x2)*temp;
	return;
}
//void discontinuous_point::calcualte_d1(vec2D& centroid_ref,double area_below)
//{
//	double temp=_length*_length*_length/(12*area_below);
//	vec2D x1(sin(_theta)*temp,-cos(_theta)*temp);
//	//vec2D x1(_s*temp,-_c*temp);
//	_d1=(_centroid_below-centroid_ref)*x1*2;
//}
double area(vec2D& p1,vec2D& p2,vec2D& p3)
{
	double nx=(p2.x-p1.x)*(p3.y-p1.y);
	double ny=(p3.x-p1.x)*(p2.y-p1.y);
	return 0.5*abs((p2.x-p1.x)*(p3.y-p1.y)-(p3.x-p1.x)*(p2.y-p1.y));
}
double root_of_cubic_interpolation(iteration_point& iteration_min,iteration_point& iteration_max)
{
	double d_min=iteration_min.theta,f_min=iteration_min.d1,df1_min=iteration_min.d2;
	double d_max=iteration_max.theta,f_max=iteration_max.d1,df1_max=iteration_max.d2;
	double a,b,c,d;
	double d1=d_min-d_max,d2=d1*d1,d3=d2*d1;

	//The cubic function can be constructed and calculate the root
	a=(df1_min+df1_max)-(2*f_min-2*f_max)/d1;
	b=-(d_min*d_min*(df1_min+2*df1_max)-d_max*d_max*(2*df1_min+df1_max)-d_min*(3*f_min-3*f_max)+d_max*(3*f_max-3*f_min+d_min*(df1_min-df1_max)))/d1;
	c=df1_max*d2+(2*d_max*(df1_min+2*df1_max))*d1-(6*d_max*d_max*(f_min-f_max))/d1 + (3*d_max*(2*f_max-2*f_min+d_max*df1_min+d_max*df1_max));
	d=(d_max*d_max*f_min*(3*d_min-d_max))/d1 - (d_min*d_max*d_max*df1_min) - (d_min*d_min*d_max*df1_max) + (d_min*d_min*f_max*(d_min-3*d_max))/d1;

	double iterator=0.5*(d_min+d_max);
	double iterator2,iterator3;
	double f,df;
	iterator2=iterator*iterator;iterator3=iterator2*iterator;
	f=a*iterator3+b*iterator2+c*iterator+d;
	df=3*a*iterator2+2*b*iterator+c;
	while (true)
	{
		iterator=iterator-f/df;
		if (iterator<d_min || iterator>d_max)
		{
			iterator=0.5*(d_max+d_min);
		}
		iterator2=iterator*iterator;iterator3=iterator2*iterator;
		f=a*iterator3+b*iterator2+c*iterator+d;
		df=3*a*iterator2+2*b*iterator+c;
		if (abs(f)<1e-10)
		{
			return iterator;
		}
		else if (d_max-d_min<1e-11)
		{
			return (d_max+d_min)*0.5;
		}
		if (f<0)
		{
			d_min=iterator;
		}
		else if (f>0)
		{
			d_max=iterator;
		}
	}
	return iterator;
}
double root_of_fifth_interpolation(iteration_point& iteration_min,iteration_point& iteration_max)
{
	double d_min=iteration_min.theta,f_min=iteration_min.d1,df1_min=iteration_min.d2,df2_min=iteration_min.d3;
	double d_max=iteration_max.theta,f_max=iteration_max.d1,df1_max=iteration_max.d2,df2_max=iteration_max.d3;
	double a,b,c,d,e,f;
	double d1=d_min-d_max,d2=d1*d1,d3=d2*d1,d4=d3*d1,d5=d4*d1;
	double d_min2=d_min*d_min,d_min3=d_min2*d_min,d_min4=d_min3*d_min,d_min5=d_min4*d_min;
	double d_max2=d_max*d_max,d_max3=d_max2*d_max,d_max4=d_max3*d_max,d_max5=d_max4*d_max;

	a=(6*f_min)/d5 - (3*df1_max)/d4 - (3*df1_min)/d4 - (6*f_max)/d5 + df2_min/(2*d_min3 - 6*d_min2*d_max + 6*d_min*d_max2 - 2*d_max3) - df2_max/(2*d_min3 - 6*d_min2*d_max + 6*d_min*d_max2 - 2*d_max3);
	b=(df2_max*(3*d_min + 2*d_max))/(2*d_min3 - 6*d_min2*d_max + 6*d_min*d_max2 - 2*d_max3) - (df2_min*(2*d_min + 3*d_max))/(2*d_min3 - 6*d_min2*d_max + 6*d_min*d_max2 - 2*d_max3) + (df1_min*(7*d_min + 8*d_max))/d4 + (df1_max*(8*d_min + 7*d_max))/d4 - (f_min*(15*d_min + 15*d_max))/d5 + (f_max*(15*d_min + 15*d_max))/d5;
	//k1=df2_min-20*a*d_min3-12*b*d_min2,k2=df2_max-20*a*d_max3-12*b*d_max2;
	c=(df2_min-20*a*d_min3-12*b*d_min2)/(6*d1) - (df2_max-20*a*d_max3-12*b*d_max2)/(6*d1);
	d=(df2_min-20*a*d_min3-12*b*d_min2-6*c*d_min)/2;
	//d=(d_min*(df2_max-20*a*d_max3-12*b*d_max2))/(2*d1) - (d_max*(df2_min-20*a*d_min3-12*b*d_min2))/(2*d1);
	e=df1_min-5*a*d_min4-4*b*d_min3-3*c*d_min2-2*d*d_min;
	f=f_min-a*d_min5-b*d_min4-c*d_min3-d*d_min2-e*d_min;
	//double iterator=0.5*(d_min+d_max);
	double iterator=d_max+f_max*d1/(f_max-f_min);
	double iterator2,iterator3,iterator4,iterator5;
	double g,dg;
	iterator2=iterator*iterator;iterator3=iterator2*iterator;iterator4=iterator3*iterator;iterator5=iterator4*iterator;
	g=a*iterator5+b*iterator4+c*iterator3+d*iterator2+e*iterator+f;
	dg=5*a*iterator4+4*b*iterator3+3*c*iterator2+2*d*iterator+e;
	while (true)
	{
		iterator=iterator-g/dg;
		if (iterator<d_min || iterator>d_max)
		{
			iterator=0.5*(d_max+d_min);
		}
	iterator2=iterator*iterator;iterator3=iterator2*iterator;iterator4=iterator3*iterator;iterator5=iterator4*iterator;
	g=a*iterator5+b*iterator4+c*iterator3+d*iterator2+e*iterator+f;
	dg=5*a*iterator4+4*b*iterator3+3*c*iterator2+2*d*iterator+e;
		if (abs(g)<1e-10)
		{
			return iterator;
		}
		else if (d_max-d_min<1e-10)
		{
			return (d_max+d_min)*0.5;
		}
		if (g<0)
		{
			d_min=iterator;
		}
		else if (g>0)
		{
			d_max=iterator;
		}
	}
	return iterator;
}
double calculate_inclination(vec2D& dp)
{
	double angle=acos(dp.x/(sqrt(dp*dp)));
	if (dp.y<0)
	{
		angle=2*pi-angle;
	}
	return angle;
}
double minval(double a,double b)
{
	if (a>=b)
	{
		return b;
	}
	else
	{
		return a;
	}
}
double maxval(double a,double b)
{
	if (a<=b)
	{
		return b;
	}
	else
	{
		return a;
	}
}
bool is_intersect(vec2D x_range_old,vec2D y_range_old,vec2D x_range_new,vec2D y_range_new)
{
	if (x_range_old.x>x_range_new.y)
	{
		return false;
	}
	if (x_range_old.y<x_range_new.x)
	{
		return false;
	}
	if (y_range_old.x>y_range_new.y)
	{
		return false;
	}
	if (y_range_old.y<y_range_new.x)
	{
		return false;
	}
	return true;
}
vec2D calculate_moment(double xl,double yl,double xr,double yr)
{
	vec2D xs,xt;
	xs.value(0.5*xl+0.5*xr,yr*0.5);
	xt.value((2*xl+xr)/3.0,(yl+2*yr)/3.0);
	return xs*(xr-xl)*yr+xt*(xr-xl)*(yl-yr)*0.5;
}