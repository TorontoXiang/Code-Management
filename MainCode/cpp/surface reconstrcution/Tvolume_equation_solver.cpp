#include "Tvolume_equation_solver.h"
#include <iostream>
#include <algorithm>
using namespace std;
void Tvolume_equation_solver::set_target_polyhedron(Tpolyhedron* poly)
{
	_poly=poly;
	_poly->calculate_piece_topology();
	_poly->calculate_polyhedron_volume();
}
double Tvolume_equation_solver::calculate_plane_constant_advanced(vec3D &n,double fraction,double predict)
{
	//return calculate_plane_constant(n,fraction,predict);
	_num_solving=_num_solving+1;
	double eps=5e-22;
	//The iterative variables
	double d_min,d_max,f_min,f_max,df_min,df_max;
	vec3D info;    //Information at iteration point
	double d_iterator,f_iterator,df_iterator;
	//[d_min,d_max] is the region where the solution locates,thus f_min<0 and f_max>0
	_poly->calculate_altitude_range(n,d_min,d_max);
	f_min=-fraction,f_max=1-fraction;df_min=df_max=0;
	//Calculate the first iteration point
	vec3D min_info,max_info;
	if (predict<d_min||predict>d_max)
	{
		//If the predict altitude is excess,using ordinary guess
		min_info.x=d_min;min_info.y=-fraction;min_info.z=0;
		max_info.x=d_max;max_info.y=1-fraction;max_info.z=0;
		d_iterator=root_of_approximate_volume_function(min_info,max_info,fraction);
		//double a=d_min,b=d_max; 
		//double fa=-fraction,fb=1-fraction,f=fa;
		//d_iterator=a-fa*(b-a)/(fb-fa);

	}
	else
	{
		//If the predict altitude is reasonable, set it as initial guess
		d_iterator=predict;
	}
	info=calculate_fraction_below(n,d_iterator);
	f_iterator=info.x-fraction;df_iterator=info.y;
	if (abs(f_iterator)<eps)
	{
		return d_iterator;
	}
	else if (f_iterator>0)
	{
		d_max=d_iterator;f_max=f_iterator;df_max=df_iterator;
	}
	else if (f_iterator<0)
	{
		d_min=d_iterator;f_min=f_iterator;df_min=df_iterator;
	}
	int num_iteration=0;
	while (abs(f_iterator)>eps && abs(f_iterator/df_iterator)>eps && num_iteration<20)
	{
		d_iterator=newton(d_iterator,f_iterator,df_iterator);
		if (d_iterator<d_max && d_iterator>d_min)
		{
			//The new iteration point is in [d_min,d_max]
			info=calculate_fraction_below(n,d_iterator);
			f_iterator=info.x-fraction;df_iterator=info.y;
			num_iteration=num_iteration+1;
			if (abs(f_iterator)<eps)
			{
				return d_iterator;
			}
			else if (f_iterator>0)
			{
				d_max=d_iterator;f_max=f_iterator;df_max=df_iterator;
			}
			else if (f_iterator<0)
			{
				d_min=d_iterator;f_min=f_iterator;df_min=df_iterator;
			}
		}
		else
		{
			min_info.x=d_min;min_info.y=f_min;min_info.z=df_min;
			max_info.x=d_max;max_info.y=f_max;max_info.z=df_max;
			d_iterator=root_of_approximate_volume_function(min_info,max_info,fraction);
			//d_iterator=0.5*(d_min+d_max);
			info=calculate_fraction_below(n,d_iterator);
			f_iterator=info.x-fraction;df_iterator=info.y;
			num_iteration=num_iteration+1;
			if (abs(f_iterator)<eps)
			{
				return d_iterator;
			}
			else if (f_iterator>0)
			{
				d_max=d_iterator;f_max=f_iterator;df_max=df_iterator;
			}
			else if (f_iterator<0)
			{
				d_min=d_iterator;f_min=f_iterator;df_min=df_iterator;
			}
		}
	}
	return d_iterator;
}
double Tvolume_equation_solver::calculate_plane_constant(vec3D &n,double fraction,double predict)
{
	double d_min,d_max;          //The altitude range
	double eps=1e-16;
	int n_secant=0,n_newton=0;
	_poly->calculate_altitude_range(n,d_min,d_max);
	double a=d_min,b=d_max;               //The endpoints of the iterative section
	double fa=-fraction,fb=1-fraction;    //The value at the endpoints
	double f=fa,df;                       //The value and its derivative at the iterative point
	vec3D info;                           //Store f and df
	fa=-fraction;fb=1-fraction;f=fa;
	double iterator;
	//Determine the initial guess
	//Newton method as first
	vec3D min_info,max_info;
	if (predict<d_min||predict>d_max)
	{
		//If the predict altitude is excess,using ordinary guess
		iterator=a-fa*(b-a)/(fb-fa);
		min_info.x=d_min;min_info.y=-fraction;min_info.z=0;
		max_info.x=d_max;max_info.y=1-fraction;max_info.z=0;

		//iterator=root_of_approximate_volume_function(min_info,max_info,fraction);
	}
	else
	{
		//If the predict altitude is reasonable, set it as initial guess
		iterator=predict;
	}
	info=calculate_fraction_below(n,iterator);
	f=info.x-fraction;df=info.y;
	//Using Newton iteration first
	while (abs(f)>eps)
	{
		iterator=newton(iterator,f,df);
		info=calculate_fraction_below(n,iterator);
		f=info.x-fraction;df=info.y;
		n_newton=n_newton+1;
		if (abs(f)<=eps)
		{
			return newton(iterator,f,df);
		}
		if (iterator<d_min||iterator>d_max)
		{
			_num_Newton_fails=_num_Newton_fails+1;
			break;
		}
	}
	//If the Newton iteration failed, using secant method subsequently
	while (abs(f)>eps)
	{
		iterator=a-fa*(b-a)/(fb-fa);
		f=calculate_fraction_below(n,iterator).x-fraction;
		if (n_secant>=10)
		{
			break;
		}
		secant(a,b,fa,fb,f);
		n_secant=n_secant+1;
	}
	//If f is convergent, return the constant
	if (abs(f)<=eps)
	{
		return a-fa*(b-a)/(fb-fa);
	}
	//If the scant method is also failed, bisection method will be used finally
	a=d_min,b=d_max;
	while (abs(a-b)>=eps)
	{
		//The terminal condition is the length of a,b and f might not less than eps due to numerical error in
		//calculating the fraction, so Newton method will be performed if f>eps
		iterator=(a+b)/2.0;
		info=calculate_fraction_below(n,iterator);
		f=info.x-fraction;
		if (abs(f)<=eps)
		{
			return (a+b)*0.5;
		}
		bisect(a,b,f);
	}
	//If the f resulting from bisection method is not less than eps, Newton method will be performed for modification 
	double iterator0=iterator;f=info.x-fraction;df=info.y;
	n_newton=0;
	while (abs(f)>eps)
	{
		iterator=newton(iterator,f,df);
		info=calculate_fraction_below(n,iterator);
		f=info.x-fraction;df=info.y;
		n_newton=n_newton+1;
		if (n_newton==10||iterator<d_min||iterator>d_max)
		{
			return iterator0;
		}
	}
	if (abs(f)<eps)
	{	
		return newton(iterator,f,df);
	}
	return iterator;
}
vec3D Tvolume_equation_solver::calculate_fraction_below(vec3D &n,double d)
{
	return _poly->calculate_fraction_below(n,d);
}