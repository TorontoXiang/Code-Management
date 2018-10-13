#include "Tviscosity_base.h"
#include <cmath>
Tsolid_viscosity::Tsolid_viscosity(string name,double c1,double c2):Tviscosity_base(name),_k1(c1),_k2(c2){}
double Tsolid_viscosity::calculate_q(double density, double vrate, double soundspeed, double clength)
{
	if (vrate<0)
	{
		return density*clength*vrate*(_k1*clength*vrate-_k2*soundspeed);
	}
	else
	{
		return 0;
	}
}
double Tsolid_viscosity::modify_soundspeed(double vrate, double soundspeed, double clength)
{
	double temp=_k2*soundspeed-_k1*clength*vrate;
	if (vrate<0)
	{
		return temp+sqrt(temp*temp+soundspeed*soundspeed);
	}
	else
	{
		return soundspeed;
	}
}