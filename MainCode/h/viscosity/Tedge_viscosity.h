#ifndef TEDGE_VISCOSITY
#define TEDGE_VISCOSITY
//Define edge artificial viscosity in fluid element
#include "Tviscosity_base.h"
class Tedge_viscosity : public Tviscosity_base
{
public:
	Tedge_viscosity(string name,double c1,double c2):Tviscosity_base(name),_k1(c1),_k2(c2){};

	virtual void calculate_viscosity_corner_force(Tcell_fluid_base* this_cell,vec3D (&assistant)[8][3],double dt);
private:
	double _k1,_k2;

};
#endif