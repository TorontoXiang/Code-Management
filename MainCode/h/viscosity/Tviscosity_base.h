#ifndef TVISCOSITY_BASE
#define TVISCOSITY_BASE
//Define the base class for artificial viscosity
#include "Tcell_fluid_base.h"
#include "data_structure.h"
#include <string>
using namespace std;
class Tcell_fluid_base;
class Tviscosity_base
{
public:
	Tviscosity_base(string name):_av_name(name){};

	//Viscosity for solid
	virtual double calculate_q(double density, double vrate, double soundspeed, double clength){return 0;};
	//Calculate the viscosity coefficient q in solid material
	virtual	double modify_soundspeed(double vrate,double soundspeed,double clength){return 0;};
	//Modify the soundspeed by artificial viscosity

	//Viscosity for fluid
	virtual void calculate_viscosity_corner_force(Tcell_fluid_base* this_cell,vec3D (&assistant)[8][3],double dt){};
	//Calculate the viscosity force for fluid material

	//Access function
	string G_name(){return _av_name;};
protected:
	string _av_name;
};
class Tsolid_viscosity : public Tviscosity_base
{
public:
	Tsolid_viscosity(string name,double c1,double c2);
	//Calculate the viscosity coefficient
	virtual double calculate_q(double density, double vrate, double soundspeed, double clength);
	//Modify the soundspeed by artificial viscosity
	virtual double modify_soundspeed(double vrate,double soundspeed, double clength);
private:
	double _k1,_k2;  
};
#endif