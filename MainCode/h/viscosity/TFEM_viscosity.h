#ifndef TFEM_VISCOSITY
#define TFEM_VISCOSITY
//Define the class for FEM artificial viscosity
#include "Tviscosity_base.h"
#include "data_structure.h"
class TFEM_viscosity :public Tviscosity_base
{
public:
	TFEM_viscosity(string name,double k1);
	virtual void calculate_viscosity_corner_force(Tcell_fluid_base* this_cell,vec3D (&assistant)[8][3],double dt);


protected:
	double _k1;

	void initialize_variables(double (&phy)[8][8],vec3D (&dphy)[8][8]);
	//Calculate the value and derivatives of shape functions on Gausspoints
	void get_key_variable(Tcell_fluid_base* this_cell,vec3D (dphy)[8][8],double (&det)[8],vec3D (&key)[8][8]);
	//Calculate the key vector and determinant of Jacobi
	double get_qs(Tcell_fluid_base* this_cell);
	//Calculate the magnitude of the artificial viscosity

};
#endif