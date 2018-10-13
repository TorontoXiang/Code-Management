#ifndef TFEM_VISCOSITY_TERM
#define TFEM_VISCOSITY_TERM
//Define the class for FEM artificial viscosity
#include "data_structure.h"
#include "Tcell_fluid_base.h"
class TFEM_fluid_viscosity
{
public:
	void calculate_fluid_viscosity_corner_force(Tcell_fluid_base* this_cell);
	//Calculate the fluid viscosity force
protected:
	double _mu;

	void initialize_variables(double (&phy)[8][8],vec3D (&dphy)[8][8]);
	//Calculate the value and derivatives of shape functions on Gausspoints
	void get_key_variable(Tcell_fluid_base* this_cell,vec3D (dphy)[8][8],double (&det)[8],vec3D (&key)[8][8]);
	//Calculate the key vector and determinant of Jacobi
};
#endif