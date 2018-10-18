#ifndef TCELL_BRICK
#define TCELL_BRICK
#include "Tcell_base.h"
#include "Tnode.h"
#include "data_structure.h"
#include <vector>
#include "Tviscosity_base.h"
//Define cell shell
class Tcell_brick : public Tcell_base
{
public:
	Tcell_brick() {};
	Tcell_brick(int cell_id, int nGauss, Tnode* (&node_ptr)[8], TMAT_base** mat);


	void calculate_corner_force(Tviscosity_base* Arti_vis, double dt);
	//Calculate the corner force
	void assemble_corner_force();
	//Assemble the corner force into nodal force and clear the corner force
	double calculate_time_step();
	//Calculate the cell time step
	void assemble_corner_inertance();
	//Assemble the corner inertance into node
	void update_strain_energy(double dt);
protected:
	Tnode* _node_ptr[8];
	double _clength; //The characteristic length of the shell
	double _volume;  //The volume of the cell
	double _soundspeed_max;    //The maximal sound speed at Gausspoints
	double _strain_energy;            //The strain energy in this cell

	void calculate_hourglass_viscosity_force(vec3D(&HG_force)[8]);
	//Calculate the hourglass viscosity term
	double calculate_characteristic_length();
	//Calculate the characteristic length
	void calculate_Jacobi(double (&Jacobi)[3][3], vec3D (&GN)[8]);
	//Calculate the Jacobi matrix at the center;
	void calculate_DN_matrix(vec3D(&DN)[8], double &detJ);
	//Calculate the derivatives of shape function with respected to x,y,z at the center
	double calculate_cell_volume();
	//Calculate the cell volume from Jacobian matrix

	//Friend Class List
	friend class Tbody_brick;

};
#endif
