#ifndef TCELL_BRICK
#define TCELL_BRICK
#include "Tnode.h"
#include "data_structure.h"
#include <vector>
#include "TMAT_base.h"
//Define cell shell
class Tcell_brick
{
public:
	Tcell_brick() {};
	Tcell_brick(int cell_id, Tnode* (&node_ptr)[8], TMAT_base* mat);


	void calculate_corner_force(double dt);
	//Calculate the corner force
	void assemble_corner_force();
	//Assemble the corner force into nodal force and clear the corner force
	double calculate_time_step();
	//Calculate the cell time step
	void assemble_corner_inertance();
	//Assemble the corner inertance into node
	void update_strain_energy(double dt);
protected:
	int _id;                    //Cell id
	TMAT_base* _mat_ptr;        //Cell material
	Tnode* _node_ptr[8];        //Cell node pointer
	vec3D _corner_force[8];     //The nodal force from this cell
	double _clength;            //The characteristic length of the shell
	double _volume;             //The volume of the cell
	double _soundspeed;         //The maximal sound speed at Gausspoints
	double _strain_energy;      //The strain energy in this cell

	void calculate_hourglass_viscosity_force(vec3D(&HG_force)[8]);
	//Calculate the hourglass viscosity term
	double calculate_characteristic_length();
	//Calculate the characteristic length
	void calculate_Jacobi(double(&Jacobi)[3][3], vec3D(&GN)[8]);
	//Calculate the Jacobi matrix at the center;
	void calculate_DN_matrix(vec3D(&DN)[8], double &detJ);
	//Calculate the derivatives of shape function with respected to x,y,z at the center
	double calculate_cell_volume();
	//Calculate the cell volume from Jacobian matrix

	//Access Functions
	double G_strain_energy() { return _strain_energy; };

	//Friend Class List
	friend class Tbody;

};
#endif
