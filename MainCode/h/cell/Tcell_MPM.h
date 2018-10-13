#ifndef TCELL_MPM
#define TCELL_MPM
//Define the class of MPM cell
#include "Tnode.h"
#include "TMPM_particle.h"
class Tcell_MPM
{
public:
	static vec3D _length;         //The length of a regular background cell
public:
	Tcell_MPM():_max_spd(-1),_clength(-1),_is_empty(true),_cell_id(-1){};
	void calculate_node_ptr(Tnode* node_ptr[8]);
	//Calculate the node_ptr of a cell
	vec3D calculate_standard_coordinate(vec3D& particle_coordinate);
	//Calculate the standard coordinate of a particle in this cell
	void distribute_particle_to_node(TMPM_particle* particle_ptr);
	//Distribute the particle mass and moment to background node
	void calculate_cell_length();
	//Calculate the character length of this cell
	double calculate_cell_time_step(double& clength,double& soundspeed,double& nodespeed);
	//Calculate the time step of this cell
	void initialize();
	//Initialize the background cell
	void distribute_particle_force(TMPM_particle* particle_ptr);
	//Distribute the particle force to the background node
	void assemble_corner_force();
	//Assemble the corner force to the node
	void update_particle_state(TMPM_particle* particle_ptr,double dt);
	//Update the particle's stress,density,volume and soundspeed in this cell
	void update_particle_DOF(TMPM_particle* particle_ptr,double dt);
	//Update the particle DOF in this cell
	
private:
	int _cell_id;
	Tnode* _node_ptr[8];
	double _max_spd;           //The maximal soundspeed in this cell
	double _clength;           //The character length of this cell
	vec3D _corner_force[8];    //The 8 corner force of a cell
	bool _is_empty;            //Whether any particle is in this cell

	void calculate_DN_matrix(vec3D (&DN)[8],double &detJ,double xi,double eta,double zeta);
	//Calculate the dervative of shape functions with respect to x,y,z at particle
	void calcualte_Jacobi(double (&Jacobi)[3][3],vec3D (&GN)[8],double xi,double eta,double zeta);
	//Calculate the GN matrix and Jacobi matrix at particle

	friend class Tbody_MPM;
	friend class TMPM_background;


};
#endif