#ifndef TCELL_PURE_FLUID_MPM
#define TCELL_PURE_FLUID_MPM
//Define the cell for ALEMPM of pure cell
#include "Tcell_pure_fluid.h"

class Tcell_pure_fluid_MPM : public Tcell_pure_fluid
{
public:
	Tcell_pure_fluid_MPM(int cell_id,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat):Tcell_pure_fluid(cell_id,node_ptr,mat){};
	virtual double calculate_time_step_particle();
	//Calculate the time step from particle
	virtual void distribute_particle_force(TMPM_particle* particle_ptr);
	//Distribute the particle force to the background node
	virtual void distribute_particle_variables(TMPM_particle* particle_ptr,vec3D interval);
	//Distribute the particle mass and moment to the cell node
	virtual void update_particle_DOF(TMPM_particle* particle_ptr,double dt);
	//Update the DOF of particles
	virtual void update_particle_stress(TMPM_particle* particle_ptr,double dt);
	//Update the state of particles
	virtual void assemble_corner_force();
	//Assemble the corner force to nodes

	virtual double G_max_spd_particle(){return _maxspd_particle;};
	virtual void S_max_spd_particle(double max_spd){_maxspd_particle=max_spd;};
protected:
	double _maxspd_particle;    //The maximal soundspeed of particles in the cell
	bool _is_empty;             //Whether a particle is in the cell
	vec3D _corner_stress_force_particle[8];    //The corner force from particles


	virtual vec3D calculate_standard_coordinate(vec3D& particle_coordinate,vec3D interval);
	//Calculate the standard coordinate of a particle in this cell
	virtual void reset_state_after_remapping();
	//Reset the variables after remapping phase
};
#endif