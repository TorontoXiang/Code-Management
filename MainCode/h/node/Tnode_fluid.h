#pragma once
#ifndef TNODE_FLUID
#define TNODE_FLUID
//Define the node in fluid simulation
#include "Tnode.h"
#include "Tcell_fluid_base.h"
#include <vector>
class Tcell_fluid_base;
using namespace std;

class Tnode_fluid : public Tnode
{
public:
		
	Tnode_fluid(){};
	Tnode_fluid(int id,vec3D coordinate);
	Tnode_fluid(int id,double x,double y,double z);

	void clear_nodal_variable();
	//Clear the nodal variable on old grid after remapping
	void predict_position(double dt);
	//Update the position in predictor step
	void add_connected_cell(Tcell_fluid_base* connected_cell);
	//Add the connected cell to node
	void calculate_nodal_variables();
	//Calculate _density,_soundspeed and _pressure by the average of _surrounding
	void reset_velocity_after_remapping(vec3D velocity_origin);
	//Reset the velocity by _moment and _mass
	//velocity_origin: the origin velocity before of the old grid for co-node
	void set_velocity_from_particle_after_remapping();
	//Set the velocity from the particle after remapping
	void perturbance_node_for_implosion(vec3D center,int dimension,double ratio,double n);
	//perturbance the node coordinate for implosion problem
	void move_node_for_implosion(vec3D center,int dimension,int N,double d0,double r_in,double r_out,double alpha);
	//move the node coordiante for implosion problem
	void set_velocity_bc(double t,int dimension);
private:
	vector<Tcell_fluid_base*> _surroundings;
	vec3D _position_temp;     //Position at predictor step
	vec3D _position_beginning;//Position at the beginning of grid smoothing
	double _density;          //Nodal density
	double _soundspeed;       //Nodal sound speed
	double _pressure;         //Nodal pressure


	friend class Tedge_viscosity;
	friend class TFEM_viscosity;
	template<class T1,class T2>
	friend class Tinteraction_conode;
	friend class Tgrid_smooth_continum;
	friend class Tcell_fluid_base;
	friend class Tcell_pure_fluid;
	friend class Tcell_mixed_fluid;
	friend class Tbody_Lagrangian;
	friend class Tbody_ALE;
	friend class Tbody_explosion;
	friend class Tbody_MPM;
	friend class Tbody_ALEMPM;
};
#endif