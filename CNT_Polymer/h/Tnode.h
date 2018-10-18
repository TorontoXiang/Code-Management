#ifndef TNODE
#define TNODE
#include "data_structure.h"
//Define the node base class with position,velocity and so on
class Tnode
{
public:
	Tnode();
	Tnode(int id, double x, double y, double z);

	void calculate_acceleration();
	//Calculate the acceleration
	void update_vel(double dt);
	//Update the velocity at the N+1/2 step
	void update_pos(double dt);
	//Update the position at N step
	void clear_force();
	//Clear the nodal force after calculating acceleration
	void apply_external_force(int direction, double magnitude);
	//Apply the external force on node
	void apply_gravity(vec3D &gravity);
	//Apply the gravity

	//Access function
	vec3D calculate_displacement() { return _position-_position0; };
	//Return the displacement of a node
	vec3D calculate_displacement_amplify(double ratio) { return _position0 + (_position - _position0)*ratio; };
	//Retrun the amplified displacement of a node in plot
	void access_nodal_bc(int(&bc_type_pos)[3], double(&bc_vel)[3]);
	//Get the nodal boundary condition

protected:
	int _id;               //node id
	vec3D _position0;      //The original position 
	vec3D _position;       //The current position
	vec3D _velocity;       //velocity
	double _mass;          //mass with co-node addition (used in computation)
	vec3D _force;          //nodal force
	vec3D _acceleration;   //acceleration
	//Boundary condition
	int _bc_type_position[3];   //The boundary for position: 0-no;1-fixed;2-velocity boundary
	double _bc_velocity[3];     //The value of velocity boundary
	vec3D _external_load;       //The external load

	int _flag;

	void accumulate_inertance(double mass);
	//Accumulate nodal inertance
	void initialize();
	//Initialize the node variables

	friend class Tcell_brick;
	friend class Tbody;
};
#endif