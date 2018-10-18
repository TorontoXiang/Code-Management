#ifndef TNODE
#define TNODE
#include "data_structure.h"
//Define the node base class with position,velocity and so on
class Tnode
{
public:
	Tnode();
	Tnode(int id,double x,double y,double z);

	void calculate_acceleration();
	//Calculate the acceleration
	void update_vel(double dt);
	//Update the velocity in dt
	void trial_velocity(double dt);
	//Calculate the trial velocity
	void back_to_origin_pos(double dt);
	//Return to the origin position after volume check
	void update_pos(double dt);
	//Update the position in dt
	void trial_pos(double dt);
	//Trial the position before the update
	void clear_force();
	//Clear the nodal force after calculating acceleration
	void apply_external_force(int direction,double magnitude);
	//Apply the external force on node
	void apply_gravity(vec3D &gravity);
	//Apply the gravity
	void apply_gravity_MPM(vec3D &grtavity);
	//Apply the gravity for particles
	void reset_position(){_position=_position0;};
	//Set the position of a node as its origin position

	//Access function
	vec3D calculate_displacement(){return _position-_position0;};
	//Return the displacement of a node
	vec3D calculate_displacement_amplify(double ratio){return _position0+(_position-_position0)*ratio;};
	//Retrun the amplified displacement of a node in plot
	void access_nodal_initial_inertance(double &mass,double &mass_rotation){mass=_mass0,mass_rotation=0;};
	//Access the nodal inertance without conode addition
	void conode_mass_addition(double conode_mass,double conode_mass_rotation){_mass=_mass+conode_mass;};
	//Add the co-node mass to this node
	void reset_nodal_inertance(){_mass=_mass0;};
	//Reset nodal inertance by its inital value before the modification after remapping
	void access_nodal_initial_force(vec3D &force,vec3D &moment){force=_force0;moment.value(0,0,0);};
	//Access the nodal force without conode addition
	void conode_force_addition(vec3D conode_force,vec3D conode_moment){_force=_force+conode_force;};
	//Add the co-node force to this node
	void access_nodal_bc(int (&bc_type_pos)[3],double (&bc_vel)[3]);
	//Get the nodal boundary condition

protected:
	int _id;               //node id
	vec3D _position0;      //The original position 
	vec3D _position;       //The current position
	vec3D _position_trial; //The trial position before update 
	vec3D _velocity;       //velpcity
	vec3D _velocity_temp;  //temporary velocity
	vec3D _velocity_trial; //The trial velocity for volume check
	double _mass;          //mass with co-node addition (used in computation)
	double _mass0;         //mass without co-node addition
	vec3D _moment;         //The nodal moment
	vec3D _force;          //nodal force
	vec3D _force0;         //nodal force without co-node addition
	vec3D _acceleration;   //acceleration
	//Boundary condition
	int _bc_type_position[3];   //The boundary for position: 0-no;1-fixed;2-velocity boundary
	double _bc_velocity[3];        //The value of velocity boundary
	vec3D _external_load;       //The external load

	int _is_conode;             //Whether the node is a conode with other body
	int _flag;

	void accumulate_inertance(double mass);
	//Accumulate nodal inertance
	void accumulate_moment(vec3D moment);
	//Accumulate nodal moment
	void initialize();
	//Initialize the node variables

	friend class Tcell_shell;
	friend class Tbody_shell;
	friend class Tcell_intersect;
	friend class Tcell_8H;
	friend class Tgrid_smooth_continum;
	friend class Tgrid_smooth_OST;
	friend class Tgrid_smooth_CST;
	friend class Tgrid_smooth_PST;
	friend class TMPM_background;
	friend class Tcell_MPM;
	friend class Tbody_MPM;
	friend class Tbody_ALEMPM;
	friend class Tbody_explosion;
	friend class Tcell_pure_fluid_MPM;
	friend class Tcell_mixed_fluid_MPM;
	friend class Tcell_brick;
	friend class Tbody_brick;
};
//--------------------------------------------------------------------------------
//Define node with rotation freedom
//--------------------------------------------------------------------------------
class Tnode_rotation : public Tnode
{
public:
	Tnode_rotation(){};
	Tnode_rotation(int id,double x,double y,double z);
	void calculate_acceleration();
	//Calculate the acceleration
	void update_vel(double dt);
	//Update velocity and Oemga in dt
	void clear_force();
	//Clear the acceleration after updating velocity
	void access_nodal_initial_inertance(double &mass,double &mass_rotation){mass=_mass0;mass_rotation=_mass_rotation0;};
	//Access the nodal inertance without conode addition
	void conode_mass_addition(double conode_mass,double conode_mass_rotation){_mass=_mass+conode_mass;_mass_rotation=_mass_rotation+conode_mass_rotation;};
	//Add the co-node mass to this node
	void reset_nodal_inertance(){_mass=_mass0;_mass_rotation=_mass_rotation0;};
	//Reset nodal inertance by its inital value before the modification after remapping
	void access_nodal_initial_force(vec3D &force,vec3D &moment){force=_force0;moment=_moment0;};
	//Access the nodal force without conode addition
	void conode_force_addition(vec3D conode_force,vec3D conode_moment){_force=_force+conode_force;_moment=_moment+conode_moment;};
	//Add the co-node force to this node
protected:
	vec3D _Omega;          //rotation rate
	vec3D _Omega_temp;     //temporary rotation rate
	vec3D _moment;         //nodal moment
	vec3D _moment0;        //nodal moment without co-node addition
	vec3D _epsilon;        //rotation acceleration
	double _mass_rotation; //mass rotation with co-node addition
	double _mass_rotation0;//mass rotation with out co-node addition
	//Boundary condition
	int _bc_type_rotation[3];    //The boundary for rotation: 0-no;1-fixed;2-Omega 
	double _bc_Omega[3];         //The value of Omega condition

	void accumulate_inertance(double mass,double mass_rot);
	//Accumulate nodal inertance

	friend class Tcell_shell;
	friend class Tbody_shell;
	template<class T1,class T2>
	friend class Tinteraction_conode;
};
#endif