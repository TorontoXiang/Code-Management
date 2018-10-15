#include "Tnode.h"
#include <iostream>
using namespace std;
//-----------------------------------------------------------------------------------------
//Functions for Tnode_base
//-----------------------------------------------------------------------------------------
Tnode::Tnode()
{
	_id=-1;
	_mass=_mass0=0;
	_is_conode=0;
	for (int i = 0; i < 3; i++)
	{
		_bc_type_position[i]=0;
		_bc_velocity[i]=0;
	}
}
Tnode::Tnode(int id,double x,double y,double z)
{
	_id=id;
	_position0.x=x;_position0.y=y;_position0.z=z;
	_position=_position0;
	_velocity.value(0,0,0);
	_velocity_temp.value(0,0,0);
	_force.value(0,0,0);_force0.value(0,0,0);
	_acceleration.value(0,0,0);
	_mass=_mass0=0;
	_is_conode=0;

	for (int i = 0; i < 3; i++)
	{
		_bc_type_position[i]=0;
		_bc_velocity[i]=0;
	}
}
void Tnode::calculate_acceleration()
{
	if (_mass>0)
	{
		_acceleration=_force/_mass;
	}
	else
	{
		_acceleration.value(0,0,0);
	}
	//Apply the boundary condition
	for (int i = 0; i < 3; i++)
	{
		if (_bc_type_position[i]!=0)
		{
			_acceleration.value(i+1,0);
		}
	}
	//Clear the force after calculating acceleration
	clear_force();
	return;
}
void Tnode::trial_velocity(double dt)
{
	_velocity_temp=_velocity;
	_velocity_trial=_velocity+_acceleration*dt;
	_velocity_temp=(_velocity+_velocity_trial)*0.5;
	for (int i = 0; i < 3; i++)
	{
		if (_bc_type_position[i]==1)
		{
			_velocity_temp.value(i+1,0);
			_velocity.value(i+1,0);
		}
	}
	return;
}
void Tnode::update_vel(double dt)
{
	_velocity_temp=_velocity;
	_velocity=_velocity+_acceleration*dt;
	_velocity_temp=(_velocity+_velocity_temp)*0.5;
	for (int i = 0; i < 3; i++)
	{
		if (_bc_type_position[i]==1)
		{
			_velocity_temp.value(i+1,0);
			_velocity.value(i+1,0);
		}
		else if (_bc_type_position[i]==2)
		{
			_velocity_temp.value(i+1,_bc_velocity[i]);
			_velocity.value(i+1,_bc_velocity[i]);
		}
	}
	return;
}
void Tnode::update_pos(double dt)
{
	_position=_position+_velocity_temp*dt;
	return;
}
void Tnode::trial_pos(double dt)
{
	_position_trial=_position+_velocity_temp*dt;
	return;
}
void Tnode::back_to_origin_pos(double dt)
{
	_position=_position-_velocity_temp*dt;
	return;
}
void Tnode::clear_force()
{
	_force.value(0,0,0);_force0.value(0,0,0);
	return;
}
void Tnode::apply_external_force(int direction,double magnitude)
{
	vec3D external_force;
	external_force.value(direction,_force.access(direction)+magnitude);
	_force=_force+external_force;
	_force0=_force0+external_force;
	return;
}
void Tnode::apply_gravity(vec3D &gravity)
{
	_force=_force+gravity*_mass;
	_force0=_force0+gravity*_mass;
}
void Tnode::apply_gravity_MPM(vec3D &gravity)
{
	if (_flag==1)
	{
		_force=_force+gravity*_mass;
		_force0=_force0+gravity*_mass;
	}
}
void Tnode::accumulate_inertance(double mass)
{
	_mass=_mass+mass;
	_mass0=_mass0+mass;
	return;
}
void Tnode::accumulate_moment(vec3D moment)
{
	_moment=_moment+moment;
	return;
}
void Tnode::initialize()
{
	_position=_position0;       
	_velocity.value(0,0,0);_velocity_temp.value(0,0,0); 
	_mass=_mass0=0;;       
	_moment.value(0,0,0);        
	_force.value(0,0,0);_force0.value(0,0,0);        
	_acceleration.value(0,0,0); 
	return;
}
void Tnode::access_nodal_bc(int (&bc_type_pos)[3],double (&bc_vel)[3])
{
	for (int i = 0; i < 3; i++)
	{
		bc_type_pos[i]=_bc_type_position[i];
		bc_vel[i]=_bc_velocity[i];
	}
	return;
}
//----------------------------------------------------------------------------------------
//Functions for Tnode_rotation
//----------------------------------------------------------------------------------------
Tnode_rotation::Tnode_rotation(int id,double x,double y,double z) : Tnode(id,x,y,z)
{
	_Omega.value(0,0,0);
	_epsilon.value(0,0,0);
	_mass_rotation=_mass_rotation0=0;
	_moment.value(0,0,0);_moment0.value(0,0,0);
	for (int i = 0; i < 3; i++)
	{
		_bc_Omega[i]=0;
		_bc_type_rotation[i]=0;
	}
}
void Tnode_rotation::calculate_acceleration()
{
	Tnode::calculate_acceleration();
	_epsilon=_moment/_mass_rotation;
	//Apply the boundary condition
	for (int i = 0; i < 3; i++)
	{
		if (_bc_type_rotation[i]!=0)
		{
			_epsilon.value(i+1,0);
		}
	}
	//Clear the nodal force after calculate acceleration
	clear_force();
	return;
}
void Tnode_rotation::update_vel(double dt)
{
	Tnode::update_vel(dt);
	_Omega_temp=_Omega;
	_Omega=_Omega+_epsilon*dt;
	_Omega_temp=(_Omega_temp+_Omega)*0.5;
	return;
}
void Tnode_rotation::clear_force()
{
	//Tnode::clear_force();
	_moment.value(0,0,0);_moment0.value(0,0,0);
	return;
}
void Tnode_rotation::accumulate_inertance(double mass,double mass_rot)
{
	Tnode::accumulate_inertance(mass);
	_mass_rotation=_mass_rotation+mass_rot;
	_mass_rotation0=_mass_rotation0+mass_rot;
	return;
}