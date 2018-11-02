#include "Tnode.h"
#include <iostream>
using namespace std;
//-----------------------------------------------------------------------------------------
//Functions for Tnode_base
//-----------------------------------------------------------------------------------------
Tnode::Tnode()
{
	_id=-1;
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
	_mass=0;
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
		cout << "Error: Zero or negative nodal mass on node "<<_id<< endl;
		system("Pause");
		exit(0);
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
void Tnode::update_vel(double dt)
{
	_velocity=_velocity+_acceleration*dt;
	for (int i = 0; i < 3; i++)
	{
		if (_bc_type_position[i]==1)
		{
			_velocity.value(i+1,0);
		}
		else if (_bc_type_position[i]==2)
		{
			_velocity.value(i+1,_bc_velocity[i]);
		}
	}
	return;
}
void Tnode::update_pos(double dt)
{
	_position=_position+_velocity*dt;
	return;
}
void Tnode::clear_force()
{
	_force.value(0, 0, 0);
	return;
}
void Tnode::apply_external_force(int direction,double magnitude)
{
	vec3D external_force;
	external_force.value(direction,_force.access(direction)+magnitude);
	_force=_force+external_force;
	return;
}
void Tnode::apply_gravity(vec3D &gravity)
{
	_force=_force+gravity*_mass;
}
void Tnode::accumulate_inertance(double mass)
{
	_mass=_mass+mass;
	return;
}
void Tnode::initialize()
{
	_position=_position0;       
	_velocity.value(0,0,0);
	_mass=0;;             
	_force.value(0,0,0);     
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
Tnode_CNT::Tnode_CNT(int id, double x, double y, double z)
{
	_id = id;
	_position.value(x, y, z);
	_position0 = _position;
	_location_id = -1;
	_ix = _iy = _iz = -2;
}