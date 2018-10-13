#include "Tnode_fluid.h"
#include "public_function.h"
Tnode_fluid::Tnode_fluid(int id,double x,double y,double z) : Tnode(id,x,y,z)
{
	_position_temp.value(0,0,0);     
	_density=_soundspeed=_pressure=0;
	_moment.value(0,0,0);
	_flag=0;
}
Tnode_fluid::Tnode_fluid(int id,vec3D coordinate) : Tnode(id,coordinate.x,coordinate.y,coordinate.z)
{
	_position_temp.value(0,0,0);     
	_density=_soundspeed=_pressure=0;
	_moment.value(0,0,0);
	_flag=0;
}
void Tnode_fluid::predict_position(double dt)
{
	_position_temp=_position+_velocity*dt*0.5;
	return;
}
void Tnode_fluid::add_connected_cell(Tcell_fluid_base* connected_cell)
{
	_surroundings.push_back(connected_cell);
	return;
}
void Tnode_fluid::calculate_nodal_variables()
{
	Tcell_fluid_base* temp;
	double av_density,av_pressure,av_soundspeed,cell_volume,total_volume=0;
	//Initialize
	_density=_soundspeed=_pressure=0;
	int num_surrounding=_surroundings.size();
	for (int i=0;i<num_surrounding;i++)
	{
		temp=_surroundings[i];
		cell_volume=temp->_cell_volume;
		temp->calculate_average_variable(av_density,av_pressure,av_soundspeed);
		_density=_density+av_density*cell_volume;
		_soundspeed=_soundspeed+av_soundspeed*cell_volume;
		_pressure=_pressure+av_pressure*cell_volume;
		total_volume=total_volume+cell_volume;
	}
	_density=_density/total_volume;
	_soundspeed=_soundspeed/total_volume;
	_pressure=_pressure/total_volume;
	return;
}
void Tnode_fluid::reset_velocity_after_remapping(vec3D velocity_origin)
{
	//If it is a co-node,the velocity will be unchanged
	if (_is_conode==1)
	{
		_velocity=velocity_origin;
	}
	else
	{
		_velocity=_moment/_mass;
		_moment.value(0,0,0);
	}
	//Apply boundary conditions in remapping
	for (int i = 0; i < 3; i++)
	{
		if (_bc_type_position[i]==1)
		{
			_velocity.value(i+1,0);
		}
		if (_bc_type_position[i]==2)
		{
			_velocity.value(i+1,_bc_velocity[i]);
		}
	}
	return;
}
void Tnode_fluid::set_velocity_from_particle_after_remapping()
{

	if (_flag==1)
	{
		_velocity=_moment/_mass;
		_moment.value(0,0,0);
	}
	//Apply boundary conditions in remapping
	for (int i = 0; i < 3; i++)
	{
		if (_bc_type_position[i]==1)
		{
			_velocity.value(i+1,0);
		}
		if (_bc_type_position[i]==2)
		{
			_velocity.value(i+1,_bc_velocity[i]);
		}
	}
	return;
}
void Tnode_fluid::clear_nodal_variable()
{
	_mass=_mass0=_density=_soundspeed=_pressure=0;
	_position_temp.value(0,0,0);     
	_moment.value(0,0,0);	
	_velocity.value(0,0,0);
	_velocity_temp.value(0,0,0);
	_acceleration.value(0,0,0);
	_force.value(0,0,0);
	return;	
}
void Tnode_fluid::perturbance_node_for_implosion(vec3D center,int dimension,double ratio,double n)
{
	if (dimension==2)
	{
		vec2D rr(_position.x-center.x,_position.y-center.y);
		double l=rr.get_length();
		double theta=asin(rr.y/l);
		double l_new=l*(1-ratio*(1-cos(n*theta)));
		//double l_new=l*(1+ratio*cos(n*theta));
		_position.value(center.x+rr.x*l_new/l,center.y+rr.y*l_new/l,_position.z);
	}
	else if (dimension==3)
	{
		vec3D rr=_position-center;
		double l=rr.get_length();
		double theta=acos(rr.z/l);
		double phy;
		if (abs(theta)<1e-10)
		{
			phy=0;
		}
		else
		{
			phy=acos(minval(1.0,rr.x/l/sin(theta)));
		}
		double l_new=l*(1-ratio*(cos(n*theta)*cos(n*phy)));
		_position.value(center.x+rr.x*l_new/l,center.y+rr.y*l_new/l,center.z+rr.z*l_new/l);
	}
}
void Tnode_fluid::move_node_for_implosion(vec3D center,int dimension,int N,double d0,double r_in,double r_out,double alpha)
{
	int n;
	double l_new=0;
	if (dimension==2)
	{
		vec2D rr(_position.x-center.x,_position.y-center.y);
		double l=rr.get_length();
		double l0=(r_out-r_in)/N;
		if (l>r_in+0.5*l0 && l<r_out+l0)
		{
			double ratio=(l-r_in)/l0;		
			if (ceil(ratio)-ratio>ratio-floor(ratio))
			{
				n=floor(ratio);
			}
			else
			{
				n=ceil(ratio);
			}
			for (int i = 0; i < n; i++)
			{
				l_new=l_new+d0*pow(alpha,i);
			}
			l_new=l_new+r_in;
			_position.value(center.x+rr.x*l_new/l,center.y+rr.y*l_new/l,_position.z);
		}
	}
	else if (dimension==3)
	{
		vec3D rr=_position-center;
		double l=rr.get_length();

	}
}
void Tnode_fluid::set_velocity_bc(double t,int dimension)
{
	double v0=4,a=40;
	if (_bc_type_position[0]==3)
	{
		double v;
		v=maxval(v0-a*t,0.0);
		if (dimension==2)
		{
			vec2D r(_position.x,_position.y);
			double l=r.get_length();
			_bc_velocity[0]=-v*r.x/l;_bc_velocity[1]=-v*r.y/l;_bc_velocity[3]=0;
		}
		else if (dimension==3)
		{
			double l=_position.get_length();
			_bc_velocity[0]=-v*_position.x/l;_bc_velocity[1]=-v*_position.y/l;_bc_velocity[2]=-v*_position.z/l;
		}
	}
}