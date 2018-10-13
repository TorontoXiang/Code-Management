#include "Tcell_mixed_fluid_MPM.h"
#include "Shape_function.h"
#include "data_structure.h"
double Tcell_mixed_fluid_MPM::calculate_time_step_particle()
{
	if (!_is_empty)
	{
		double max_node_speed=0;
		for (int i = 0; i < 8; i++)
		{
			max_node_speed=maxval(_node_ptr[i]->_velocity.get_length(),max_node_speed);
		}
		return _cell_length/(_maxspd_particle+max_node_speed);
	}
	else
	{
		return 1e10;
	}
}
void Tcell_mixed_fluid_MPM::distribute_particle_force(TMPM_particle* particle_ptr)
{
	//Apply the internal force
	double stress[6],fluid_pressure;
	particle_ptr->_gausspoint->G_stress(stress);
	fluid_pressure=_gausspoint[0]->G_pressure();
	double sxx=stress[0]-fluid_pressure,syy=stress[1]-fluid_pressure,szz=stress[2]-fluid_pressure,sxy=stress[3],sxz=stress[4],syz=stress[5];
	vec3D standard_coordinate=particle_ptr->_standard_coordinate;
	//Calculate the shape function at the particle
	vec3D GNshape[8];
	double detJ;
	calculate_DN_matrix(GNshape,detJ,standard_coordinate.x,standard_coordinate.y,standard_coordinate.z);
	double particle_volume=particle_ptr->_volume;
	for (int i = 0; i < 8; i++)
	{
		//_corner_stress_force_particle[i].x=_corner_stress_force_particle[i].x-(sxx*GNshape[i].x+sxy*GNshape[i].y+sxz*GNshape[i].z)*particle_volume;
		//_corner_stress_force_particle[i].y=_corner_stress_force_particle[i].y-(sxy*GNshape[i].x+syy*GNshape[i].y+syz*GNshape[i].z)*particle_volume;
		//_corner_stress_force_particle[i].z=_corner_stress_force_particle[i].z-(sxz*GNshape[i].x+syz*GNshape[i].y+szz*GNshape[i].z)*particle_volume;
		_corner_stress_force[i].x=_corner_stress_force[i].x-(sxx*GNshape[i].x+sxy*GNshape[i].y+sxz*GNshape[i].z)*particle_volume;
		_corner_stress_force[i].y=_corner_stress_force[i].y-(sxy*GNshape[i].x+syy*GNshape[i].y+syz*GNshape[i].z)*particle_volume;
		_corner_stress_force[i].z=_corner_stress_force[i].z-(sxz*GNshape[i].x+syz*GNshape[i].y+szz*GNshape[i].z)*particle_volume;
	}
	//Apply the external force
	double Nshape[8];
	for (int i = 0; i < 8; i++)
	{
		Nshape[i]=N_brick_8(i,standard_coordinate.x,standard_coordinate.y,standard_coordinate.z);
		//_corner_stress_force_particle[i]=_corner_stress_force_particle[i]+particle_ptr->_external_force*Nshape[i];
		_corner_stress_force[i]=_corner_stress_force[i]+particle_ptr->_external_force*Nshape[i];
	}
	return;
}
void Tcell_mixed_fluid_MPM::update_particle_DOF(TMPM_particle* particle_ptr,double dt)
{
	vec3D standard_coordiante=particle_ptr->_standard_coordinate;
	double Nshape[8];
	//Update the particle position
	particle_ptr->_position.value(0,0,0);
	//particle_ptr->_velocity.value(0,0,0);
	vec3D acc;
	for (int i = 0; i < 8; i++)
	{
		Nshape[i]=N_brick_8(i,standard_coordiante.x,standard_coordiante.y,standard_coordiante.z);
		particle_ptr->_position=particle_ptr->_position+_node_ptr[i]->_position*Nshape[i];
		//particle_ptr->_velocity=particle_ptr->_velocity+_node_ptr[i]->_velocity*Nshape[i];
		acc=acc+_node_ptr[i]->_acceleration*Nshape[i];
	}
	//Update the particle velocity
	particle_ptr->_velocity=particle_ptr->_velocity+acc*dt;
}
void Tcell_mixed_fluid_MPM::update_particle_stress(TMPM_particle* particle_ptr,double dt)
{
	//Calculate the shape function at the particle
	vec3D standard_coordinate=particle_ptr->_standard_coordinate;
	vec3D GNshape[8];
	double detJ;
	calculate_DN_matrix(GNshape,detJ,standard_coordinate.x,standard_coordinate.y,standard_coordinate.z);
	double de[6]={0},vort[3]={0};
	//de[0]:ux,x;de[1]:uy,y;de[2]:uz,z;de[3]:uy,z+uz,y;de[4]:uz,x+ux,z;de[5]:ux,y+uy,x;
	//vort[0]:uz,y-uy,z;vort[1]:ux,z-uz,x;vort[2]:uy,x-ux,y
	//Calculate the strain rate and the vortex rate
	for (int i = 0; i < 8; i++)
	{
		de[0]=de[0]+_node_ptr[i]->_velocity.x*GNshape[i].x;
		de[1]=de[1]+_node_ptr[i]->_velocity.y*GNshape[i].y;
		de[2]=de[2]+_node_ptr[i]->_velocity.z*GNshape[i].z;
		de[3]=de[3]+_node_ptr[i]->_velocity.y*GNshape[i].z+_node_ptr[i]->_velocity.z*GNshape[i].y;
		de[4]=de[4]+_node_ptr[i]->_velocity.z*GNshape[i].x+_node_ptr[i]->_velocity.x*GNshape[i].z;
		de[5]=de[5]+_node_ptr[i]->_velocity.x*GNshape[i].y+_node_ptr[i]->_velocity.y*GNshape[i].x;
	}
	for (int i = 0; i < 6; i++)
	{
		particle_ptr->_de[i]=de[i];
	}
	particle_ptr->_detJ=detJ;
	//cout<<"de: "<<de[0]<<" "<<de[1]<<" "<<de[2]<<" "<<de[3]<<" "<<de[4]<<" "<<de[5]<<endl;
	for (int i = 0; i < 8; i++)
	{
		vort[0]=vort[0]+_node_ptr[i]->_velocity.z*GNshape[i].y-_node_ptr[i]->_velocity.y*GNshape[i].z;
		vort[1]=vort[1]+_node_ptr[i]->_velocity.x*GNshape[i].z-_node_ptr[i]->_velocity.z*GNshape[i].x;
		vort[2]=vort[2]+_node_ptr[i]->_velocity.y*GNshape[i].x-_node_ptr[i]->_velocity.x*GNshape[i].y;
	}
	vort[0]=vort[0]/2.0;vort[1]=vort[1]/2.0;vort[2]=vort[2]/2.0;
	//Update the stress, density and soundspeed on the particle
	particle_ptr->_gausspoint->update_state(de,vort,dt);
	//Update the particle volume
	particle_ptr->_volume=particle_ptr->_mass/particle_ptr->_gausspoint->G_density();
	return;
}
void Tcell_mixed_fluid_MPM::distribute_particle_variables(TMPM_particle* particle_ptr,vec3D interval)
{
	_is_empty=false;
	vec3D standard_coordinate=calculate_standard_coordinate(particle_ptr->_position,interval);
	particle_ptr->_standard_coordinate=standard_coordinate;
	//Calculate the shape function at the particle
	double Nshape[8];
	for (int i = 0; i < 8; i++)
	{
		Nshape[i]=N_brick_8(i,standard_coordinate.x,standard_coordinate.y,standard_coordinate.z);
	}
	//Distribute the particle mass and moment
	for (int i = 0; i < 8; i++)
	{
		double mass=particle_ptr->_mass-_gausspoint[0]->G_density()*particle_ptr->_volume;
		_node_ptr[i]->accumulate_inertance(mass*Nshape[i]);
		_node_ptr[i]->accumulate_moment(particle_ptr->_velocity*mass*Nshape[i]);
		_node_ptr[i]->_flag=1;
	}
	return;
}
vec3D Tcell_mixed_fluid_MPM::calculate_standard_coordinate(vec3D& particle_coordinate,vec3D interval)
{
	vec3D standard_coordinate;
	standard_coordinate.x=2*(particle_coordinate-_node_ptr[0]->_position).x/interval.x-1;
	standard_coordinate.y=2*(particle_coordinate-_node_ptr[0]->_position).y/interval.y-1;
	standard_coordinate.z=2*(particle_coordinate-_node_ptr[0]->_position).z/interval.z-1;
	if (abs(standard_coordinate.x)>1 || abs(standard_coordinate.y)>1 || abs(standard_coordinate.z)>1)
	{
		cout<<particle_coordinate.x<<" "<<particle_coordinate.y<<" "<<particle_coordinate.z<<endl;
		cout<<standard_coordinate.x<<" "<<standard_coordinate.y<<" "<<standard_coordinate.z<<endl;
		cout<<_node_ptr[0]->_position.x<<" "<<_node_ptr[0]->_position.y<<" "<<_node_ptr[0]->_position.z<<endl;
		cout<<"Error: The standard coordiante of a particle excess!"<<endl;
		system("Pause");
		exit(0);
	}
	return standard_coordinate;
}
void Tcell_mixed_fluid_MPM::reset_state_after_remapping()
{
	Tcell_mixed_fluid::reset_state_after_remapping();
	_is_empty=true;
	return;
}
void Tcell_mixed_fluid_MPM::assemble_corner_force()
{
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->_force=_node_ptr[i]->_force+_corner_stress_force[i]+_corner_hourglass_force[i]+_corner_viscosity_force[i]+_corner_stress_force_particle[i];
		_node_ptr[i]->_force0=_node_ptr[i]->_force0+_corner_stress_force[i]+_corner_hourglass_force[i]+_corner_viscosity_force[i]+_corner_stress_force_particle[i];
		_corner_stress_force_particle[i].value(0,0,0);
	}
	return;
}