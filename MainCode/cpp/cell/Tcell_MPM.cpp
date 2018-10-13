#include "Tcell_MPM.h"
#include "Shape_function.h"
#include "public_function.h"
#include <cmath>
vec3D Tcell_MPM::_length;
void Tcell_MPM::calculate_node_ptr(Tnode* node_ptr[8])
{
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]=node_ptr[i];
	}
	return;
}
void Tcell_MPM::distribute_particle_to_node(TMPM_particle* particle_ptr)
{
	_is_empty=false;
	vec3D standard_coordinate=calculate_standard_coordinate(particle_ptr->_position);
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
		_node_ptr[i]->accumulate_inertance(particle_ptr->_mass*Nshape[i]);
		_node_ptr[i]->accumulate_moment(particle_ptr->_velocity*particle_ptr->_mass*Nshape[i]);
	}
	return;
}
vec3D Tcell_MPM::calculate_standard_coordinate(vec3D& particle_coordinate)
{
	vec3D standard_coordinate;
	standard_coordinate.x=2*(particle_coordinate-_node_ptr[0]->_position).x/_length.x-1;
	standard_coordinate.y=2*(particle_coordinate-_node_ptr[0]->_position).y/_length.y-1;
	standard_coordinate.z=2*(particle_coordinate-_node_ptr[0]->_position).z/_length.z-1;
	if (abs(standard_coordinate.x)>1 || abs(standard_coordinate.y)>1 || abs(standard_coordinate.z)>1)
	{
		cout<<"Error: The standard coordiante of a particle excess!"<<endl;
		system("Pause");
		exit(0);
	}
	return standard_coordinate;
}
void Tcell_MPM::initialize()
{
	_max_spd=-1;
	for (int i = 0; i < 8; i++)
	{
		_corner_force[i].value(0,0,0);
	}
	_clength=0;
	_is_empty=true;
	return;
}
void Tcell_MPM::calculate_cell_length()
{
	double l,l_min=100000;
	for (int i = 0; i < 4; i++)
	{
		l=(_node_ptr[i]->_position-_node_ptr[(i+1)%4]->_position).self_multuply();
		l_min=minval(l,l_min);
		l=(_node_ptr[i+4]->_position-_node_ptr[(i+1)%4+4]->_position).self_multuply();
		l_min=minval(l,l_min);
		l=(_node_ptr[i]->_position-_node_ptr[i+4]->_position).self_multuply();
		l_min=minval(l,l_min);
	}
	_clength=sqrt(l_min);
	return;
}
double Tcell_MPM::calculate_cell_time_step(double& clength,double& soundspeed,double& nodespeed)
{
	if (!_is_empty)
	{
		double max_node_speed=0;
		for (int i = 0; i < 8; i++)
		{
			max_node_speed=maxval(_node_ptr[i]->_velocity.get_length(),max_node_speed);
		}
		clength=_clength;soundspeed=_max_spd;nodespeed=max_node_speed;
		return _clength/(_max_spd+max_node_speed);
	}
	else
	{
		return 1e10;
	}
}
void Tcell_MPM::distribute_particle_force(TMPM_particle* particle_ptr)
{
	//Apply the internal force
	double stress[6];
	particle_ptr->_gausspoint->G_stress(stress);
	double sxx=stress[0],syy=stress[1],szz=stress[2],sxy=stress[3],sxz=stress[4],syz=stress[5];
	vec3D standard_coordinate=particle_ptr->_standard_coordinate;
	//Calculate the shape function at the particle
	vec3D GNshape[8];
	double detJ;
	calculate_DN_matrix(GNshape,detJ,standard_coordinate.x,standard_coordinate.y,standard_coordinate.z);
	double particle_volume=particle_ptr->_volume;
	for (int i = 0; i < 8; i++)
	{
		_corner_force[i].x=_corner_force[i].x-(sxx*GNshape[i].x+sxy*GNshape[i].y+sxz*GNshape[i].z)*particle_volume;
		_corner_force[i].y=_corner_force[i].y-(sxy*GNshape[i].x+syy*GNshape[i].y+syz*GNshape[i].z)*particle_volume;
		_corner_force[i].z=_corner_force[i].z-(sxz*GNshape[i].x+syz*GNshape[i].y+szz*GNshape[i].z)*particle_volume;
	}
	//Apply the external force
	double Nshape[8];
	for (int i = 0; i < 8; i++)
	{
		Nshape[i]=N_brick_8(i,standard_coordinate.x,standard_coordinate.y,standard_coordinate.z);
		_corner_force[i]=_corner_force[i]+particle_ptr->_external_force*Nshape[i];
	}
	return;
}
void Tcell_MPM::assemble_corner_force()
{
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->_force=_node_ptr[i]->_force+_corner_force[i];
		_corner_force[i].value(0,0,0);
	}
	return;
}
void Tcell_MPM::update_particle_state(TMPM_particle* particle_ptr,double dt)
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
	//cout<<"de: "<<de[0]<<" "<<de[1]<<" "<<de[2]<<" "<<de[3]<<" "<<de[4]<<" "<<de[5]<<endl;
	for (int i = 0; i < 8; i++)
	{
		vort[0]=vort[0]+_node_ptr[i]->_velocity.z*GNshape[i].y-_node_ptr[i]->_velocity.y*GNshape[i].z;
		vort[1]=vort[1]+_node_ptr[i]->_velocity.x*GNshape[i].z-_node_ptr[i]->_velocity.z*GNshape[i].x;
		vort[2]=vort[2]+_node_ptr[i]->_velocity.y*GNshape[i].x-_node_ptr[i]->_velocity.x*GNshape[i].y;
	}
	vort[0]=vort[0]/2.0;vort[1]=vort[1]/2.0;vort[2]=vort[2]/2.0;
	//Update the stress, density and soundspeed on the particle
	particle_ptr->_gausspoint->update_state(de,vort,dt,particle_ptr->_volume);
	//Update the particle volume
	particle_ptr->_volume=particle_ptr->_mass/particle_ptr->_gausspoint->G_density();
	return;
}
void Tcell_MPM::update_particle_DOF(TMPM_particle* particle_ptr,double dt)
{
	vec3D standard_coordiante=particle_ptr->_standard_coordinate;
	double Nshape[8];
	//Update the particle position
	particle_ptr->_position.value(0,0,0);
	vec3D acc;
	for (int i = 0; i < 8; i++)
	{
		Nshape[i]=N_brick_8(i,standard_coordiante.x,standard_coordiante.y,standard_coordiante.z);
		particle_ptr->_position=particle_ptr->_position+_node_ptr[i]->_position*Nshape[i];
		acc=acc+_node_ptr[i]->_acceleration*Nshape[i];
	}
	//Update the particle velocity
	particle_ptr->_velocity=particle_ptr->_velocity+acc*dt;
	return;
}
void Tcell_MPM::calculate_DN_matrix(vec3D (&DN)[8],double &detJ,double xi,double eta,double zeta)
{
	double Jacobi[3][3];
	vec3D GN[8];
	bool degenerate;
	calcualte_Jacobi(Jacobi,GN,xi,eta,zeta);
	inverse(Jacobi,detJ,degenerate);
	for (int i = 0; i < 8; i++)
	{
		DN[i]=GN[i].multiply_by_matrix(Jacobi);
	}
	return;
}
void Tcell_MPM::calcualte_Jacobi(double (&Jacobi)[3][3],vec3D (&GN)[8],double xi,double eta,double zeta)
{
	for (int i = 0; i < 8; i++)
	{
		GN[i]=GN_brick_8(i,xi,eta,zeta);
	}
	vec3D temp[3];
	//Calculate the Jacobi matrix with respected to the current position
	for (int i=0;i<8;i++)
	{
		for (int j = 0; j < 3; j++)
		{
			temp[j]=temp[j]+GN[i]*_node_ptr[i]->_position.access(j+1);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Jacobi[i][j]=temp[j].access(i+1);
		}
	}
	return;
}