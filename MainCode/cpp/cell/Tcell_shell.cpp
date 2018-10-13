#include "Tcell_shell.h"
#include "public_function.h"
#include <cmath>
#include "public_function.h"
#include "Gauss_weight.h"

Tcell_shell::Tcell_shell(int id,int nGauss,Tnode_rotation* (&node_ptr)[4],double h,TMAT_base** mat) : Tcell_base(id,nGauss)
{
	//Allocate memory
	_corner_stress_force.resize(4);
	//Set cell node
	for (int i = 0; i < 4; i++)
	{
		_node_ptr[i]=node_ptr[i];
	}
	//Set Gauss point material
	for (int i = 0; i < nGauss; i++)
	{
		_gausspoint[i]=mat[i];
	}
	//Set geometry information
	_h=h;
	_area=calculate_area();
	_clength=calculate_characteristic_length();
	_soundspeed_max=_gausspoint[0]->G_soundspeed();
	//Initialize variables
	for (int i = 0; i < 5; i++)
	{
		_HGstress[i]=0;
	}
	_strain_energy=0;
	return;
}
double Tcell_shell::calculate_area()
{
	double triangle_area,total_area=0,temp;
	for (int i = 0; i < 2; i++)
	{
		area(_node_ptr[0]->_position,_node_ptr[i+1]->_position,_node_ptr[i+2]->_position,triangle_area,temp);
		total_area=total_area+triangle_area;
	}
	return total_area;
}
double Tcell_shell::calculate_characteristic_length()
{
	double side_length[4],diagonal_length[2];
	for (int i = 0; i < 4; i++)
	{
		side_length[i]=(_node_ptr[(i+1)%4]->_position-_node_ptr[i]->_position).get_length();
	}
	diagonal_length[0]=(_node_ptr[0]->_position-_node_ptr[2]->_position).get_length();
	diagonal_length[1]=(_node_ptr[1]->_position-_node_ptr[3]->_position).get_length();
	double l_max;
	l_max=maxval(maxval(side_length[0],side_length[1]),maxval(side_length[2],side_length[3]));
	l_max=maxval(l_max,maxval(diagonal_length[0],diagonal_length[1]));
	return _area/l_max;
}
void Tcell_shell::assemble_corner_inertance()
{
	double mass=_area*_h*_gausspoint[0]->G_density();
	double mass_rot=maxval(_h*_h/12.0,_area/8.0)*mass;
	for (int i = 0; i < 4; i++)
	{
		_node_ptr[i]->accumulate_inertance(mass/4.0,mass_rot/4.0);
	}
	return;
}
void Tcell_shell::calculate_corner_force(Tviscosity_base* Arti_Vis,double dt)
{
	auxiliary_solver tool;
	tool.calculate_local_system(this);
	tool.calculate_local_coordinate(this);
	tool.calculate_geometry_property(this);
	tool.calculate_component_velocity_gratient();
	//Update the state at Gauss point and calculate nodal force by stress and artificial viscosity
	int num_Gauss=_gausspoint.size();
	double zeta,weight;
	double de[6]={0},vort[3]={0};         //The strain rate rotation rate
	double h_new=0;        //The new cell thickness
	_soundspeed_max=0;   //Initialize the maximal soundspeed
	for (int i = 0; i < num_Gauss; i++)
	{
		zeta=Gauss_coordinate[num_Gauss-1][i];
		weight=Gauss_weight[num_Gauss-1][i];
		double z=_h*0.5*zeta;
		//Calculate the strain rate
		calculate_strain_rate(z,tool,de);
		//Update stress at Gausspoint
		_gausspoint[i]->update_state(de,vort,dt);
		//Calculate the artificial viscosity term and modify soundspeed
		double vrate=de[0]+de[1]+de[2];
		double q=Arti_Vis->calculate_q(_gausspoint[i]->G_density(),vrate,_gausspoint[i]->G_soundspeed(),_clength);
		//q=0;
		double soundspeed=Arti_Vis->modify_soundspeed(vrate,_gausspoint[i]->G_soundspeed(),_clength);
		_gausspoint[i]->S_soundspeed(soundspeed);
		//Calculate the integral on stress
		_gausspoint[i]->calculate_stress_integral_for_shell(tool.integral_sigma,tool.integral_z_sigma,weight*_h*0.5,z,q);
		//Update the cell thickness
		h_new=h_new+_h*(1+de[2]*dt)*weight/2.0;
		//Calculate the maximal soundspeed
		if (_gausspoint[i]->G_soundspeed()>_soundspeed_max)
		{
			_soundspeed_max=_gausspoint[i]->G_soundspeed();
		}
		//Calculate the strain energy
		double gauss_volume=_area*_h*(1+de[2]*dt)*weight/2.0;
		_strain_energy=_strain_energy+_gausspoint[i]->calculate_strain_energy_increment_for_shell(de,q,gauss_volume,dt);
	}
	//Update the shell thickness
	_h=h_new;
	//Calculate the corner force and moment at local coordinate system
	vec3D local_force[4],local_moment[4];
	calculate_local_corner_force(tool,local_force,local_moment);
	//Add the hourglass force
	calculate_hourglass_viscosity(tool,dt,local_force,local_moment);
	//Calcualte the corner in global system
	calculate_global_corner_force(tool,local_force,local_moment);
	return;
}
void Tcell_shell::calculate_strain_rate(double z,const auxiliary_solver &tool, double (&de)[6])
{
	de[0]=tool.Bvelocity[0][0]+z*tool.BOmega[0][1];
	de[1]=tool.Bvelocity[1][1]-z*tool.BOmega[1][0];
	de[3]=tool.Bvelocity[1][2]-tool.NOemga[0];
	de[4]=tool.Bvelocity[0][2]+tool.NOemga[1];
	de[5]=tool.Bvelocity[1][0]+tool.Bvelocity[0][1]+z*(tool.BOmega[1][1]-tool.BOmega[0][0]);
	//de[2] will be calculated to ensure sigma[z][z]=0 in stress update
	de[2]=0;
	return;
}
void Tcell_shell::calculate_local_corner_force(auxiliary_solver &tool,vec3D (&local_force)[4],vec3D (&local_moment)[4])
{
	double y24=tool.y24,y31=tool.y31,x42=tool.x42,x13=tool.x13;
	double mx0=x42*tool.integral_z_sigma[1]+y24*tool.integral_z_sigma[2];
	double mx1=x13*tool.integral_z_sigma[1]+y31*tool.integral_z_sigma[2];
	double my0=-y24*tool.integral_z_sigma[0]-x42*tool.integral_z_sigma[2];
	double my1=-y31*tool.integral_z_sigma[0]-x13*tool.integral_z_sigma[2];
	double fx=0.25*tool.integral_sigma[4]*_area;
	double fy=0.25*tool.integral_sigma[5]*_area;
	//Calculate local corner force
	local_force[0].x=(y24*tool.integral_sigma[0]+x42*tool.integral_sigma[3])*0.5;
	local_force[0].y=(x42*tool.integral_sigma[1]+y24*tool.integral_sigma[3])*0.5;
	local_force[0].z=(y24*tool.integral_sigma[5]+x42*tool.integral_sigma[4])*0.5;

	local_force[1].x=(y31*tool.integral_sigma[0]+x13*tool.integral_sigma[3])*0.5;
	local_force[1].y=(x13*tool.integral_sigma[1]+y31*tool.integral_sigma[3])*0.5;
	local_force[1].z=(y31*tool.integral_sigma[5]+x13*tool.integral_sigma[4])*0.5;

	local_force[2].x=-local_force[0].x;
	local_force[2].y=-local_force[0].y;
	local_force[2].z=-local_force[0].z;

	local_force[3].x=-local_force[1].x;
	local_force[3].y=-local_force[1].y;
	local_force[3].z=-local_force[1].z;

	//Calculate local corner moment
	local_moment[0].x=mx0*0.5-fx;local_moment[0].y=my0*0.5+fy;
	local_moment[1].x=mx1*0.5-fx;local_moment[1].y=my1*0.5+fy;
	local_moment[2].x=-mx0*0.5-fx;local_moment[2].y=-my0*0.5+fy;
	local_moment[3].x=-mx1*0.5-fx;local_moment[3].y=-my1*0.5+fy; 
	local_moment[0].z=local_moment[1].z=local_moment[2].z=local_moment[3].z=0;
	return;
}
void Tcell_shell::calculate_hourglass_viscosity(auxiliary_solver &tool,double dt,vec3D (&local_force)[4],vec3D (&local_moment)[4])
{
	//Calculate the coefficient
	double E=_gausspoint[0]->G_Youngs(),G=_gausspoint[0]->G_G();
	double tmode,wmode,mmode;
	tmode=0.02*E/8.0;wmode=0.02*G/12.0;mmode=0.02*E/192.0;
	//Calculate gamma[i]
	double Invarea=1.0/_area;
	double htx=Invarea*(-tool.local_position[1].x+tool.local_position[2].x-tool.local_position[3].x);
	double hty=Invarea*(-tool.local_position[1].y+tool.local_position[2].y-tool.local_position[3].y);
	double gamma[4];
	gamma[0]=1-htx*tool.y24-hty*tool.x42;
	gamma[1]=-1-htx*tool.y31-hty*tool.x13;
	gamma[2]=2-gamma[0];gamma[3]=-2-gamma[1];
	//Calculate the hourglass rate
	double HGrate[5]={0};
	for (int i = 0; i < 4; i++)
	{
		HGrate[0]=HGrate[0]+gamma[i]*tool.local_velocity[i].x;
		HGrate[1]=HGrate[1]+gamma[i]*tool.local_velocity[i].y;
		HGrate[2]=HGrate[2]+gamma[i]*tool.local_velocity[i].z;
		HGrate[3]=HGrate[3]+gamma[i]*tool.local_Omega[i].x;
		HGrate[4]=HGrate[4]+gamma[i]*tool.local_Omega[i].y;
	}
	//Calculate hourglass stress
	double coef[5];
	double BB=tool.calculate_sumb();
	double fac1=Invarea*BB*_h*2.0;
	double fac2=fac1*_h*_h;
	coef[0]=coef[1]=tmode*fac1;
	coef[2]=wmode*fac2;
	coef[3]=coef[4]=mmode*fac2;
	//Update the hourglass stress
	for (int i = 0; i < 5; i++)
	{
		_HGstress[i]=_HGstress[i]+coef[i]*HGrate[i]*dt;
	}
	//Add hourglass force into corner force
	for (int i = 0; i < 4; i++)
	{
		local_force[i].x=local_force[i].x+gamma[i]*_HGstress[0];
		local_force[i].y=local_force[i].y+gamma[i]*_HGstress[1];
		local_force[i].z=local_force[i].z+gamma[i]*_HGstress[2];
		local_moment[i].x=local_moment[i].x+gamma[i]*_HGstress[3];
		local_moment[i].y=local_moment[i].y+gamma[i]*_HGstress[4];
	}
	return;
}
void Tcell_shell::calculate_global_corner_force(auxiliary_solver &tool,vec3D (&local_force)[4],vec3D (&local_moment)[4])
{
	for (int i = 0; i < 4; i++)
	{
		_corner_stress_force[i]=-local_force[i].multiply_by_matrix(tool.transfer_matrix);
		_corner_stress_moment[i]=-local_moment[i].multiply_by_matrix(tool.transfer_matrix);
	}
	return;
}
void Tcell_shell::assemble_corner_force()
{
	for (int i = 0; i < 4; i++)
	{
		_node_ptr[i]->_force=_node_ptr[i]->_force+_corner_stress_force[i];
		_node_ptr[i]->_force0=_node_ptr[i]->_force0+_corner_stress_force[i];
		//_corner_stress_force[i].value(0,0,0);
		_node_ptr[i]->_moment=_node_ptr[i]->_moment+_corner_stress_moment[i];
		_node_ptr[i]->_moment0=_node_ptr[i]->_moment0+_corner_stress_moment[i];
		//_corner_stress_moment[i].value(0,0,0);
	}
	return;
}
double Tcell_shell::calculate_time_step()
{
	return _clength/_soundspeed_max;
}
void Tcell_shell::update_strain_energy(double dt)
{
	double de=0;
	for (int i = 0; i < 4; i++)
	{
		de=de+_corner_stress_force[i]*_node_ptr[i]->_velocity+_corner_stress_moment[i]*_node_ptr[i]->_Omega;
		_corner_stress_force[i].value(0,0,0);
	}
	_strain_energy=_strain_energy+de*dt;
	return;
}
//------------------------------------------------------------------------------
//Functions for auxiliary solver
//------------------------------------------------------------------------------
Tcell_shell::auxiliary_solver::auxiliary_solver()
{
	for (int i = 0; i < 6; i++)
	{
		integral_sigma[i]=0;
	}
	for (int i = 0; i < 3; i++)
	{
		integral_z_sigma[i]=0;
	}
}
void Tcell_shell::auxiliary_solver::calculate_local_system(Tcell_shell*	this_cell)
{
	r21=this_cell->_node_ptr[1]->_position-this_cell->_node_ptr[0]->_position;
	r31=this_cell->_node_ptr[2]->_position-this_cell->_node_ptr[0]->_position;
	r41=this_cell->_node_ptr[3]->_position-this_cell->_node_ptr[0]->_position;
	r42=this_cell->_node_ptr[3]->_position-this_cell->_node_ptr[1]->_position;
	r34=this_cell->_node_ptr[2]->_position-this_cell->_node_ptr[3]->_position;
	//Calculate e3
	vec3D v1,v2,v_temp1,v_temp2;
	v1=r21+r34;v2=r31+r42;
	v1.normalize();v2.normalize();
	e3=v1%v2;e3.normalize();
	//Calculate e1 and e2
	v_temp1=(v1+v2)*0.5;v_temp1.normalize();
	v_temp2=e3%v_temp1;
	e1=(v_temp1-v_temp2)/sqrt(2);
	e2=(v_temp1+v_temp2)/sqrt(2);
	transfer_matrix[0][0]=e1.x;transfer_matrix[0][1]=e2.x;transfer_matrix[0][2]=e3.x;
	transfer_matrix[1][0]=e1.y;transfer_matrix[1][1]=e2.y;transfer_matrix[1][2]=e3.y;
	transfer_matrix[2][0]=e1.z;transfer_matrix[2][1]=e2.z;transfer_matrix[2][2]=e3.z;
	return;
}
void Tcell_shell::auxiliary_solver::calculate_local_coordinate(Tcell_shell* this_cell)
{
	//Calculate node coordinate in local system
	//The z-component of local coordinate is zero
	local_position[0].value(0,0,0);
	local_position[1].value(e1*r21,e2*r21,0);
	local_position[2].value(e1*r31,e2*r31,0);
	local_position[3].value(e1*r41,e2*r41,0);
	//Calculate node velocity in local system
	for (int i = 0; i < 4; i++)
	{
		vec3D temp=this_cell->_node_ptr[i]->_velocity;
		local_velocity[i].value(e1*temp,e2*temp,e3*temp);
	}
	//Calculate the node rotation rate in local system
	//The z-component of local Omega is zero
	for (int i = 0; i < 4; i++)
	{
		vec3D temp=this_cell->_node_ptr[i]->_Omega;
		local_Omega[i].value(e1*temp,e2*temp,0);
	}
	//Calculate the component for B materix
	y24=local_position[1].y-local_position[3].y;
	y31=local_position[2].y;
	x42=local_position[3].x-local_position[1].x;
	x13=-local_position[2].x;
	return;
}
void Tcell_shell::auxiliary_solver::calculate_geometry_property(Tcell_shell* this_cell)
{
	//Calculate the area of cell
	this_cell->_area=0.5*(x13*y24-x42*y31);
	invJ=0.5/this_cell->_area;
	//Calculate the diagonal and side length of the cell
	double diag[2],side[4];
	diag[0]=x13*x13+y31*y31;diag[1]=x42*x42+y24*y24;
	for (int i = 0; i < 4; i++)
	{
		side[i]=(local_position[(i+1)%4]-local_position[i]).self_multuply();
	}
	double side_max=maxval(maxval(side[0],side[1]),maxval(side[2],side[3]));
	double diag_max=maxval(diag[0],diag[1]);
	double l_max=maxval(side_max,diag_max);
	this_cell->_clength=this_cell->_area/sqrt(l_max);
	return;
}
void Tcell_shell::auxiliary_solver::calculate_component_velocity_gratient()
{
	vec3D v13=(local_velocity[0]-local_velocity[2])*invJ;
	vec3D v24=(local_velocity[1]-local_velocity[3])*invJ;
	vec3D O13=(local_Omega[0]-local_Omega[2])*invJ;
	vec3D O24=(local_Omega[1]-local_Omega[3])*invJ;
	//Calculate the component for velocity gradient
	Bvelocity[0][0]=y24*v13.x+y31*v24.x;
	Bvelocity[0][1]=y24*v13.y+y31*v24.y;
	Bvelocity[0][2]=y24*v13.z+y31*v24.z;
	Bvelocity[1][0]=x42*v13.x+x13*v24.x;
	Bvelocity[1][1]=x42*v13.y+x13*v24.y;
	Bvelocity[1][2]=x42*v13.z+x13*v24.z;

	BOmega[0][0]=y24*O13.x+y31*O24.x;
	BOmega[0][1]=y24*O13.y+y31*O24.y;
	BOmega[1][0]=x42*O13.x+x13*O24.x;
	BOmega[1][1]=x42*O13.y+x13*O24.y;

	vec3D temp=(local_Omega[0]+local_Omega[1]+local_Omega[2]+local_Omega[3])*0.25;
	NOemga[0]=temp.x;NOemga[1]=temp.y;

	return;
}