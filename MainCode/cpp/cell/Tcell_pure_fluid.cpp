#include "Tcell_pure_fluid.h"
Tcell_pure_fluid::Tcell_pure_fluid(int cell_id,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat)
	:Tcell_fluid_base(cell_id,1,node_ptr,mat)
{
	if (mat!=NULL)
	{
		calculate_corner_mass();
		_material_mass=_gausspoint[0]->G_density()*_cell_volume;
	}
	_material_internal_energy=0;
	_flag = 0;
}
void Tcell_pure_fluid::reset_conservative_variable_initial_condition()
{
	calculate_corner_mass();
	_material_mass=_gausspoint[0]->G_density()*_cell_volume;
	_material_internal_energy=0;
}
bool Tcell_pure_fluid::calculate_corner_stress_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3],double dt)
{
	double temp_volume;                     //The volume of cell at predict step
	double temp_pressure;                   //The pressure at predict step
	bool negative;
	temp_volume=calculate_cell_volume("Temporary",negative);
	double ro,e;
	ro=_gausspoint[0]->G_density();e=_gausspoint[0]->G_internal_energy();
	temp_pressure=_gausspoint[0]->G_EOS()->predict_p(ro,e,temp_volume,_cell_volume);
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<3;j++)
		{
			_corner_stress_force[i]=_corner_stress_force[i]+surfacial[i][j]*temp_pressure;
		}
	}
	//calculate_viscosity_force();
	return negative;
}
void Tcell_pure_fluid::calculate_corner_hourglass_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3])
{
	double dp[8];                    //The difference of pressure due to hourglass
	for (int i = 0; i < 8; i++)
	{
		dp[i]=0;
	}
	double ro=_gausspoint[0]->G_density();
	double p=_gausspoint[0]->G_pressure();
	double e=_gausspoint[0]->G_internal_energy();
	double c=_gausspoint[0]->G_soundspeed();
	for (int i=0;i<8;i++)
	{
		dp[i]=_gausspoint[0]->G_EOS()->ro_e_p(_corner_mass[i]/_corner_volume[i],e)-p;
	}
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<3;j++)
		{
 			_corner_hourglass_force[i]=_corner_hourglass_force[i]+surfacial[i][j]*dp[i]
 			+assistant[i][j]*(dp[i]-dp[Tcell_fluid_base::cell_corner_v[i][j]])*0.5;
		}
	}
	return;
}
void Tcell_pure_fluid::calculate_corner_mass()
{
	for (int i = 0; i < 8; i++)
	{
		_corner_mass[i]=_corner_volume[i]*_gausspoint[0]->G_density();
	}
	return;
}
bool Tcell_pure_fluid::update_state(double dt)
{
	//Calculate the work by pressure and viscosity force
	double de_p=0,de_v=0,de_g=0;	
	for (int i=0;i<8;i++)
	{
		de_p=de_p+_node_ptr[i]->_velocity_temp*_corner_stress_force[i];
		de_v=de_v+_node_ptr[i]->_velocity_temp*_corner_viscosity_force[i];
		de_g=de_g+_node_ptr[i]->_velocity_temp*_corner_hourglass_force[i];
	}
	if (de_v>1e-5)
	{
		//cout<<"Error in Aritficitial force work!!!"<<endl;
		//cout<<de_v<<endl;
	}
	//Calculate internal energy before cell volume is updated
	double total_internal_energy=_gausspoint[0]->G_density()*_gausspoint[0]->G_internal_energy()*_cell_volume;
	//Update geometry variables
	bool negative;
	calculate_cell_volume("Full",negative);
	//Update density after cell volume is updated
	_gausspoint[0]->S_density(_material_mass/_cell_volume);
	//Update energy by compatible scheme	
	total_internal_energy=total_internal_energy-(de_p+de_v+de_g)*dt;
	bool failed;
	_gausspoint[0]->update_state(_material_mass,_cell_volume,total_internal_energy,failed);
	if (failed)
	{
		cout<<"Error:soundspeed in cell "<<_id<<" is negative"<<endl;
		system("Pause");
		exit(0);
	}

	//The corner force must be set to zero when the energy is calculated
	for (int i=0;i<8;i++)
	{
		_corner_stress_force[i].value(0,0,0);
		_corner_viscosity_force[i].value(0,0,0);
		_corner_hourglass_force[i].value(0,0,0);
	}
	return negative;
}
void Tcell_pure_fluid::reset_state_after_remapping()
{
	bool failed;
	_gausspoint[0]->update_state(_material_mass,_cell_volume,_material_internal_energy,failed);
	assemble_nodal_mass_moment();
	calculate_corner_mass();
	return;
}
void Tcell_pure_fluid::calculate_average_variable(double  &av_density,double &av_pressure,double &av_soundspeed)
{
	av_density=_gausspoint[0]->G_density();
	av_pressure=_gausspoint[0]->G_pressure();
	av_soundspeed=_gausspoint[0]->G_soundspeed();
	return;
}
void Tcell_pure_fluid::assemble_corner_inertance()
{
	double mass=_cell_volume*_gausspoint[0]->G_density();
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->accumulate_inertance(mass/8.0);
	}
	return;
}
void Tcell_pure_fluid::calculate_gradient_in_cell()
{
	double k[3][3];
	vec3D b,c;
	bool degenerate;
	double det;
	for (int j=0;j<3;j++)
	{
		for (int l=0;l<3;l++)
		{
			k[j][l]=0;
		}
	}
	vec3D cell_centroid=_cell_polyhedron.G_centroid();
	for (unsigned int j=0;j<_around.size();j++)
	{
		double ro,roe,ro0,roe0;
		ro=_around[j]->G_gausspoint(0)->G_density();ro0=_gausspoint[0]->G_density();
		roe=_around[j]->G_gausspoint(0)->G_roe();roe0=_gausspoint[0]->G_roe();
		vec3D centroid_around=_around[j]->_cell_polyhedron.G_centroid();
		b=b+(centroid_around-cell_centroid)*(ro-ro0);
		c=c+(centroid_around-cell_centroid)*(roe-roe0);
		for (int s = 0; s < 3; s++)
		{
			for (int t = s; t < 3; t++)
			{
				k[s][t]=k[s][t]+(centroid_around-cell_centroid).access(s+1)*(centroid_around-cell_centroid).access(t+1);
			}
		}
	}
	k[1][0]=k[0][1];k[2][0]=k[0][2];k[2][1]=k[1][2];
	inverse(k,det,degenerate);
	if (degenerate)
	{
		//If k is degenerated, set the gradient to be zero
		_g_ro.x=_g_ro.y=_g_ro.z=0;
		_g_roe.x=_g_roe.y=_g_roe.z=0;
	}
	else
	{
		_g_ro=b.multiply_by_matrix(k);_g_roe=c.multiply_by_matrix(k);
	}
	return;
}
void Tcell_pure_fluid::limit_gradient_in_cell()
{
	double ro_min,ro_max,roe_min,roe_max;           //The minimum and maximum density and energy in surrounding cells
	double ro_cell,roe_cell;
	double ro_around,roe_around;                                                                
	double ro_vertex,roe_vertex;                     
	double limiter_ro,limiter_roe;                  //limiter factor at each vertexes
	double limiter_ro_min,limiter_roe_min;          //The minimum limiter factor which is used in gradient limitation
	vec3D vertex;

	vec3D cell_centroid=_cell_polyhedron.G_centroid();
	ro_max=ro_min=_gausspoint[0]->G_density();roe_max=roe_min=_gausspoint[0]->G_roe();
	//Calculate the minimum and maximum density and energy in surrounding cells
	for (unsigned int j=0;j<_around.size();j++)
	{
		ro_around=_around[j]->G_gausspoint(0)->G_density();
		roe_around=_around[j]->G_gausspoint(0)->G_roe();
		if (ro_around>ro_max) ro_max=ro_around;
		if (ro_around<ro_min) ro_min=ro_around;
		if (roe_around>roe_max) roe_max=roe_around;
		if (roe_around<roe_min) roe_min=roe_around;
	}
	//Calculate the minimum limiter factor
	limiter_ro_min=limiter_roe_min=1e50;		
	ro_cell=_gausspoint[0]->G_density();
	roe_cell=_gausspoint[0]->G_roe();
	for (int i=0;i<_cell_polyhedron.G_numv();i++)
	{
		vertex=_cell_polyhedron.G_vertex(i);

		//Calculate the limiter for ro
		ro_vertex=ro_cell+_g_ro*(vertex-cell_centroid);
		if (ro_vertex>ro_cell)
		{
			limiter_ro=minval(1.0,(ro_max-ro_cell)/(ro_vertex-ro_cell));
		}
		else if (ro_vertex<ro_cell)
		{
			limiter_ro=minval(1.0,(ro_min-ro_cell)/(ro_vertex-ro_cell));
		}
		else
		{
			limiter_ro=1;
		}
		if (limiter_ro<limiter_ro_min) limiter_ro_min=limiter_ro;

		//Calculate the limiter for roe
		roe_vertex=roe_cell+_g_roe*(vertex-cell_centroid);
		if (roe_vertex>roe_cell)
		{
			limiter_roe=minval(1.0,(roe_max-roe_cell)/(roe_vertex-roe_cell));
		}
		else if (roe_vertex<roe_cell)
		{
			limiter_roe=minval(1.0,(roe_min-roe_cell)/(roe_vertex-roe_cell));
		}
		else
		{
			limiter_roe=1;
		}
		if (limiter_roe<limiter_roe_min) limiter_roe_min=limiter_roe;
	}
	_g_ro=_g_ro*limiter_ro_min;
	_g_roe=_g_roe*limiter_roe_min;
	return;
}
double Tcell_pure_fluid::calculate_cell_internal_energy()
{
	return _material_mass*_gausspoint[0]->G_internal_energy();
}
void Tcell_pure_fluid::clear_cell_variable()
{
	Tcell_fluid_base::clear_cell_variable();
	_material_mass=0;
	_material_internal_energy=0;
	return;
}
void Tcell_pure_fluid::accumulate_variables(int mid,double volume,vec3D int_x,double mass,double internal_energy)
{
	_material_mass=_material_mass+mass;
	_material_internal_energy=_material_internal_energy+internal_energy;
	return;
}