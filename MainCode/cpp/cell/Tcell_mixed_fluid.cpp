#include "Tcell_mixed_fluid.h"
Tcell_mixed_fluid::Tcell_mixed_fluid(int cell_id,int num_material,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat)
	:Tcell_fluid_base(cell_id,num_material,node_ptr,mat)
{
	_material_mass.resize(num_material);_material_internal_energy.resize(num_material);
	_material_centroid.resize(num_material);_material_fraction.resize(num_material);
	_material_centroid_std.resize(num_material);
	_fraction_increment.resize(num_material);
	_g_ro.resize(num_material);_g_roe.resize(num_material);
	_material_polyhedron.resize(num_material);
	for (int i = 0; i < num_material; i++)
	{
		_material_internal_energy[i]=_material_fraction[i]=_material_internal_energy[i]=_fraction_increment[i]=0;
		_material_polyhedron[i]=NULL;
	}
	_average_pressure=0;
}
Tcell_mixed_fluid::Tcell_mixed_fluid(int cell_id,int num_material,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat[])
	:Tcell_fluid_base(cell_id,num_material,node_ptr,mat)
{
	_material_mass.resize(num_material);_material_internal_energy.resize(num_material);
	_material_centroid.resize(num_material);_material_fraction.resize(num_material);
	_material_centroid_std.resize(num_material);
	_fraction_increment.resize(num_material);
	_g_ro.resize(num_material);_g_roe.resize(num_material);
	_material_polyhedron.resize(num_material);
	for (int i = 0; i < num_material; i++)
	{
		_material_internal_energy[i]=_material_fraction[i]=_material_internal_energy[i]=_fraction_increment[i]=0;
		_material_polyhedron[i]=NULL;
	}
	_average_pressure=0;
}
void Tcell_mixed_fluid::clear_cell_variable()
{
	Tcell_fluid_base::clear_cell_variable();
	int num_material=_gausspoint.size();
	for (int i = 0; i < num_material; i++)
	{
		_material_mass[i]=_material_internal_energy[i]=0;
		_material_fraction[i]=_fraction_increment[i]=0;
		_material_centroid[i].value(0,0,0);
		_material_centroid_std[i].value(0,0,0);
	}
	delete_material_polyhedron();
	_average_pressure=0;
	return;
}
void Tcell_mixed_fluid::reset_geometry_initial_condition()
{
	int num_material=_gausspoint.size();
	double fraction_immersed=0;
	vec3D centroid_immersed;
	for (int i = 1; i < num_material; i++)
	{
		if (_material_fraction[i]<1e-10) _material_fraction[i]=0;
		if (abs(_material_fraction[i]-1)<1e-10) _material_fraction[i]=1;
	}
	for (int i = 1; i < num_material; i++)
	{
		fraction_immersed=fraction_immersed+_material_fraction[i];
		centroid_immersed=centroid_immersed+_material_centroid[i]*_material_fraction[i];
	}
	_material_fraction[0]=1-fraction_immersed;
	if (_material_fraction[0]<1e-10)
	{
		_material_fraction[0]=0;
	}
	else
	{
		_material_centroid[0]=(_cell_polyhedron.G_centroid()-centroid_immersed)/_material_fraction[0];
	}
	for (int i = 0; i < num_material; i++)
	{
		if (_material_fraction[i]>0 && _material_fraction[i]<1)
		{
			_material_centroid_std[i]=calculate_std_coordinate(_material_centroid[i]);
		}
	}
	return;
}
void Tcell_mixed_fluid::reset_conservative_variable_initial_condition()
{
	int num_material=_gausspoint.size();
	calculate_corner_mass();
	for (int i = 0; i < num_material; i++)
	{
		_material_mass[i]=_material_fraction[i]*_cell_volume*_gausspoint[i]->G_density();
		_material_internal_energy[i]=_material_mass[i]*_gausspoint[i]->G_internal_energy();
	}
	return;
}
void Tcell_mixed_fluid::reset_material_initial_state(int mid)
{
	if (_material_fraction[mid]>1e-10)
	{
		_material_centroid[mid]=_material_centroid[mid]/_material_fraction[mid];
		_material_fraction[mid]=_material_fraction[mid]/_cell_volume;
	}
	else
	{
		//Ignore the material with very small volume fraction
		_material_fraction[mid]=0;
	}
	return;
}
void Tcell_mixed_fluid::calculate_corner_mass()
{
	double average_density=0;
	int num_material=_gausspoint.size();
	for (int i = 0; i < num_material; i++)
	{
		average_density=average_density+_material_fraction[i]*_gausspoint[i]->G_density();
	}
	for (int i = 0; i < 8; i++)
	{
		_corner_mass[i]=_corner_volume[i]*average_density;
	}
	return;
}
void Tcell_mixed_fluid::assemble_corner_inertance()
{
	int num_material=_gausspoint.size();
	double mass=0;
	for (int i = 0; i < num_material; i++)
	{
		mass=mass+_cell_volume*_gausspoint[i]->G_density()*_material_fraction[i];
	}
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->accumulate_inertance(mass/8.0);
	}
	return;
}
bool Tcell_mixed_fluid::calculate_corner_stress_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3],double dt)
{
	double temp_pressure;        //The pressure at predict step
	bool negative;
	bool is_pure=false;
	double temp_volume=calculate_cell_volume("Temporary",negative);
	int num_material=_gausspoint.size();
	for (int i = 0; i < num_material; i++)
	{
		if (_material_fraction[i]==1)
		{
			//This cell is a pure cell
			//if (i!=2)
			//{
				double ro,e;
				ro=_gausspoint[i]->G_density();e=_gausspoint[i]->G_internal_energy();
				temp_pressure=_gausspoint[i]->G_EOS()->predict_p(ro,e,temp_volume,_cell_volume);
				_average_pressure=temp_pressure;
				is_pure=true;
			//}
			//else
			//{
			//	_average_pressure=temp_pressure=_gausspoint[i]->G_pressure();
			//	is_pure=true;
			//}
		}
	}
	if (!is_pure)
	{
		//This cell is a mixed cell,using pressure relxison model
		//temp_pressure=0;
		//for (int i = 0; i < num_material; i++)
		//{
		//	temp_pressure=temp_pressure+_material_fraction[i]*_gausspoint[i]->G_pressure();
		//	_fraction_increment[i]=0;
		//}
		//_average_pressure=temp_pressure;
		double B[10],avp,avB;
		double ro,soundspeed,pressure,e;
		double temp1,temp2;
		temp1=temp2=0;
		double p[2];
		for (int i=0;i<num_material;i++)
		{
			if (_material_fraction[i]>0)
			{
				ro=_gausspoint[i]->G_density();soundspeed=_gausspoint[i]->G_soundspeed();e=_gausspoint[i]->G_internal_energy();
				B[i]=ro*soundspeed*soundspeed+ro*soundspeed*_cell_length/dt;
				pressure=_gausspoint[i]->calculate_pressure();
				p[i]=_gausspoint[i]->G_EOS()->predict_p(ro,e,temp_volume,_cell_volume);
				temp1=temp1+_material_fraction[i]*pressure/B[i];
				temp2=temp2+_material_fraction[i]/B[i];
			}	
		}
		avp=temp1/temp2;avB=1.0/temp2;
		temp_pressure=avp-avB*(temp_volume-_cell_volume)/_cell_volume;
		_average_pressure=temp_pressure;
		//Update the material fraction
		double trial_fraction_increment[10];      //The trial fraction increment without limit
		double fraction_new[10];                  //The material fraction in predict step
		double ratio=1e10;
		//Calulate original fraction increments and limite them
		for (int i=0;i<num_material;i++)
		{
			if (_material_fraction[i]>0)
			{
				pressure=_gausspoint[i]->calculate_pressure();
				trial_fraction_increment[i]=_material_fraction[i]*((pressure-avp)/B[i]+(avB/B[i]-1)*(temp_volume-_cell_volume)/_cell_volume);
				trial_fraction_increment[i]=sign(trial_fraction_increment[i])*minval(abs(trial_fraction_increment[i]),ratio*_material_fraction[i]);
			}
			else
			{
				trial_fraction_increment[i]=0;
			}
		}
		//Renormalize the fraction increment
		double df_positive,df_negative;
		df_positive=df_negative=0;
		for (int i = 0; i < num_material; i++)
		{
			if (trial_fraction_increment[i]>0)
			{
				df_positive=df_positive+trial_fraction_increment[i];
			}
			else if (trial_fraction_increment[i]<0)
			{
				df_negative=df_negative+abs(trial_fraction_increment[i]);
			}
		}
		if (df_negative<=df_positive)
		{
			for (int i = 0; i < num_material; i++)
			{
				if (trial_fraction_increment[i]>0)
				{
					trial_fraction_increment[i]=trial_fraction_increment[i]*df_negative/df_positive;
				}
			}
		}
		else if (df_negative>df_positive)
		{
			for (int i = 0; i < num_material; i++)
			{
				if (trial_fraction_increment[i]<0)
				{
					trial_fraction_increment[i]=trial_fraction_increment[i]*df_positive/df_negative;
				}
			}
		}
		//Set the fraction increment in the target cell
		for (int i = 0; i < num_material; i++)
		{
			_fraction_increment[i]=2*trial_fraction_increment[i];
		}
		//Exam of the fraction increment
		double total_increment=0;
		for (int i = 0; i < num_material; i++)
		{
			total_increment=total_increment+_fraction_increment[i];
		}
		if (abs(total_increment)>1e-10)
		{
			cout<<"The total increment of fraction is not equal to zero!!!"<<endl;
			system("Pause");
			exit(0);
		}
		for (int i = 0; i < num_material; i++)
		{
			fraction_new[i]=_material_fraction[i]+_fraction_increment[i];
			if (fraction_new[i]<0 || fraction_new[i]>1)
			{
				cout<<"The volume fraction excess!!!"<<endl;
				system("Pause");
				exit(0);
			}
		}
	}
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<3;j++)
		{
			_corner_stress_force[i]=_corner_stress_force[i]+surfacial[i][j]*temp_pressure;
		}
	}
	return negative;
}
void Tcell_mixed_fluid::calculate_corner_hourglass_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3])
{
	int num_material=_gausspoint.size();
	double dp[8];
	for (int i = 0; i < 8; i++)
	{
		dp[i]=0;
	}
	for (int mid = 0; mid < num_material; mid++)
	{
		//The hourglass force only applied on the pure cells
		if (_material_fraction[mid]==1)
		{
				double ro=_gausspoint[mid]->G_density();
				double p=_gausspoint[mid]->G_pressure();
				double e=_gausspoint[mid]->G_internal_energy();
				double c=_gausspoint[mid]->G_soundspeed();
				for (int i=0;i<8;i++)
				{
					dp[i]=_gausspoint[mid]->G_EOS()->ro_e_p(_corner_mass[i]/_corner_volume[i],e)-p;
				}
		}
	}
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<3;j++)
		{
 			_corner_hourglass_force[i]=_corner_hourglass_force[i]+(surfacial[i][j]*dp[i]
 			+assistant[i][j]*(dp[i]-dp[Tcell_fluid_base::cell_corner_v[i][j]])*0.5);
		}
	}
	return;
}
bool Tcell_mixed_fluid::update_state(double dt)
{
	//Calculate the work by pressure and viscosity force
	double de_p=0,de_v=0,de_g=0;	
	for (int i=0;i<8;i++)
	{
		de_p=de_p+_node_ptr[i]->_velocity_temp*_corner_stress_force[i];
		de_v=de_v+_node_ptr[i]->_velocity_temp*_corner_viscosity_force[i];
		de_g=de_g+_node_ptr[i]->_velocity_temp*_corner_hourglass_force[i];
	}
	int num_material=_gausspoint.size();
	int num_different_material;
	int id[10];
	//Identify the materials in cell
	num_different_material=classify_mixed_cell(id);
	bool negative,failed;
	if (num_different_material==1)
	{
		//For pure cells
		for (int i = 0; i < num_material; i++)
		{
			if (_material_fraction[i]==1)
			{
				//Calculate internal energy before cell volume is updated
				double total_internal_energy=_gausspoint[i]->G_density()*_gausspoint[i]->G_internal_energy()*_cell_volume;
				//Update geometry variables
				calculate_cell_volume("Full",negative);
				//Update density after cell volume is updated
				_gausspoint[i]->S_density(_material_mass[i]/_cell_volume);
				//Update energy by compatible scheme	
				total_internal_energy=total_internal_energy-(de_p+de_v+de_g)*dt;
				_gausspoint[i]->update_state(_material_mass[i],_cell_volume,total_internal_energy,failed);
				if (failed)
				{
					cout<<"Error:soundspeed in cell "<<_id<<",material "<<i<<" is negative"<<endl;
					system("Pause");
					exit(0);
				}
				break;
			}
		}
	}
	else
	{
		//For multi-material cells
		//Calculate internal energy before cell volume and fraction is updated
		double internal_energy[10];
		double volume_pre;
		volume_pre=_cell_volume;       //The volume at the previous time step
		for (int i=0;i<num_different_material;i++)
		{
			internal_energy[id[i]]=_gausspoint[id[i]]->G_roe()*_cell_volume*_material_fraction[id[i]];
		}
		//Update geometry variables
		calculate_cell_volume("Full",negative);
		//Update volume fraction
		for (int i=0;i<num_different_material;i++)
		{
			_material_fraction[id[i]]=_material_fraction[id[i]]+_fraction_increment[id[i]];
		}
		//Update the state of each material
		double total_mass=0;
		for (int i = 0; i < num_different_material; i++)
		{
			total_mass=total_mass+_material_mass[id[i]];
		}
		int mid;
		for (int j=0;j<num_different_material;j++)
		{
			mid=id[j];
			internal_energy[mid]=internal_energy[mid]-(_material_fraction[mid]*de_p+_material_mass[mid]*(de_v+de_g)/total_mass)*dt-_fraction_increment[mid]*volume_pre*_average_pressure;
			//internal_energy[mid]=internal_energy[mid]-(_material_fraction[mid]*de_p+_material_fraction[mid]*(de_v+de_g))*dt-_fraction_increment[mid]*volume_pre*_average_pressure;
			_gausspoint[mid]->update_state(_material_mass[mid],_cell_volume*_material_fraction[mid],internal_energy[mid],failed);
			if (failed)
			{
				cout<<"Error:soundspeed in cell "<<_id<<",material "<<mid<<" is negative"<<endl;
				system("Pause");
				exit(0);
			}
		}
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
void Tcell_mixed_fluid::accumulate_variables(int mid,double volume,vec3D int_x,double mass,double internal_energy)
{
	_material_mass[mid]=_material_mass[mid]+mass;
	_material_internal_energy[mid]=_material_internal_energy[mid]+internal_energy;
	_material_fraction[mid]=_material_fraction[mid]+volume;
	_material_centroid[mid]=_material_centroid[mid]+int_x;
	return;
}
void Tcell_mixed_fluid::reset_state_after_remapping()
{
	int num_material=_gausspoint.size();
	//Update the volume fraction first
	for (int i = 0; i < num_material; i++)
	{
		_material_fraction[i]=_material_fraction[i]/_cell_volume;
	}
	//Filter the fraction
	filter_fraction();
	//Update physical variables
	for (int i = 0; i < num_material; i++)
	{
		if (_material_fraction[i]>0)
		{
			bool failed;
			_gausspoint[i]->update_state(_material_mass[i],_cell_volume*_material_fraction[i],_material_internal_energy[i],failed);
			_material_centroid[i]=_material_centroid[i]/(_cell_volume*_material_fraction[i]);
			_material_centroid_std[i]=calculate_std_coordinate(_material_centroid[i]);
		}

	}
	assemble_nodal_mass_moment();
	calculate_corner_mass();
	return;
}
void Tcell_mixed_fluid::reset_state_after_overlapping_exp()
{
	int num_material = _gausspoint.size();
	//Update the volume fraction first
	for (int i = 0; i < num_material; i++)
	{
		_material_fraction[i] = _material_fraction[i] / _cell_volume;
	}
	//Filter the fraction
	filter_fraction();
	//Update physical variables
	for (int i = 0; i < num_material; i++)
	{
		if (_material_fraction[i] > 0)
		{
			bool failed;
			_gausspoint[i]->update_state(_material_mass[i], _cell_volume*_material_fraction[i], _material_internal_energy[i], failed);
			_material_centroid[i] = _material_centroid[i] / (_cell_volume*_material_fraction[i]);
			_material_centroid_std[i] = calculate_std_coordinate(_material_centroid[i]);
		}
	}
	return;
}
void Tcell_mixed_fluid::calculate_average_variable(double  &av_density,double &av_pressure,double &av_soundspeed)
{
	int num_material=_gausspoint.size();
	av_density=av_pressure=av_soundspeed=0;
	for (int i = 0; i < num_material; i++)
	{
		av_density=av_density+_gausspoint[i]->G_density()*_material_fraction[i];
		av_pressure=av_pressure+_gausspoint[i]->G_pressure()*_material_fraction[i];
		av_soundspeed=av_soundspeed+_gausspoint[i]->G_soundspeed()*_material_fraction[i];
	}
	return;
}
double Tcell_mixed_fluid::calculate_cell_internal_energy()
{
	int num_material=_gausspoint.size();
	double internal_energy=0;
	for (int i = 0; i < num_material; i++)
	{
		internal_energy=internal_energy+_material_mass[i]*_gausspoint[i]->G_internal_energy();
	}
	return internal_energy;
}
void Tcell_mixed_fluid::calculate_gradient_in_cell()
{
	if (_id==114)
	{
		double a=0;
	}
	double k[3][3];
	vec3D b,c;
	bool degenerate;
	double det;
	int num_material=_gausspoint.size();
	for (int mid = 0; mid < num_material; mid++)
	{
		if (_material_fraction[mid]>0)
		{
			b.value(0,0,0);c.value(0,0,0);
			for (int j=0;j<3;j++)
			{
				for (int l=0;l<3;l++)
				{
					k[j][l]=0;
				}
			}
			vec3D material_centroid=_material_centroid[mid];
			double ro0=_gausspoint[mid]->G_density();
			double roe0=_gausspoint[mid]->G_roe();
			vec3D g_ref;
			g_ref.value(5,6,7);g_ref=g_ref*(mid+1);
			double error0=ro0-material_centroid*g_ref-1;
			for (unsigned int j = 0; j < _around.size(); j++)
			{
				if (_around[j]->G_material_fraction(mid)>0)
				{
					double ro,roe;
					ro=_around[j]->G_gausspoint(mid)->G_density();		
					roe=_around[j]->G_gausspoint(mid)->G_roe();
					vec3D centroid_around=_around[j]->G_material_centroid(mid);
					double error1=ro-centroid_around*g_ref-1;
					b=b+(centroid_around-material_centroid)*(ro-ro0);
					c=c+(centroid_around-material_centroid)*(roe-roe0);
					for (int s = 0; s < 3; s++)
					{
						for (int t = s; t < 3; t++)
						{
							k[s][t]=k[s][t]+(centroid_around-material_centroid).access(s+1)*(centroid_around-material_centroid).access(t+1);
						}
					}	
				}	
			}
			k[1][0]=k[0][1];k[2][0]=k[0][2];k[2][1]=k[1][2];
			inverse(k,det,degenerate);
			if (degenerate)
			{
				//If k is degenerated, set the gradient to be zero
				_g_ro[mid].value(0,0,0);_g_roe[mid].value(0,0,0);
			}
			else
			{
				_g_ro[mid]=b.multiply_by_matrix(k);_g_roe[mid]=c.multiply_by_matrix(k);
			}
		}
	}
	return;
}
void Tcell_mixed_fluid::limit_gradient_in_cell()
{
	double ro_min,ro_max,roe_min,roe_max;           //The minimum and maximum density and energy in surrounding cells
	double ro_material,roe_material;
	double ro_around,roe_around;                                                                
	double ro_vertex,roe_vertex;                     
	double limiter_ro,limiter_roe;                  //limiter factor at each vertexes
	double limiter_ro_min,limiter_roe_min;          //The minimum limiter factor which is used in gradient limitation
	vec3D vertex;
	int num_material=_gausspoint.size();
	for (int mid = 0; mid < num_material; mid++)
	{
		if (_material_fraction[mid]>0)
		{
			vec3D material_centroid=_material_centroid[mid];
			ro_max=ro_min=_gausspoint[mid]->G_density();roe_max=roe_min=_gausspoint[mid]->G_roe();
			for (unsigned int j=0;j<_around.size();j++)
			{
				if (_around[j]->G_material_fraction(mid)>0)
				{
					ro_around=_around[j]->G_gausspoint(mid)->G_density();
					roe_around=_around[j]->G_gausspoint(mid)->G_roe();
					if (ro_around>ro_max) ro_max=ro_around;
					if (ro_around<ro_min) ro_min=ro_around;
					if (roe_around>roe_max) roe_max=roe_around;
					if (roe_around<roe_min) roe_min=roe_around;
				}
			}
			limiter_ro_min=limiter_roe_min=1e50;		
			ro_material=_gausspoint[mid]->G_density();
			roe_material=_gausspoint[mid]->G_roe();
			for (int i=0;i<_material_polyhedron[mid]->G_numv();i++)
			{
				vertex=_material_polyhedron[mid]->G_vertex(i);
				//Calculate the limiter for ro
				ro_vertex=ro_material+_g_ro[mid]*(vertex-material_centroid);
				if (ro_vertex>ro_material)
				{
					limiter_ro=minval(1.0,(ro_max-ro_material)/(ro_vertex-ro_material));
				}
				else if (ro_vertex<ro_material)
				{
					limiter_ro=minval(1.0,(ro_min-ro_material)/(ro_vertex-ro_material));
				}
				else
				{
					limiter_ro=1;
				}
				if (limiter_ro<limiter_ro_min) limiter_ro_min=limiter_ro;

				//Calculate the limiter for roe
				roe_vertex=roe_material+_g_roe[mid]*(vertex-material_centroid);
				if (roe_vertex>roe_material)
				{
					limiter_roe=minval(1.0,(roe_max-roe_material)/(roe_vertex-roe_material));
				}
				else if (roe_vertex<roe_material)
				{
					limiter_roe=minval(1.0,(roe_min-roe_material)/(roe_vertex-roe_material));
				}
				else
				{
					limiter_roe=1;
				}
				if (limiter_roe<limiter_roe_min) limiter_roe_min=limiter_roe;
			}
			_g_ro[mid]=_g_ro[mid]*limiter_ro_min;
			_g_roe[mid]=_g_roe[mid]*limiter_roe_min;
		}
	}
	return;
}
int Tcell_mixed_fluid::classify_mixed_cell(int (&id)[10])
{
	int num_different_material=0;
	int num_material=_gausspoint.size();
	for (int i = 0; i < num_material; i++)
	{
		id[i]=0;
	}
	for (int i = 0; i < num_material; i++)
	{
		if (_material_fraction[i]>0)
		{
			id[num_different_material]=i;
			num_different_material=num_different_material+1;
		}
	}
	return num_different_material;
}
void Tcell_mixed_fluid::filter_fraction()
{
	int num_material=_gausspoint.size();
	for (int i = 0; i < num_material; i++)
	{
		if (_material_fraction[i]<1e-4)
		{
			_material_fraction[i]=0;
		}
		if (abs(_material_fraction[i]-1)<1e-4)
		{
			_material_fraction[i]=1;
		}
	}
	return;
}
int Tcell_mixed_fluid::is_surface_cell(int mid)
{
	int num_surface_piece=0;
	for (int i=0;i<6;i++)
	{
		if (_adjacent[i]!=NULL)
		{
			if (_adjacent[i]->G_material_fraction(mid)==0)
			{
				//If this cell is adjacent to a cell which does not contain material mid
				//This face is a material surface of material mid
				num_surface_piece=num_surface_piece+4;
			}
		}
	}
	return num_surface_piece;
}
void Tcell_mixed_fluid::add_surface_cell_surface(int mid,int** piece_list,int vertex_id_start,int &piece_id)
{
	//Define the id of center vertex of each cell face in cell polyhedron
	int id[6];
	id[0]=12;id[5]=13;
	for (int i=1;i<5;i++)
	{
		id[i]=i+7;
	}
	for (int i=0;i<6;i++)
	{
		if (_adjacent[i]!=NULL)
		{
			if (_adjacent[i]->G_material_fraction(mid)==0)
			{
				for (int j=0;j<4;j++)
				{
					piece_list[piece_id][0]=vertex_id_start+id[i];
					piece_list[piece_id][1]=vertex_id_start+Tcell_fluid_base::cell_face[i][j];
					piece_list[piece_id][2]=vertex_id_start+Tcell_fluid_base::cell_face[i][(j+1)%4];
					piece_id=piece_id+1;
				}
			}
		}
	}
	return;
}
void Tcell_mixed_fluid::delete_material_polyhedron()
{
	int num_material=_gausspoint.size();
	for (int i = 0; i < num_material; i++)
	{
		if (_material_polyhedron[i]!=NULL)
		{
			delete _material_polyhedron[i];
			_material_polyhedron[i]=NULL;
		}
	}
	return;
}
void Tcell_mixed_fluid::update_material_centroid()
{
	int num_material=_gausspoint.size();
	for (int i = 0; i < num_material; i++)
	{
		if (_material_fraction[i]>0 && _material_fraction[i]<1)
		{
			_material_centroid[i]=shape_function_interpolate(_material_centroid_std[i],"Position");
		}
	}
	return;
}
void Tcell_mixed_fluid::S_material_base(TMAT_base* mat[])
{
	int num_material=_gausspoint.size();
	for (int i = 0; i < num_material; i++)
	{
		_gausspoint[i]=mat[i];
	}
	return;
}
void Tcell_mixed_fluid::output_cell_monitor(ofstream& output,double t)
{
	output<<t<<" "<<_material_fraction[0]*_gausspoint[0]->G_pressure()+_material_fraction[1]*_gausspoint[1]->G_pressure()<<endl;
	return;
}
void Tcell_mixed_fluid::set_pressure(double t)
{
	double p;
	if (t<=0.04)
	{
		p=13;
	}
	else if (t<0.04 && t<=0.125)
	{
		p=13-12.5*(t-0.04)/(0.125-0.04);
	}
	else
	{
		p=0.5;
	}
	_gausspoint[2]->S_pressure(p);
}