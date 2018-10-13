#include "Tbody_MPM.h"
#include <iomanip>
void Tbody_MPM::input_body(ifstream& input)
{
	Skeyword keyword;
	read_in_keyword_file(input,keyword);
	//Read the body global information from keyword
	_nume=keyword.cell_8_list.size();
	_nump=keyword.node_list.size();
	_endtime=keyword.time_control.endtime;
	_CFL=keyword.time_control.CFL;
	//----------------------------------------------------
	//Read particle group list
	//----------------------------------------------------
	_particle_group_list=new vector<TMPM_particle>[_num_part];
	int num_particle;
	num_particle=keyword.particle_list.size();        //The total number of particles
	vec3D x_min_particle(1e10,1e10,1e10),x_max_particle(-1e10,-1e10,-1e10);
	for (int i = 0; i < num_particle; i++)
	{
		Sparticle_temp this_particle=keyword.particle_list[i];
		int group_id=this_particle.group_id;
		int part_id=this_particle.part_id;
		if (group_id>_num_part-1)
		{
			cout<<"Error: The particle group id excess!"<<endl;
			system("Pause");
			exit(0);
		}
		int id_in_part=_particle_group_list[group_id].size();
		TMPM_particle new_particle(this_particle.position,this_particle.volume,group_id,part_id,id_in_part);
		x_min_particle.x=minval(x_min_particle.x,this_particle.position.x);x_min_particle.y=minval(x_min_particle.y,this_particle.position.y);
		x_min_particle.z=minval(x_min_particle.z,this_particle.position.z);x_max_particle.x=maxval(x_max_particle.x,this_particle.position.x);
		x_max_particle.y=maxval(x_max_particle.y,this_particle.position.y);x_max_particle.z=maxval(x_max_particle.z,this_particle.position.z);
		_particle_group_list[part_id].push_back(new_particle);
	}
	//-----------------------------------------------------
	//Read material porperty for particle
	//-----------------------------------------------------
	int material_id,EOS_id;
	TMPM_particle* particle_ptr;
	for (int i = 0; i < _num_part; i++)
	{
		material_id=keyword.part_list[i].material_id;
		EOS_id=keyword.part_list[i].EOS_id;
		Smaterial this_mat=keyword.material_list[material_id-1];
		SEOS this_EOS;
		if (keyword.EOS_list.size()>0)
		{
			this_EOS=keyword.EOS_list[EOS_id-1];
		}
		int num_particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < num_particles_in_part; j++)
		{
			TMAT_base* mat_ptr;
			particle_ptr=&_particle_group_list[i][j];
			mat_ptr=generate_material(this_mat,this_EOS,0);
			particle_ptr->_gausspoint=mat_ptr;
			particle_ptr->_mass=particle_ptr->_volume*particle_ptr->_gausspoint->G_density();
		}
	}
	//-------------------------------------------------------
	//Generate the background
	//-------------------------------------------------------
	vec3D x_min=keyword.background.x_min;
	vec3D x_max=keyword.background.x_max;
	int nx=keyword.background.nx;
	int ny=keyword.background.ny;
	int nz=keyword.background.nz;
	for (int i = 1; i <= 3; i++)
	{
		if (x_min_particle.access(i)<x_min.access(i) || x_max_particle.access(i)>x_max.access(i))
		{
			cout<<"Error: The initial particle model is out of the background grid!"<<endl;
			system("Pause");
			exit(0);
		}
	}
	_background.generate_background(x_min,x_max,nx,ny,nz);
	_nump=_background._total_node;
	_nume=_background._total_cell;
	//----------------------------------------------------------------
	//Read initial velocity
	//----------------------------------------------------------------
	int num_initial_velocity=keyword.initial_velocity_list.size();
	for (int i = 0; i < num_initial_velocity; i++)
	{
		int part_id=keyword.initial_velocity_list[i].part_id;
		int particles_in_part=_particle_group_list[part_id].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			_particle_group_list[i][j]._velocity=keyword.initial_velocity_list[i].velocity;
		}
	}
	//-------------------------------------------------------
	//Read boundary condition
	//-------------------------------------------------------
	int num_boundary_condition=keyword.MPM_boundary_condition_list.size();
	for (int n = 0; n < num_boundary_condition; n++)
	{
		int i_min=keyword.MPM_boundary_condition_list[n].strat_subscript[0];
		int j_min=keyword.MPM_boundary_condition_list[n].strat_subscript[1];
		int k_min=keyword.MPM_boundary_condition_list[n].strat_subscript[2];
		int i_max=keyword.MPM_boundary_condition_list[n].terminate_subscript[0];
		int j_max=keyword.MPM_boundary_condition_list[n].terminate_subscript[1];
		int k_max=keyword.MPM_boundary_condition_list[n].terminate_subscript[2];
		int index;
		for (int i = i_min; i < i_max+1; i++)
		{
			for (int j = j_min; j < j_max+1; j++)
			{
				for (int k = k_min; k < k_max+1; k++)
				{
					index=(ny+1)*(nz+1)*i+(nz+1)*j+k;
					for (int s = 0; s < 3; s++)
					{
						_background._node_list[index]._bc_type_position[s]=keyword.MPM_boundary_condition_list[n].constraint[s];
					}
				}
			}
		}
	}
	//---------------------------------------------------------------------------------
	//Read external load on background
	//---------------------------------------------------------------------------------
	int num_background_load=keyword.MPM_background_load_list.size();
	for (int n = 0; n < num_background_load; n++)
	{
		int i_min=keyword.MPM_background_load_list[n].strat_subscript[0];
		int j_min=keyword.MPM_background_load_list[n].strat_subscript[1];
		int k_min=keyword.MPM_background_load_list[n].strat_subscript[2];
		int i_max=keyword.MPM_background_load_list[n].terminate_subscript[0];
		int j_max=keyword.MPM_background_load_list[n].terminate_subscript[1];
		int k_max=keyword.MPM_background_load_list[n].terminate_subscript[2];
		int index;
		for (int i = i_min; i < i_max+1; i++)
		{
			for (int j = j_min; j < j_max+1; j++)
			{
				for (int k = k_min; k < k_max+1; k++)
				{
					index=(ny+1)*(nz+1)*i+(nz+1)*j+k;
					for (int s = 0; s < 3; s++)
					{
						_background._node_list[index]._external_load=keyword.MPM_background_load_list[n].load;
					}
				}
			}
		}
	}
	//-------------------------------------------------------------------------------
	//Read external load on particle
	//--------------------------------------------------------------------------------
	int num_particle_load=keyword.MPM_particle_load_list.size();
	int part_id,particle_id;
	vec3D load;
	for (int i = 0; i < num_particle_load; i++)
	{
		part_id=keyword.MPM_particle_load_list[i].part_id;
		particle_id=keyword.MPM_particle_load_list[i].particle_id;
		load=keyword.MPM_particle_load_list[i].load;
		_particle_group_list[part_id][particle_id]._external_force=load;
	}
	return;
}
void Tbody_MPM::start_up_step()
{
	//Map the particle to background to generate the initial background grid
	map_particle_to_background();
	//Calculate the initial nodal force
	calculate_force_in_body();
	return;
}
bool Tbody_MPM::corrector_step()
{
	//Update the DOF of the background grid
	_background.update_background_DOF(_dt);
	//Update the stress on the particle
	TMPM_particle* particle_ptr;
	for (int i = 0; i < _num_part; i++)
	{
		int particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			//_background.update_particle_state(particle_ptr,_dt);
			_background.update_particle_DOF(particle_ptr,_dt);
		}
	}
	//Create a new background and map particle to the background
	map_particle_to_background();
	for (int i = 0; i < _num_part; i++)
	{
		int particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			_background.update_particle_state(particle_ptr,_dt);
			//_background.update_particle_DOF(particle_ptr,_dt);
		}
	}
	//Calculate the nodal force
	calculate_force_in_body();
	return false;
}
void Tbody_MPM::map_particle_to_background()
{
	_background.initialize();
	TMPM_particle* particle_ptr;
	for (int i = 0; i < _num_part; i++)
	{
		int particle_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particle_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			_background.map_particle_to_cell(particle_ptr);
		}
	}
	_background.reset_nodal_velocity();
	return;
}
double Tbody_MPM::calculate_time_step()
{
	_background.calculate_character_length();
	TMPM_particle* particle_ptr;
	int particles_in_part,cell_id;
	double max_spd,this_spd;
	for (int i = 0; i < _num_part; i++)
	{
		particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			cell_id=particle_ptr->_cell_index;
			max_spd=_background._cell_list[cell_id]._max_spd;
			this_spd=particle_ptr->_gausspoint->G_soundspeed();
			_background._cell_list[cell_id]._max_spd=maxval(max_spd,this_spd);
		}
	}
	double time_step=100000;
	double clength,soundspeed,nodespeed;
	int id=-1;
	for (int i = 0; i < _background._total_cell; i++)
	{
		double cell_clength,cell_soundspeed,cell_nodespeed;
		double cell_timespeed=_background._cell_list[i].calculate_cell_time_step(cell_clength,cell_soundspeed,cell_nodespeed);
		if (time_step>cell_timespeed)
		{
			time_step=cell_timespeed;
			clength=cell_clength;soundspeed=cell_soundspeed;nodespeed=cell_nodespeed;
			id=i;
		}
	}
	//cout<<"clength: "<<clength<<" cps: "<<soundspeed<<" nodespd: "<<nodespeed<<" cell_id: "<<id<<endl;
	return _CFL*time_step;
	return 0.01;
}
void Tbody_MPM::calculate_force_in_body()
{
	calculate_corner_force();
	apply_external_force();
	_background.assemble_nodal_force();
}
void Tbody_MPM::calculate_corner_force()
{
	TMPM_particle* particle_ptr;
	for (int i = 0; i < _num_part; i++)
	{
		int particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			_background.distribute_particle_force(particle_ptr);
		}
	}
	return;
}
void Tbody_MPM::calculate_final_acceleration()
{
	for (int i = 0; i < _background._total_node; i++)
	{
		_background._node_list[i].calculate_acceleration();
	}
	return;
}
void Tbody_MPM::apply_external_force()
{
	_background.apply_background_load();
	return;
}
void Tbody_MPM::output_tecplot(ofstream& output,double ratio)
{
	//_background.plot_background_cell(output,0);
	//_background.plot_background(output);
	TMPM_particle* particle_ptr;
	output<<"TITLE=\"Particle coordinate\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"Z\",\"pressure\""<<endl;
	output<<"Zone T=\""<<_current_time<<"\""<<endl;
	for (int i = 0; i < _num_part; i++)
	{
		int particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			output<<particle_ptr->_position.x<<" "<<particle_ptr->_position.y<<" "<<particle_ptr->_position.z;
			output<<" "<<particle_ptr->_gausspoint->G_seqv()<<endl;
		}
	}
	return;
}
void Tbody_MPM::output_curve()
{
	Tbody_base::output_curve();
	int num_curve=_curve_out_list.size();
	int type,id;
	Tnode* node_ptr;
	TMPM_particle* particle_ptr;
	for (int i = 0; i < num_curve; i++)
	{
		type=_curve_out_list[i]->_type;
		id=_curve_out_list[i]->_id-1;
		_curve_out_list[i]->output.precision(10);
		if (type==0)
		{
			particle_ptr=&_particle_group_list[0][0];
			_curve_out_list[i]->output<<setw(20)<<_current_time<<" "<<particle_ptr->_position.x-0.5<<endl;
			//double stress[6];
			//particle_ptr->_gausspoint->G_stress(stress);
			//for (int j = 0; j < 6; j++)
			//{
			//	_curve_out_list[i]->output<<setw(20)<<stress[j];
			//}
			//_curve_out_list[i]->output<<endl;
			//node_ptr=&_background._node_list[id];
			//vec3D dis=node_ptr->calculate_displacement();
			//_curve_out_list[i]->output<<setw(20)<<_current_time<<" ";
			//_curve_out_list[i]->output<<setw(20)<<dis.x<<setw(20)<<dis.y<<setw(20)<<dis.z<<endl;
		}
		else if (type==1)
		{
		}
	}
	return;
}
