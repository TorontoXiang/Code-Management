#include "Tbody_ALEMPM.h"
#include "Tcell_mixed_fluid.h"
Tbody_ALEMPM::Tbody_ALEMPM(int i,int num_material,string body_name):Tbody_ALE(i,num_material,body_name)
{
	_num_particle_part=0;
	output_MPM.open("MPM_particle.dat");
	output_temp.open("temp.txt");
}
void Tbody_ALEMPM::input_body(ifstream& input)
{
	Skeyword keyword;
	read_in_keyword_file(input,keyword);
	//Read the body global information from keyword
	_endtime=keyword.time_control.endtime;
	_CFL=keyword.time_control.CFL;
	//----------------------------------------------------
	//Read the MPM particle list
	//----------------------------------------------------
	_particle_group_list=new vector<TMPM_particle>[_num_particle_part];
	int num_particle=keyword.particle_list.size();        //The total number of particles
	vec3D x_min_particle(1e10,1e10,1e10),x_max_particle(-1e10,-1e10,-1e10);
	for (int i = 0; i < num_particle; i++)
	{
		Sparticle_temp this_particle=keyword.particle_list[i];
		int group_id=this_particle.group_id;
		int part_id=this_particle.part_id;
		if (group_id>_num_particle_part-1)
		{
			cout<<"Error: The particle group id excess in ALEMPM solver!"<<endl;
			system("Pause");
			exit(0);
		}
		int id_in_part=_particle_group_list[group_id].size();
		TMPM_particle new_particle(this_particle.position,this_particle.volume,group_id,part_id,id_in_part);
		x_min_particle.x=minval(x_min_particle.x,this_particle.position.x);x_min_particle.y=minval(x_min_particle.y,this_particle.position.y);
		x_min_particle.z=minval(x_min_particle.z,this_particle.position.z);x_max_particle.x=maxval(x_max_particle.x,this_particle.position.x);
		x_max_particle.y=maxval(x_max_particle.y,this_particle.position.y);x_max_particle.z=maxval(x_max_particle.z,this_particle.position.z);
		_particle_group_list[group_id].push_back(new_particle);
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
	generate_back_ground(x_min,x_max,nx,ny,nz);
	//-------------------------------------------------------
	//Read the material properties for particles
	//-------------------------------------------------------
	int part_id,material_id,EOS_id;
	TMPM_particle* particle_ptr;
	for (int i = 0; i < _num_particle_part; i++)
	{
		part_id=_particle_group_list[i][0]._part_id-1;
		material_id=keyword.part_list[part_id].material_id;
		EOS_id=keyword.part_list[part_id].EOS_id;
		Smaterial this_mat=keyword.material_list[material_id-1];
		SEOS this_EOS;
		if (EOS_id>0)
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
	//Read the material properties for the background fluid
	//-------------------------------------------------------
	int num_fluid_property=keyword.ALEMPM_background_property.size();
	if (num_fluid_property==0)
	{
		cout<<"The fluid property must be defined in ALEMPM solver"<<endl;
		system("Pause");
		exit(0);
	}
	for (int i = 0; i < num_fluid_property; i++)
	{
		int index;
		part_id=keyword.ALEMPM_background_property[i].part_id-1;
		material_id=keyword.part_list[part_id].material_id;
		EOS_id=keyword.part_list[part_id].EOS_id;
		Smaterial this_mat=keyword.material_list[material_id-1];
		SEOS this_EOS=keyword.EOS_list[EOS_id-1];
		for (int i1 = keyword.ALEMPM_background_property[i].nx_begin; i1 < keyword.ALEMPM_background_property[i].nx_end+1; i1++)
		{
			for (int i2 = keyword.ALEMPM_background_property[i].ny_begin; i2 < keyword.ALEMPM_background_property[i].ny_end+1; i2++)
			{
				for (int i3 = keyword.ALEMPM_background_property[i].nz_begin; i3 < keyword.ALEMPM_background_property[i].nz_end+1; i3++)
				{
					index=_ny*_nz*i1+_nz*i2+i3;
					TMAT_base *mat_ptr1,*mat_ptr2;
					mat_ptr1=generate_material(this_mat,this_EOS,0);
					mat_ptr2=generate_material(this_mat,this_EOS,0);
					Tcell_fluid_base *cell1,*cell2;
					cell1=_grid1._cell_list[index];
					cell2=_grid2._cell_list[index];
					cell1->S_material_base(mat_ptr1);
					cell2->S_material_base(mat_ptr2);
				}
			}
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
						_grid1._node_list[index]->_bc_type_position[s]=maxval(keyword.MPM_boundary_condition_list[n].constraint[s],_grid1._node_list[index]->_bc_type_position[s]);
						_grid2._node_list[index]->_bc_type_position[s]=maxval(keyword.MPM_boundary_condition_list[n].constraint[s],_grid2._node_list[index]->_bc_type_position[s]);
					}
				}
			}
		}
	}
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


	//_grid_old point to grid1 and _grid_new point to grid2 initially
	_grid_old=&_grid1;_grid_new=&_grid2;

	//Get the cell connection for grid2
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=_grid2._cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=_grid2._node_list[i];
	}
	set_grid_topology();

	//Get the cell connection for grid1
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=_grid1._cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=_grid1._node_list[i];
	}
	set_grid_topology();

	//Input the multi-material if it exists
	if (_num_material>1)
	{
		input_immersed_material();
	}
	else
	{
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->reset_conservative_variable_initial_condition();
		}
	}
	//Initialize the nodal variables
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->calculate_nodal_variables();
	}
	//Clear the variables on grid2
	for (int i = 0; i < _nume; i++)
	{
		_grid2._cell_list[i]->clear_cell_variable();
	}
	for (int i = 0; i < _nump; i++)
	{
		_grid2._node_list[i]->clear_nodal_variable();
	}

	//_grid_old point to grid1 and _grid_new point to grid2 initially
	_grid_old=&_grid1;_grid_new=&_grid2;

	//Open the file for remapping information
	_iremapping.open("Remapping_info.txt");
	return;
}
void Tbody_ALEMPM::generate_back_ground(vec3D x_min,vec3D x_max,int nx,int ny, int nz)
{
	_x_min=x_min;_x_max=x_max;
	_nx=nx;_ny=ny;_nz=nz;
	//Allocate the node and cell information in sutrcuture mesh
	vec3D*** coordinate_list;
	coordinate_list=new vec3D**[_nx+1];
	for (int i = 0; i < _nx+1; i++)
	{
		coordinate_list[i]=new vec3D*[_ny+1];
	}
	for (int i = 0; i < _nx+1; i++)
	{
		for (int j = 0; j < _ny+1; j++)
		{
			coordinate_list[i][j]=new vec3D[_nz+1];
		}
	}
	cell_topology8H*** cell_list;
	cell_list=new cell_topology8H**[_nx];
	for (int i = 0; i < _nx; i++)
	{
		cell_list[i]=new cell_topology8H*[_ny];
	}
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			cell_list[i][j]=new cell_topology8H[_nz];
		}
	}
	_interval.x=(_x_max-_x_min).x/_nx;
	_interval.y=(_x_max-_x_min).y/_ny;
	_interval.z=(_x_max-_x_min).z/_nz;
	vec3D dx;
	//Generate the coordinate list
	for (int i = 0; i < _nx+1; i++)
	{
		for (int j = 0; j < _ny+1; j++)
		{
			for (int k = 0; k < _nz+1; k++)
			{
				dx.value(i*_interval.x,j*_interval.y,k*_interval.z);
				coordinate_list[i][j][k]=_x_min+dx;
			}
		}
	}
	_nump=(_nx+1)*(_ny+1)*(_nz+1);
	_nume=_nx*_ny*_nz;
	_grid1._node_list.resize(_nump);
	_grid2._node_list.resize(_nump);
	_nodeptr_list.resize(_nump);
	int index;
	for (int i = 0; i < _nx+1; i++)
	{
		for (int j = 0; j < _ny+1; j++)
		{
			for (int k = 0; k < _nz+1; k++)
			{
				index=(_ny+1)*(_nz+1)*i+(_nz+1)*j+k;
				_grid1._node_list[index]=new Tnode_fluid(index+1,coordinate_list[i][j][k]);
				_grid2._node_list[index]=new Tnode_fluid(index+1,coordinate_list[i][j][k]);
			}
		}
	}
	_grid1._cell_list.resize(_nume);
	_grid2._cell_list.resize(_nume);	
	_cellptr_list.resize(_nume);
	//Generate the cell list
	int index_node[8];
	Tnode_fluid *node_ptr1[8],*node_ptr2[8];
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			for (int k = 0; k < _nz; k++)
			{
				index=_ny*_nz*i+_nz*j+k; //The index convertion for cell
				//_cell_list[index]=cell_list[i][j][k];
				index_node[0]=(_ny+1)*(_nz+1)*i+(_nz+1)*j+k;index_node[1]=(_ny+1)*(_nz+1)*(i+1)+(_nz+1)*j+k;
				index_node[2]=(_ny+1)*(_nz+1)*(i+1)+(_nz+1)*(j+1)+k;index_node[3]=(_ny+1)*(_nz+1)*i+(_nz+1)*(j+1)+k;
				index_node[4]=(_ny+1)*(_nz+1)*i+(_nz+1)*j+k+1;index_node[5]=(_ny+1)*(_nz+1)*(i+1)+(_nz+1)*j+k+1;
				index_node[6]=(_ny+1)*(_nz+1)*(i+1)+(_nz+1)*(j+1)+k+1;index_node[7]=(_ny+1)*(_nz+1)*i+(_nz+1)*(j+1)+k+1;
				for (int n = 0; n < 8; n++)
				{
					node_ptr1[n]=_grid1._node_list[index_node[n]];
					node_ptr2[n]=_grid2._node_list[index_node[n]];
				}
				if (_num_material==1)
				{
					_grid1._cell_list[index]=new Tcell_pure_fluid_MPM(index+1,node_ptr1,NULL);
					_grid2._cell_list[index]=new Tcell_pure_fluid_MPM(index+1,node_ptr2,NULL);
				}
				else
				{
					_grid1._cell_list[index]=new Tcell_mixed_fluid_MPM(index,_num_material,node_ptr1,NULL);
					_grid2._cell_list[index]=new Tcell_mixed_fluid_MPM(index,_num_material,node_ptr2,NULL);
				}
			}
		}
	}
	//Delete the temporary mesh
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			delete cell_list[i][j];
		}
	}
	for (int i = 0; i < _nx; i++)
	{
		delete cell_list[i];
	}
	delete cell_list;
	for (int i = 0; i < _nx+1; i++)
	{
		for (int j = 0; j < _ny+1; j++)
		{
			delete coordinate_list[i][j];
		}
	}
	for (int i = 0; i < _nx+1; i++)
	{
		delete coordinate_list[i];
	}
	delete coordinate_list;
	return;
}
double Tbody_ALEMPM::calculate_time_step()
{
	double step_fluid=Tbody_ALE::calculate_time_step();
	double step_particle=calculate_time_step_particle();
	return minval(step_fluid,step_particle);
}
double Tbody_ALEMPM::calculate_time_step_particle()
{
	TMPM_particle* particle_ptr;
	int particles_in_part,cell_id;
	double max_spd,this_spd;
	for (int i = 0; i < _num_particle_part; i++)
	{
		particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			cell_id=particle_ptr->_cell_index;
			max_spd=_cellptr_list[cell_id]->G_max_spd_particle();
			this_spd=particle_ptr->_gausspoint->G_soundspeed();
			_cellptr_list[cell_id]->S_max_spd_particle(maxval(max_spd,this_spd));
		}
	}
	double time_step=100000;
	for (int i = 0; i < _nume; i++)
	{
		double cell_timestep=_cellptr_list[i]->calculate_time_step_particle();
		time_step=minval(time_step,cell_timestep);
	}
	//cout<<"clength: "<<clength<<" cps: "<<soundspeed<<" nodespd: "<<nodespeed<<" cell_id: "<<id<<endl;
	return _CFL*time_step;
	if (_current_step==0)
	{
		return _CFL*time_step;
	}
	else
	{
		return minval(_CFL*time_step,1.05*_dt);
	}
}
bool Tbody_ALEMPM::predictor_step()
{
	bool negative;
	negative=volume_check("Predict");
	if (negative)
	{
		return negative;
	}
	//Calculate the corner from fluid
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->predict_position(_dt);
	}
	calculate_corner_force(_dt);
	//Calculate the corner force from particles
	TMPM_particle* particle_ptr;
	Tcell_fluid_base* cell_ptr;
	for (int i = 0; i < _num_particle_part; i++)
	{
		int particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			cell_ptr=_cellptr_list[particle_ptr->_cell_index];
			cell_ptr->distribute_particle_force(particle_ptr);
		}
	}
	//Calculate the final nodal force
	assemble_nodal_force();
	return false;
}
bool Tbody_ALEMPM::corrector_step()
{
	bool negative=false,remapping=false;
	negative=volume_check("Correct");
	if (!negative)
	{
		#pragma omp parallel for
		for (int i = 0; i < _nump; i++)
		{
			//Calculate the velocity at the next time step
			_nodeptr_list[i]->update_vel(_dt);
			//Calculate the position at the next time point
			_nodeptr_list[i]->update_pos(_dt);
		}
		#pragma omp parallel for
		//Update cell geometry and material state
		for (int i = 0; i < _nume; i++)
		{
			if (_cellptr_list[i]->update_state(_dt))
			{
				cout<<"Error:Negative volume in the corrector step"<<endl;
				cout<<"body id= "<<_id<<" "<<"cell id= "<<_cellptr_list[i]->_id<<endl;
				system("Pause");
				exit(0);
			}
		}
		#pragma omp parallel for
		//Update the nodal variables
		for (int i = 0; i < _nump; i++)
		{
			_nodeptr_list[i]->calculate_nodal_variables();
		}
		//Update the DOF and the stress of particles
		TMPM_particle* particle_ptr;
		Tcell_fluid_base* cell_ptr;
		for (int i = 0; i < _num_particle_part; i++)
		{
			int particles_in_part=_particle_group_list[i].size();
			for (int j = 0; j < particles_in_part; j++)
			{
				particle_ptr=&_particle_group_list[i][j];
				cell_ptr=_cellptr_list[particle_ptr->_cell_index];
				cell_ptr->update_particle_DOF(particle_ptr,_dt);
				cell_ptr->update_particle_stress(particle_ptr,_dt);
			}
		}
		remapping=is_remapping();
	}
	else
	{
		cout<<"Remapping for the negative update"<<endl;
	}
	//Remap the variable if the grid is distorted
	if (remapping || negative)
	{
		output_MPM<<"remap"<<endl;
		cout<<"Remapping at step "<<_current_step<<",current time= "<<_current_time<<endl;
		_iremapping<<"Remapping at step "<<_current_step<<",current time= "<<_current_time<<endl;
		double t_begin=GetTickCount();
		remapping_variables();
		double t_end=GetTickCount();
		cout<<"This remapping phase cost: "<<(t_end-t_begin)/1000.0<<" seconds"<<endl;
		_iremapping<<"This remapping phase cost: "<<(t_end-t_begin)/1000.0<<" seconds"<<endl;
		remapping=true;
	}
	return remapping;
}
void Tbody_ALEMPM::calculate_particle_force()
{
	TMPM_particle* particle_ptr;
	Tcell_fluid_base* cell_ptr;
	for (int i = 0; i < _num_particle_part; i++)
	{
		int particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			cell_ptr=_cellptr_list[particle_ptr->_cell_index];
			cell_ptr->distribute_particle_force(particle_ptr);
		}
	}
	return;
}
void Tbody_ALEMPM::calculate_nodal_inertance()
{
	//Calculate the inertance and moment from the particles
	TMPM_particle* particle_ptr;
	Tcell_fluid_base* cell_ptr;
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->_flag=0;
	}
	for (int i = 0; i < _num_particle_part; i++)
	{
		int particle_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particle_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			find_particle_location(particle_ptr);
			cell_ptr=_cellptr_list[particle_ptr->_cell_index];
			cell_ptr->distribute_particle_variables(particle_ptr,_interval);
		}
	}
	//Calculate the initial velocity
	for (int i = 0; i < _nump; i++)
	{
		if (_nodeptr_list[i]->_mass>0)
		{
			_nodeptr_list[i]->_velocity=_nodeptr_list[i]->_moment/_nodeptr_list[i]->_mass;
			for (int j = 0; j < 3; j++)
			{
				if (_nodeptr_list[i]->_bc_type_position[j]==1)
				{
					_nodeptr_list[i]->_velocity.value(j+1,0);
				}
			}
			
		}
	}
	//Calculate the inertance from the fluid
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->assemble_corner_inertance();
	}
}
int Tbody_ALEMPM::find_particle_location(TMPM_particle* particle_ptr)
{
	vec3D particle_coordinate=particle_ptr->_position;
	for (int i = 1; i <= 3; i++)
	{
		if (particle_coordinate.access(i)<_x_min.access(i) || particle_coordinate.access(i)>_x_max.access(i))
		{
			cout<<"The particle is out of the background ground"<<endl;
			cout<<"Particleis located in "<<particle_ptr->_part_id<<"th part "<<particle_ptr->_id_in_part<<"th particle"<<endl;
			cout<<"The coordinate is: "<<particle_coordinate.x<<" "<<particle_coordinate.y<<" "<<particle_coordinate.z<<endl;
			system("Pause");
			exit(0);
		}
	}
	int ix,iy,iz;
	ix=int((particle_coordinate-_x_min).x/_interval.x);
	iy=int((particle_coordinate-_x_min).y/_interval.y);
	iz=int((particle_coordinate-_x_min).z/_interval.z);
	particle_ptr->_cell_index=(_ny*_nz)*ix+_nz*iy+iz;
	return (_ny*_nz)*ix+_nz*iy+iz;
}
void Tbody_ALEMPM::output_tecplot(ofstream& output,double ratio)
{
	Tbody_ALE::output_tecplot(output,ratio);
	MPM_tec_output();
	//TMPM_particle* particle_ptr;
	//double mass=0;
	//vec3D pos,vel;
	//for (int i = 0; i < _num_particle_part; i++)
	//{
	//	int particle_in_part=_particle_group_list[i].size();
	//	for (int j = 0; j < particle_in_part; j++)
	//	{
	//		particle_ptr=&_particle_group_list[i][j];
	//		mass=mass+particle_ptr->_mass;
	//		pos=pos+particle_ptr->_position*particle_ptr->_mass;
	//		vel=vel+particle_ptr->_velocity*particle_ptr->_mass;
	//	}
	//}
	//pos=pos/mass;vel=vel/mass;
	//output_temp<<_current_time<<" "<<pos.z<<" "<<vel.z<<endl;
	return;
}
void Tbody_ALEMPM::MPM_tec_output()
{
	TMPM_particle* particle_ptr;
	output_MPM<<"TITLE=\"Particle coordinate\""<<endl;
	output_MPM<<"Variables=\"X\",\"Y\",\"Z\",\"pressure\""<<endl;
	output_MPM<<"Zone T=\""<<_current_time<<"\""<<endl;
	for (int i = 0; i < _num_particle_part; i++)
	{
		int particles_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particles_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			output_MPM<<particle_ptr->_position.x<<" "<<particle_ptr->_position.y<<" "<<particle_ptr->_position.z;
			output_MPM<<" "<<particle_ptr->_gausspoint->G_seqv()<<endl;
		}
	}
	//double de[6];
	//for (int i = 0; i < 6; i++)
	//{
	//	de[i]=_particle_group_list[0][0]._de[i];
	//}
	//output_MPM<<_current_time<<" "<<de[0]<<" "<<de[1]<<" "<<de[2]<<" "<<de[3]<<" "<<de[4]<<" "<<de[5]<<endl;
	//output_MPM<<_current_time<<" "<<_particle_group_list[0][0]._gausspoint->G_seqv()<<endl;;
	//output_MPM<<_current_time<<" "<<_particle_group_list[0][0]._velocity.x<<_particle_group_list[0][0]._velocity.y<<" "<<_particle_group_list[0][0]._velocity.z<<endl;
	//output_MPM<<_current_time<<" "<<_cellptr_list[321]->G_gausspoint(0)->G_pressure()<<endl;
	//output_MPM<<_current_time<<" ";
	//for (int i = 0; i < 8; i++)
	//{
	//	output_MPM<<_cellptr_list[384]->_node_ptr[i]->_mass<<" ";
	//}
	//output_MPM<<endl;
	//output_MPM<<_current_time<<" "<<_particle_group_list[0][0]._gausspoint->G_seqv()<<endl;
	return;
}
void Tbody_ALEMPM::remapping_variables()
{
	//----------------------------------------prepare for remapping---------------------------------------------
	vec3D coor_min,coor_max;
	double cell_edge_max;
	//output_mesh("Final");

	_grid_smooth->generate_new_grid();
	//Smooth the old grid to generate the new grid
	reset_new_grid_geometry();
	//Reset the geometry information of the new grid
	update_old_grid_geometry();
	//Update the geometry information of the old grid
	old_grid_reconstruction();
	//Reconstruct the variable gradient and material surface on old grid
	calculate_old_grid_size(coor_min,coor_max,cell_edge_max);
	//Prepare for the bucket generation
	Tbucket_searching buckets(coor_min,coor_max,cell_edge_max);
	//Put the cells in old grid into the bucket
	bool access;
	for (int i=0;i<_nume;i++)
	{
		int ix,iy,iz;
		vec3D centroid=_grid_old->_cell_list[i]->_cell_polyhedron.G_centroid();
		buckets.search_bucket_ijk(centroid,ix,iy,iz,access);
		if (access)
		{
			cout<<"Error: The cell centroid is out of body range"<<endl;
			system("Pause");
			exit(0);
		}
		buckets.add_element(i,ix,iy,iz);
	}	
	//------------------------------remapping begin-------------------------------------------
	//----------------------------------------------------------------------------------------
	//Remap the fluid variables: calculate corner mass,corner moment,energy
	//----------------------------------------------------------------------------------------
	Tcell_fluid_base* new_cell;Tcell_fluid_base* old_cell;
	int nx_min,ny_min,nz_min,nx_max,ny_max,nz_max;
	double error=0,error_max=0;int n_error=0;
	vec3D vertex_list[8];
	vector<int> error_id;
	Tcell_intersect intersection;
	cout<<"remapping beginning!"<<endl;
	cout<<"Completing ratio: ";
	int cell_interval=_nume/minval(1000,_nume);
	bool first=true;
	double num_intersection=0;
	#pragma omp parallel for schedule(dynamic) private(intersection,nx_min,ny_min,nz_min,nx_max,ny_max,nz_max,old_cell,new_cell,vertex_list)
	for (int i=0;i<_nume;i++)
	{
		vec3D coor_min_new,coor_max_new;
		new_cell=_grid_new->_cell_list[i];
		//new_cell->calculate_coordinate_range(coor_min_new,coor_max_new);
		if (i%cell_interval==0)
		{
			if (first)
			{
				cout<<setw(4)<<100*i/_nume<<"%";
				first=false;
			}
			else
			{
				cout<<'\b'<<'\b'<<'\b'<<'\b'<<'\b';
				cout<<setw(4)<<100*i/_nume<<"%";
			}
		}
		for (int j = 0; j < 8; j++)
		{
			vertex_list[j]=new_cell->_node_ptr[j]->_position;
		}
		buckets.search_bucket_range(vertex_list,8,nx_min,ny_min,nz_min,nx_max,ny_max,nz_max,access);
		if (access)
		{
			cout<<"Error: The new grid is out of the range of old grid"<<endl;
			system("Pause");
			exit(0);
		}
		intersection.S_new_cell(new_cell);
		intersection.create_tet_list();
		for (int i1=nx_min-1;i1<=nx_max+1;i1++)
		{
			for (int i2=ny_min-1;i2<=ny_max+1;i2++)
			{
				for (int i3=nz_min-1;i3<=nz_max+1;i3++)
				{
					Sbucket* the_bucket=&buckets._bucket_list[i1][i2][i3];	
					for (int j=0;j<the_bucket->_nume;j++)
					{
						old_cell=_grid_old->_cell_list[the_bucket->eleid[j]];
						vec3D coor_min_old,coor_max_old;
						//old_cell->calculate_coordinate_range(coor_min_old,coor_max_old);
						if (!comparison_cell_coordinate(old_cell->_coor_min,old_cell->_coor_max,new_cell->_coor_min,new_cell->_coor_max))
						{
							intersection.S_old_cell(old_cell);
							intersection.intersect();
							num_intersection=num_intersection+1;
						}
					}
				}
			}
		}
	}
	cout<<endl;
	//---------------------------Calculate the fluid variables on the new grid----------------------
	double volume_error,av_volume_error=0,max_volume_error=0;
	Tcell_fluid_base* newcell_ptr=NULL;
	bool remapping_failed=false;
	for (int i=0;i<_nume;i++)
	{
		newcell_ptr=_grid_new->_cell_list[i];
		volume_error=abs(newcell_ptr->_volume_intersection-newcell_ptr->_cell_volume);
		if (volume_error>1e-10)
		{
			cout<<"Error:Remapping fails in cell  "<<i<<" with volume error= "<<volume_error<<endl;
			_iremapping<<"Error:Remapping fails in cell  "<<i<<" with volume error= "<<volume_error<<endl;
			remapping_failed=true;
		}	
		if (volume_error>max_volume_error)
		{
			max_volume_error=volume_error;
		}
		av_volume_error=av_volume_error+volume_error;
		//Update the variables in new grid
		newcell_ptr->reset_state_after_remapping();    //Nodal mass and moment from fluid are assembled here
	}
	cout<<"------Remapping result---------"<<endl;
	cout<<"The average volume error in this remapping phase is: "<<av_volume_error/_nume<<endl;
	_iremapping<<"------Remapping result---------"<<endl;
	_iremapping<<"The average volume error in this remapping phase is: "<<av_volume_error/_nume<<endl;
	_iremapping<<"The maximal volume error in this remapping phase is: "<<max_volume_error<<endl;
	_iremapping<<"The number of intersection in this remapping phase is: "<<num_intersection<<endl;
	_iremapping<<endl;
	if (remapping_failed)
	{
		system("Pause");
	}
	//----------------------------------------------------------------------------------
	//Calculate the new velocity from fluid
	//----------------------------------------------------------------------------------
	for (int i=0;i<_nump;i++)
	{
		vec3D velocity_origin=_grid_old->_node_list[i]->_velocity;
		_grid_new->_node_list[i]->reset_velocity_after_remapping(velocity_origin);
		_grid_new->_node_list[i]->calculate_nodal_variables();
		_grid_new->_node_list[i]->_flag=0;
	//	_grid_new->_node_list[i]->_moment.value(0,0,0);
	}
	//------------------------------------------------------------------------------------------
	//Accumulate the nodal mass and moment from particles
	//------------------------------------------------------------------------------------------
	TMPM_particle* particle_ptr;
	Tcell_fluid_base* cell_ptr;
	for (int i = 0; i < _num_particle_part; i++)
	{
		int particle_in_part=_particle_group_list[i].size();
		for (int j = 0; j < particle_in_part; j++)
		{
			particle_ptr=&_particle_group_list[i][j];
			find_particle_location(particle_ptr);
			cell_ptr=_grid_new->_cell_list[particle_ptr->_cell_index];
			cell_ptr->distribute_particle_variables(particle_ptr,_interval);
		}
	}
	//------------------------------------------------------------------------------------------
	//Reset the nodal velocity and apply boundary conditions form particles
	//------------------------------------------------------------------------------------------
	for (int i=0;i<_nump;i++)
	{
		_grid_new->_node_list[i]->set_velocity_from_particle_after_remapping();
	}
	//Clear the variables on the old grid
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i]->clear_cell_variable();
	}
	for (int i = 0; i < _nump; i++)
	{
		_grid_old->_node_list[i]->clear_nodal_variable();
	}

	//Change the grid pointer
	Sfluid_grid* ptr_temp=NULL;
	ptr_temp=_grid_old;_grid_old=_grid_new;_grid_new=ptr_temp;

	//Change the _cellptr_list and _nodeptr_list
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=_grid_old->_cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=_grid_old->_node_list[i];
	}
	//output_mesh("Final");
	//system("Pause");
	return;
}