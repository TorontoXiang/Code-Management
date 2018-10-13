#include "Tbody_explosion.h"
#include "Tcell_pure_fluid.h"
#include "Tcell_mixed_fluid.h"
#include "Tgrid_smooth.h"
#include "Tcell_intersect.h"
#include <omp.h>
#include <Windows.h>
#include <iomanip>
#include <cmath>
#include <iostream>
using namespace std;
Tbody_explosion::Tbody_explosion(int i,int num_material,string body_name):Tbody_ALEMPM(i,num_material,body_name)
{
	_num_remapping=0;
	_num_deactive_particle_part=0;
}

void Tbody_explosion::input_body(ifstream& input)
{
	//---------------------------------------------------------
	//Read the basic information
	//---------------------------------------------------------
	Skeyword keyword;
	read_in_keyword_file(input,keyword);
	input.close();
	//input.open("xxsample_particle1.k");
	//read_in_keyword_file(input,keyword);
	//input.close();
	_endtime=keyword.time_control.endtime;
	_CFL=keyword.time_control.CFL;
	//---------------------------------------------------------
	//Input the initial grid
	//---------------------------------------------------------
	_nume_initial=keyword.cell_8_list.size();
	_nump_initial=keyword.node_list.size();
	//Read node list
	int id;
	double x,y,z;
	_grid_initial._node_list.resize(_nump_initial);
	//_grid1._node_list.resize(_nump);
	//_grid2._node_list.resize(_nump);
	//_nodeptr_list.resize(_nump);
	cout<<"generate initial node"<<endl;
	for (int i = 0; i < _nump_initial; i++)
	{
		id=keyword.node_list[i].id;
		x=keyword.node_list[i].x;
		y=keyword.node_list[i].y;
		z=keyword.node_list[i].z;
		Tnode_fluid new_node(id,x,y,z);
		_grid_initial._node_list[i]=new Tnode_fluid(id,x,y,z);
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	cout<<"generate initial cell"<<endl;
	_grid_initial._cell_list.resize(_nume_initial);
	int part_id,material_id,EOS_id;
	Tnode_fluid* node_ptr_initial[8];
	int num_part=keyword.part_list.size()-keyword.num_particle_part;
	if (num_part!=_num_material)
	{
		cout<<"Error: The number of part is not equal to the number of material in Tbody_explosion"<<endl;
		system("Pause");
		exit(0);
	}
	for (int i = 0; i < _nume_initial; i++)
	{
		//Read cell id and IEN
		id=keyword.cell_8_list[i].cell_id;
		for (int j = 0; j < 8; j++)
		{
			int IENj=keyword.cell_8_list[i].IEN[j];
			node_ptr_initial[j]=_grid_initial._node_list[IENj-1];
		}
		//Allocate the materials in this problem
		TMAT_base** mat;
		mat=new TMAT_base*[_num_material];
		for (int j = 0; j < _num_material; j++)
		{
			EOS_id=keyword.part_list[j].EOS_id;
			if (EOS_id==0 || keyword.EOS_list.size()==0)
			{
				cout<<"Error:There must be an EOS in body_explosion"<<endl;
				system("Pause");
				exit(0);
			}
			material_id=keyword.part_list[j].material_id;
			Smaterial this_mat=keyword.material_list[material_id-1];
			SEOS this_EOS=keyword.EOS_list[EOS_id-1];
			mat[j]=generate_material(this_mat,this_EOS,1);
		}
		//Create a new cell
		_grid_initial._cell_list[i]=new Tcell_mixed_fluid(id,_num_material,node_ptr_initial,mat);


		//Initialize the volume fraction and material centroid
		part_id=keyword.cell_8_list[i].part_id;
		material_id=keyword.part_list[part_id-1].material_id;
		_grid_initial._cell_list[i]->S_material_fraction(1,material_id-1);
		_grid_initial._cell_list[i]->S_material_centroid(_grid_initial._cell_list[i]->G_cell_polyhedron()->G_centroid(),material_id-1);
	}
	//-------------------------------------------------------
	//Read boundary condition
	//-------------------------------------------------------
	cout<<"generate initial bd"<<endl;
	int node_group_id;   //Boundary condition (load) is applied on node_group_id

	int num_bd_group;    //Number of boundary condition group
	int num_node;        //Number of nodes in the node group
	num_bd_group=keyword.boundary_list.size();
	for (int i = 0; i < num_bd_group; i++)
	{
		node_group_id=keyword.boundary_list[i].id;
		num_node=keyword.node_group_list[node_group_id-1].node_id.size();
		for (int j = 0; j < 3; j++)
		{
			if (keyword.boundary_list[i].pos[j]==1)
			{
				for (int k = 0; k < num_node; k++)
				{
					_grid_initial._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[j]=1;
				}
			}
		}
	}
	//keyword.clear_keyword();
	//input.close();    //End the input of the initial grid
	//------------------------------------------------------------------
	//Generate the background grid
	//------------------------------------------------------------------
	cout<<"generate background grid"<<endl;
	vec3D x_min=keyword.background.x_min;
	vec3D x_max=keyword.background.x_max;
	int nx=keyword.background.nx;
	int ny=keyword.background.ny;
	int nz=keyword.background.nz;
	generate_back_ground(x_min,x_max,nx,ny,nz);
	//------------------------------------------------------------------
	//Set the material for the background grid
	//------------------------------------------------------------------
	for (int i = 0; i < _nume_com; i++)
	{
		TMAT_base** mat1;
		TMAT_base** mat2;
		mat1=new TMAT_base*[_num_material];
		mat2=new TMAT_base*[_num_material];
		for (int j = 0; j < _num_material; j++)
		{
			EOS_id=keyword.part_list[j].EOS_id;
			if (EOS_id==0 || keyword.EOS_list.size()==0)
			{
				cout<<"Error:There must be an EOS in body_explosion"<<endl;
				system("Pause");
				exit(0);
			}
			material_id=keyword.part_list[j].material_id;
			Smaterial this_mat=keyword.material_list[material_id-1];
			SEOS this_EOS=keyword.EOS_list[EOS_id-1];
			mat1[j]=generate_material(this_mat,this_EOS,1);
			mat2[j]=generate_material(this_mat,this_EOS,1);
			_grid1._cell_list[i]->S_material_base(mat1);
			_grid2._cell_list[i]->S_material_base(mat2);
		}
	}
	//----------------------------------------------------------------------
	//Read the particle information
	//----------------------------------------------------------------------
	cout<<"generate particles"<<endl;
	_num_particle_part=keyword.num_particle_part;
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
	//-------------------------------------------------------------------------
	//Read the particle material
	//-------------------------------------------------------------------------
//	int part_id,material_id,EOS_id;
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
	_num_deactive_particle_part=_num_particle_part;
	_num_particle_part=0;
	//-----------------------------------------------------------------------
	//Read the boundary condition of the background grid
	//-----------------------------------------------------------------------
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

	//Allocate the _cellptr_list and _nodeptr_list
	_cellptr_list.resize(maxval(_nume_com,_nume_initial));
	_nodeptr_list.resize(maxval(_nump_com,_nump_initial));

	//Calculate the topology of the computational grid
	_nume=_nume_com;_nump=_nump_com;
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=_grid1._cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=_grid1._node_list[i];
	}
	set_grid_topology();
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=_grid2._cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=_grid2._node_list[i];
	}
	set_grid_topology();

	vec2D center;
	center.value(0,0);
	generate_cylinder_particles(0.25,0.1,0.11,center);
	//Clear the variables on the computational grid
	for (int i = 0; i < _nume; i++)
	{
		_grid1._cell_list[i]->clear_cell_variable();
		_grid2._cell_list[i]->clear_cell_variable();
	}
	for (int i = 0; i < _nump; i++)
	{
		_grid1._node_list[i]->clear_nodal_variable();
		_grid2._node_list[i]->clear_nodal_variable();
	}

	//Calculate the topology of the initial grid
	_nume=_nume_initial;_nump=_nump_initial;
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=_grid_initial._cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=_grid_initial._node_list[i];
	}
	set_grid_topology();	

	//Initialize the initial grid
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->reset_conservative_variable_initial_condition();
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->calculate_nodal_variables();
	}

	//Point the grid_old and grid_new
	_grid_old=&_grid_initial;_grid_new=&_grid1;

	//Open the file for remapping information
	_iremapping.open("Remapping_info.txt");
	return;
}
void Tbody_explosion::remapping_variables()
{
	//----------------------------------------prepare for remapping---------------------------------------------
	_num_remapping=_num_remapping+1;
	vec3D coor_min,coor_max;
	double cell_edge_max;
	//output_mesh("Initial");
	//system("Pause");
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
	Tcell_fluid_base* new_cell;Tcell_fluid_base* old_cell;
	int nx_min,ny_min,nz_min,nx_max,ny_max,nz_max;
	double error=0,error_max=0;int n_error=0;
	vec3D vertex_list[8];
	vector<int> error_id;
	Tcell_intersect intersection;
	cout<<"remapping beginning!"<<endl;
	cout<<"Completing ratio: ";	
	bool first=true;
	double num_intersection=0;
	_nume=_grid_new->_cell_list.size();
	_nump=_grid_new->_node_list.size();
	int cell_interval=maxval(1,_nume/1000);
	_grid_smooth->generate_new_grid();
	//Smooth the old grid to generate the new grid
	reset_new_grid_geometry();
	//Reset the geometry information of the new grid
	double time_begin=GetTickCount();
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
	double time_end=GetTickCount();
	cout<<endl;
	//------------------Update the variables on new grid and out put the remapping result------------
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
		newcell_ptr->reset_state_after_remapping();
	}
	cout<<"------Remapping result---------"<<endl;
	cout<<"The average volume error in this remapping phase is: "<<av_volume_error/_nume<<endl;
	cout<<"The time consuming is: "<<(time_end-time_begin)/1000<<"seconds"<<endl;
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
	_nume=_grid_old->_cell_list.size();
	_nump=_grid_old->_node_list.size();
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
	if (_num_remapping==1)
	{
		ptr_temp=_grid_old;_grid_old=_grid_new;_grid_new=&_grid2;
	}
	else
	{
		ptr_temp=_grid_old;_grid_old=_grid_new;_grid_new=ptr_temp;
	}
	//Change the _cellptr_list and _nodeptr_list
	_nume=_grid_old->_cell_list.size();
	_nump=_grid_old->_node_list.size();
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=_grid_old->_cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=_grid_old->_node_list[i];
	}
	//system("Pause");
	//output_mesh("Initial");
	return;	
}
bool Tbody_explosion::is_remapping()
{
	if (_current_time<_detonation_time && _current_time+_dt>=_detonation_time)
	{
		_num_particle_part=_num_deactive_particle_part;
		return true;
	}
	if (_current_time<_detonation_time)
	{
		return false;
	}
	Tbody_ALE::is_remapping();
}
void Tbody_explosion::generate_back_ground(vec3D x_min,vec3D x_max,int nx,int ny, int nz)
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
	_nump_com=(_nx+1)*(_ny+1)*(_nz+1);
	_nume_com=_nx*_ny*_nz;
	_grid1._node_list.resize(_nump_com);
	_grid2._node_list.resize(_nump_com);
	_nodeptr_list.resize(_nump_com);
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
	_grid1._cell_list.resize(_nume_com);
	_grid2._cell_list.resize(_nume_com);	
	//_cellptr_list.resize(_nume);
	//Generate the cell list
	int index_node[8];
	Tnode_fluid *node_ptr1[8],*node_ptr2[8];
	TMAT_base* matptr=NULL;
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
				_grid1._cell_list[index]=new Tcell_mixed_fluid_MPM(index,_num_material,node_ptr1,matptr);
				_grid2._cell_list[index]=new Tcell_mixed_fluid_MPM(index,_num_material,node_ptr2,matptr);
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
void Tbody_explosion::S_remapping_scheme(Sremapping_scheme &remapping_scheme)
{
	_remapping_shceme=remapping_scheme;
	if (_remapping_shceme.name=="Times")
	{
		int times=_remapping_shceme.num_remapping;
		double dt=(_endtime-_detonation_time)/(times+1);
		_remapping_shceme.remapping_time_point.resize(times+1);
		for (int i = 0; i < times; i++)
		{
			_remapping_shceme.remapping_time_point[i]=_detonation_time+(i+1)*dt;
		}
	}
	return;
}
void Tbody_explosion::generate_cylinder_particles(double r,double z_min,double z_max,vec2D center)
{
	Tcell_fluid_base* cell_ptr;
	Tcell_fluid_base* adjacent_cell_ptr;
	vec3D cell_center;
	vec2D xy;
	double z;
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->_flag=0;
	}
	for (int i = 0; i < _nume; i++)
	{
		cell_ptr=_cellptr_list[i];
		cell_ptr->_flag=0;
		cell_center=cell_ptr->_cell_polyhedron.calculate_polyhedron_centroid();
		xy.value(cell_center.x,cell_center.y);
		z=cell_center.z;
		if ((xy-center).get_length()<=r && z>=z_min && z<=z_max)
		{
			cell_ptr->_flag=1;
		}
	}
	for (int i = 0; i < _nume; i++)
	{
		cell_ptr=_cellptr_list[i];
		if (cell_ptr->_flag==1)
		{
			for (int j = 0; j < 6; j++)
			{
				adjacent_cell_ptr=cell_ptr->_adjacent[j];
				if (adjacent_cell_ptr!=NULL)
				{
					cell_center=adjacent_cell_ptr->_cell_polyhedron.G_centroid();
					xy.value(cell_center.x,cell_center.y);
					if ((xy-center).get_length()>r)
					{
						for (int k = 0; k < 4; k++)
						{
							cell_ptr->_node_ptr[Tcell_fluid_base::cell_face[j][k]]->_flag=1;
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < _nume; i++)
	{
		cell_ptr=_cellptr_list[i];
		cell_ptr->_flag=0;
		cell_center=cell_ptr->_cell_polyhedron.G_centroid();
		xy.value(cell_center.x,cell_center.y);
		z=cell_center.z;
		if ((xy-center).get_length()>r && z>=z_min && z<=z_max)
		{
			for (int j = 0; j < 8; j++)
			{
				//cell_ptr->_node_ptr[j]->_flag=1;
			}
		}
	}
	vector<int> fixed_node_list;
	for (int i = 0; i < _nump; i++)
	{
		if (_nodeptr_list[i]->_flag==1)
		{
			fixed_node_list.push_back(_nodeptr_list[i]->_id);
		}
	}
	int num_fixed_node=fixed_node_list.size();
	for (int i = 0; i < num_fixed_node; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_grid1._node_list[fixed_node_list[i]]->_bc_type_position[j]=1;
			_grid2._node_list[fixed_node_list[i]]->_bc_type_position[j]=1;
		}
	}
	return;
	ofstream output;
	output.open("fixed_node_test.dat");
	output<<"TITLE=\"fixed node coordinate\""<<endl;
	output<<"Variables=\"X\",\"Y\",\"Z\""<<endl;
	for (int i = 0; i < num_fixed_node; i++)
	{
//		if (_nodeptr_list[i]->_flag==1)
//		{
			vec3D coor=_grid1._node_list[fixed_node_list[i]]->_position;
			output<<coor.x<<" "<<coor.y<<" "<<coor.z<<endl;
//		}
	}
	system("Pause");
	exit(0);
}