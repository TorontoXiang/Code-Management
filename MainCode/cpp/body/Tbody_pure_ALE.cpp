#include "Tbody_pure_ALE.h"
#include "Tgrid_smooth.h"
#include "Tgrid_smooth_continum.h"
#include "Tcell_intersect.h"
#include <omp.h>
#include <Windows.h>
#include <iomanip>
#include <cmath>
#include <iostream>
using namespace std;
Tbody_pure_ALE::Tbody_pure_ALE(int i,string body_name):Tbody_pure_Lagrangian(i,body_name)
{
	_grid_new=_grid_old=NULL;
};
bool Tbody_pure_ALE::corrector_step()
{
	Tbody_pure_Lagrangian::corrector_step();
	bool remapping=is_remapping();
	if (remapping)
	{
		cout<<"Remapping at step "<<_current_step<<",current time= "<<_current_time<<endl;
		_iremapping<<"Remapping at step "<<_current_step<<",current time= "<<_current_time<<endl;
		double t_begin=GetTickCount();
		remapping_variables();
		double t_end=GetTickCount();
		cout<<"This remapping phase cost: "<<(t_end-t_begin)/1000.0<<" seconds"<<endl;
		_iremapping<<"This remapping phase cost: "<<(t_end-t_begin)/1000.0<<" seconds"<<endl;
	}
	return remapping;
}
void Tbody_pure_ALE::S_remapping_scheme(Sremapping_scheme &remapping_scheme)
{
	_remapping_shceme=remapping_scheme;
	if (_remapping_shceme.name=="Times")
	{
		int times=_remapping_shceme.num_remapping;
		double dt=_endtime/(times+1);
		_remapping_shceme.remapping_time_point.resize(times);
		for (int i = 0; i < times; i++)
		{
			_remapping_shceme.remapping_time_point[i]=(i+1)*dt;
		}
	}
	return;
}
void Tbody_pure_ALE::S_remesh_scheme(string name)
{
	if (name=="Origin")
	{
		_grid_smooth=new Tgrid_smooth_origin(this);
	}
	else if (name=="Continum_analogy")
	{
		_grid_smooth=new Tgrid_smooth_continum(this);
	}
	else
	{
		cout<<"Error:Invalid remesh scheme name"<<endl;
		system("Pause");
		exit(0);
	}
	return;
}
bool Tbody_pure_ALE::is_remapping()
{
	if (_remapping_shceme.name=="None")
	{
		return false;
	}
	else if (_remapping_shceme.name=="Times")
	{
		int times=_remapping_shceme.num_remapping;
		for (int i = 0; i <_remapping_shceme.num_remapping; i++)
		{
			double remapping_time=_remapping_shceme.remapping_time_point[i];
			if (_current_time<remapping_time && _current_time+_dt>remapping_time)
			{
				_remapping_shceme.remapping_time_point[i]=-1;
				return true;
			}
		}
		return false;
	}
	else
	{
		double ratio_angle,ratio_length;
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->calculate_distorted_ratio(ratio_angle,ratio_length);
			if (ratio_angle>=_remapping_shceme.tolerance_angle || ratio_length>=_remapping_shceme.tolerance_length)
			{
				return true;
			}
		}
		return false;
	}
}
void Tbody_pure_ALE::calculate_gradient()
{
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i].calculate_gradient_in_cell();
		//_grid_old->_cell_list[i].limit_gradient_in_cell();
	}
}
void Tbody_pure_ALE::remapping_variables()
{
	//----------------------------------------prepare for remapping---------------------------------------------
	vec3D coor_min,coor_max;
	double cell_edge_max;
	output_mesh("Final");
	_grid_smooth->generate_new_grid();
	//Smooth the old grid to generate the new grid
	plot_grid_new();
	reset_new_grid_geometry();
	//Reset the geometry information of the new grid
	update_old_grid_geometry();
	//Update the geometry information of the old grid
	old_grid_reconstruction();
	//Reconstruct the variable gradient on old grid
	calculate_old_grid_size(coor_min,coor_max,cell_edge_max);
	//Prepare for the bucket generation
	Tbucket_searching buckets(coor_min,coor_max,cell_edge_max);
	//Put the cells in old grid into the bucket
	bool access;
	for (int i=0;i<_nume;i++)
	{
		int ix,iy,iz;
		vec3D centroid=_grid_old->_cell_list[i]._cell_centroid;
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
	Tcell_pure_ALE* new_cell;Tcell_pure_ALE* old_cell;
	int nx_min,ny_min,nz_min,nx_max,ny_max,nz_max;
	double error=0,error_max=0;int n_error=0;
	vec3D vertex_list[8];
	vector<int> error_id;
	Tcell_intersect intersection;
	cout<<"remapping beginning!"<<endl;
	cout<<"Completing ratio: ";
	int cell_interval=_nume/20;
	bool first=true;
	#pragma omp parallel for schedule(dynamic) private(intersection,nx_min,ny_min,nz_min,nx_max,ny_max,nz_max,old_cell,new_cell,vertex_list)
	for (int i=0;i<_nume;i++)
	{
		new_cell=&_grid_new->_cell_list[i];
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
						old_cell=&_grid_old->_cell_list[the_bucket->eleid[j]];
						if (!comparison_cell_coordinate(old_cell->_coor_min,old_cell->_coor_max,new_cell->_coor_min,new_cell->_coor_max))
						{
							intersection.S_old_cell(old_cell);
							intersection.intersect();
						}
					}
				}
			}
		}
	}
	cout<<endl;
	//------------------Update the variables on new grid and out put the remapping result------------
	double volume_error,av_volume_error=0;
	Tcell_pure_ALE* newcell_ptr=NULL;
	bool remapping_failed=false;
	for (int i=0;i<_nume;i++)
	{
		newcell_ptr=&_grid_new->_cell_list[i];
		volume_error=abs(newcell_ptr->_volume_intersection-newcell_ptr->_cell_volume);
		if (volume_error>1e-10)
		{
			cout<<"Error:Remapping fails in cell  "<<i<<" with volume error= "<<volume_error<<endl;
			_iremapping<<"Error:Remapping fails in cell  "<<i<<" with volume error= "<<volume_error<<endl;
			remapping_failed=true;
		}		
		av_volume_error=av_volume_error+volume_error;
		//Update the variables in new grid
		double density=newcell_ptr->_cell_mass/newcell_ptr->_cell_volume;
		double internal_energy=newcell_ptr->_cell_energy/(newcell_ptr->_cell_volume*density);
		newcell_ptr->_gausspoint[0]->S_density(density);
		newcell_ptr->_gausspoint[0]->S_internal_energy(internal_energy);
		double pressure=newcell_ptr->_gausspoint[0]->calculate_pressure();
		double soundspeed=newcell_ptr->_gausspoint[0]->calculate_sound_speed();
		newcell_ptr->_gausspoint[0]->S_pressure(pressure);
		newcell_ptr->_gausspoint[0]->S_soundspeed(soundspeed);
		newcell_ptr->assemble_nodal_mass_moment();
		newcell_ptr->calculate_corner_mass();
	}
	cout<<"------Remapping result---------"<<endl;
	cout<<"The average volume error in this remapping phase is: "<<av_volume_error/_nume<<endl;
	_iremapping<<"------Remapping result---------"<<endl;
	_iremapping<<"The average volume error in this remapping phase is: "<<av_volume_error/_nume<<endl;
	if (remapping_failed)
	{
		system("Pause");
	}
	//Reset the nodal velocity and apply boundary conditions on new grid
	//Calculate the nodal variables
	for (int i=0;i<_nump;i++)
	{
		vec3D velocity_origin=_grid_old->_node_list[i]._velocity;
		_grid_new->_node_list[i].reset_velocity_after_remapping(velocity_origin);
		_grid_new->_node_list[i].calculate_nodal_variables();
	}

	//Clear the variables on the old grid
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i].clear_cell_variable();
	}
	for (int i = 0; i < _nump; i++)
	{
		_grid_old->_node_list[i].clear_nodal_variable();
	}

	//Change the grid pointer
	Sgrid<Tcell_pure_ALE,Tnode_fluid>* ptr_temp=NULL;
	ptr_temp=_grid_old;_grid_old=_grid_new;_grid_new=ptr_temp;

	//Change the _cellptr_list and _nodeptr_list
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=&_grid_old->_cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=&_grid_old->_node_list[i];
	}
	return;
}
void Tbody_pure_ALE::reset_new_grid_geometry()
{
	bool negative;
	for (int i = 0; i < _nume; i++)
	{
		_grid_new->_cell_list[i].calculate_coordinate_range();
		//Reset the cell coordinate range of new grid
		_grid_new->_cell_list[i].calculate_cell_volume("Full",negative);
		if (negative)
		{
			cout<<"Error: Negative volume in calling reset_new_grid_geometry()"<<endl;
			cout<<"body id= "<<_id<<" "<<"cell id= "<<_grid_new->_cell_list[i]._id<<endl;
		}
		//Reset cell volume and corner volume of new grid
		_grid_new->_cell_list[i].update_cell_polyhedron();
		//Reset the cell polyhedron of new grid
		_grid_new->_cell_list[i]._cell_centroid=_grid_new->_cell_list[i]._cell_polyhedron.G_centroid();
		//Reset the cell centroid of new grid
	}
	return;
}
void Tbody_pure_ALE::update_old_grid_geometry()
{
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i].calculate_coordinate_range();
		//Update the cell coordinate range of old grid
		_grid_old->_cell_list[i].update_cell_polyhedron();
		//Update the cell polyhedron of old grid
		_grid_old->_cell_list[i]._cell_centroid=_grid_old->_cell_list[i]._cell_polyhedron.G_centroid();
		_grid_old->_cell_list[i]._cell_volume=_grid_old->_cell_list[i]._cell_polyhedron.G_volume();
		//Update the cell centroid of old grid
	}
	return;
}
void Tbody_pure_ALE::old_grid_reconstruction()
{
	//Reconstruct the gradient of variables
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i].calculate_gradient_in_cell();
		_grid_old->_cell_list[i].limit_gradient_in_cell();
	}
	return;
}
void Tbody_pure_ALE::calculate_old_grid_size(vec3D &coor_min,vec3D &coor_max,double &cell_edge_max)
{
	cell_edge_max=0;
	double l_ele;
	for (int i=0;i<_nume;i++)
	{
		l_ele=_grid_old->_cell_list[i].calculate_maximal_edge();
		if (l_ele>cell_edge_max)
		{
			cell_edge_max=l_ele;
		}
	}
	coor_min.value(1e10,1e10,1e10);coor_max.value(-1e10,-1e10,-1e10);
	for (int i=0;i<_nump;i++)
	{
		for (int j = 1; j < 4; j++)
		{
			if (_grid_old->_node_list[i]._position.access(j)<coor_min.access(j)) coor_min.value(j,_grid_old->_node_list[i]._position.access(j));
			if (_grid_old->_node_list[i]._position.access(j)>coor_max.access(j)) coor_max.value(j,_grid_old->_node_list[i]._position.access(j));
		}
	}
	return;
}
void Tbody_pure_ALE::input_body(ifstream &input)
{
	Skeyword keyword;
	read_in_keyword_file(input,keyword);
	//Read the body global information from keyword
	_nume=keyword.cell_8_list.size();
	_nump=keyword.node_list.size();
	_endtime=keyword.time_control.endtime;
	_CFL=keyword.time_control.CFL;
	//----------------------------------------------------
	//Read node list
	//----------------------------------------------------
	int id;
	double x,y,z;
	_grid1._node_list.resize(_nump);
	_grid2._node_list.resize(_nump);
	_nodeptr_list.resize(_nump);
	for (int i = 0; i < _nump; i++)
	{
		id=keyword.node_list[i].id;
		x=keyword.node_list[i].x;
		y=keyword.node_list[i].y;
		z=keyword.node_list[i].z;
		Tnode_fluid new_node(id,x,y,z);
		_grid1._node_list[i]=new_node;
		_grid2._node_list[i]=new_node;
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	_grid1._cell_list.resize(_nume);
	_grid2._cell_list.resize(_nume);
	_cellptr_list.resize(_nume);
	int nGauss;
	int part_id,material_id,EOS_id;
	Tnode_fluid* node_ptr1[8];
	Tnode_fluid* node_ptr2[8];
	for (int i = 0; i < _nume; i++)
	{
		//Read cell id and IEN
		id=keyword.cell_8_list[i].cell_id;
		for (int j = 0; j < 8; j++)
		{
			int IENj=keyword.cell_8_list[i].IEN[j];
			node_ptr1[j]=&_grid1._node_list[IENj-1];
			node_ptr2[j]=&_grid2._node_list[IENj-1];
		}
		//Read cell properties (nGauss and thickness)
		part_id=keyword.cell_8_list[i].part_id;
		nGauss=1;
		//Read cell material
		material_id=keyword.part_list[part_id-1].material_id;
		EOS_id=keyword.part_list[part_id-1].EOS_id;
		if (EOS_id==0 || keyword.EOS_list.size()==0)
		{
			cout<<"Error:There must be an EOS in body_pure_ALE"<<endl;
			system("Pause");
			exit(0);
		}
		TMAT_base** mat1;
		TMAT_base** mat2;
		mat1=new TMAT_base*[nGauss];mat2=new TMAT_base*[nGauss];
		for (int k = 0; k < nGauss; k++)
		{
			Smaterial this_mat=keyword.material_list[material_id-1];
			SEOS this_EOS=keyword.EOS_list[EOS_id-1];
			mat1[k]=generate_material(this_mat,this_EOS,1);
			mat2[k]=generate_material(this_mat,this_EOS,1);
		}
		//Create a new cell
		Tcell_pure_ALE new_cell1(id,nGauss,node_ptr1,mat1);
		Tcell_pure_ALE new_cell2(id,nGauss,node_ptr2,mat2);
		_grid1._cell_list[i]=new_cell1;
		_grid2._cell_list[i]=new_cell2;
	}
	//-------------------------------------------------------
	//Read boundary condition
	//-------------------------------------------------------
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
					_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]._bc_type_position[j]=1;
					_grid2._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]._bc_type_position[j]=1;
				}
			}
		}
	}
	//---------------------------------------------------------------------------------
	//Read external load
	//---------------------------------------------------------------------------------
	int num_load_group;          //Number of load group
	int curve_id;                //The curve id of the load
	int direction;               //The load direction
	double magnitude;            //The load magnitude
	Sexternal_load new_load;
	num_load_group=keyword.load_list.size();
	for (int i = 0; i < num_load_group; i++)
	{
		node_group_id=keyword.load_list[i].id;
		new_load.node_id=keyword.node_group_list[node_group_id-1].node_id;
		new_load.num_node=new_load.node_id.size();
		direction=keyword.load_list[i].dof;
		new_load.direction=direction;
		curve_id=keyword.load_list[i].curve_id;
		magnitude=keyword.curve_list[curve_id-1].v2[0];      //Only for constant load
		new_load.magnitude=magnitude;
		_external_load_list.push_back(new_load);
	}

	//Get the cell connection for grid2
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=&_grid2._cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=&_grid2._node_list[i];
	}
	set_grid_topology();

	//Get the cell connection for grid1
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=&_grid1._cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=&_grid1._node_list[i];
	}
	set_grid_topology();

	//Initialize the nodal variables
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->calculate_nodal_variables();
	}

	//Clear the variables on grid2
	for (int i = 0; i < _nume; i++)
	{
		_grid2._cell_list[i].clear_cell_variable();
	}
	for (int i = 0; i < _nump; i++)
	{
		_grid2._node_list[i].clear_nodal_variable();
	}

	//_grid_old point to grid1 and _grid_new point to grid2 initially
	_grid_old=&_grid1;_grid_new=&_grid2;

	//Open the file for remapping information
	_iremapping.open("Remapping_info.txt");
	return;
}
//-------------------------------------------------------------------------------
//Functions for test
//-------------------------------------------------------------------------------

void Tbody_pure_ALE::set_velocity(int type)
{
	if (type==1)
	{
		for (int i = 0; i < _nump; i++)
		{
			vec3D pos=_grid_old->_node_list[i]._position;
			//_grid_old->_node_list[i]._velocity.x=0.25-(pos.x-0.5)*(pos.x-0.5);
			//_grid_old->_node_list[i]._velocity.y=0.25-(pos.y-0.5)*(pos.y-0.5);
			_grid_old->_node_list[i]._velocity.z=pos.x+2*pos.y+3*pos.z;
		}
	}
	else if (type==2)
	{
		vec3D temp1,temp2;
		double l_max,l_min;
		calculate_body_size(temp1,temp2,l_max,l_min);
		for (int i = 0; i < _nump; i++)
		{
			vec3D pos=_grid_old->_node_list[i]._position;
			_grid_old->_node_list[i]._velocity.x=random(-0.5,0.5)*l_max;
			_grid_old->_node_list[i]._velocity.y=random(-0.5,0.5)*l_max;
			_grid_old->_node_list[i]._velocity.z=random(-0.5,0.5)*l_max;
			if (abs((pos.x-1)*pos.x)<1e-10)
			{
				_grid_old->_node_list[i]._velocity.x=0;
			}
			if (abs((pos.y-1)*pos.y)<1e-10)
			{
				_grid_old->_node_list[i]._velocity.y=0;
			}
			if (abs((pos.z-1)*pos.z)<1e-10)
			{
				_grid_old->_node_list[i]._velocity.z=0;
			}
		}
	}
	return;
}
void Tbody_pure_ALE::update_position(double dt)
{
	for (int i = 0; i < _nump; i++)
	{
		vec3D pos=_grid_old->_node_list[i]._position;
		vec3D vel=_grid_old->_node_list[i]._velocity;
		_grid_old->_node_list[i]._position=pos+vel*dt;
	}
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i].update_cell_polyhedron();
		_grid_old->_cell_list[i]._cell_centroid=_grid_old->_cell_list[i]._cell_polyhedron.G_centroid();
	}
	return;
}
void Tbody_pure_ALE::set_variables()
{
	for (int i = 0; i < _nume; i++)
	{
		vec3D cell_entroid=_grid_old->_cell_list[i]._cell_centroid;
		vec3D g_ro(5,6,7),g_roe(5,6,7);
		double ro=1+g_ro*cell_entroid;
		double roe=1+g_roe*cell_entroid;
		_grid_old->_cell_list[i]._gausspoint[0]->S_density(ro);
		_grid_old->_cell_list[i]._gausspoint[0]->S_internal_energy(roe/ro);
	}
	return;
}
void Tbody_pure_ALE::iteration_test()
{
	for (int i = 0; i < _nume; i++)
	{
		vec3D coordinate_std(random(-1,1),random(-1,1),random(-1,1));
		vec3D coordinate_physical=_grid_old->_cell_list[i].shape_function_interpolate(coordinate_std,"Position");
		vec3D result=_grid_old->_cell_list[i].calculate_std_coordinate(coordinate_physical);
		double error=sqrt((result-coordinate_std).self_multuply());
		if (error>1e-10)
		{
			cout<<i<<" "<<error<<endl;
		}
	}
	return;
}
void Tbody_pure_ALE::set_displacemant_at_conode(double dl)
{
	vec3D n;
	output_mesh("Final");
	for (int i = 0; i < _nump; i++)
	{
		if (_nodeptr_list[i]->_is_conode==1)
		{
			n=_nodeptr_list[i]->_position0/_nodeptr_list[i]->_position.get_length();
			_nodeptr_list[i]->_position=_nodeptr_list[i]->_position0+n*dl;
		}
	}
	output_mesh("Final");
	_grid_smooth->generate_new_grid();
	plot_grid_new();
	return;
}
void Tbody_pure_ALE::plot_grid_new()
{
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=&_grid_new->_cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=&_grid_new->_node_list[i];
	}
	output_mesh("Final");
	return;
}
void Tbody_pure_ALE::test_for_fraction_below()
{
	Tvolume_equation_solver volume_solver1;
	Tvolume_equation_solver volume_solver2;
	set_velocity(2);
	update_position(0.1);
	Tpolyhedron* polyhedron_ptr;
	double theta,phy,d,d_min,d_max,ratio;
	vec3D n;
	vec3D result1,result2;
	double fraction_error,derivative_error;
	double fraction_error_max=0,derivative_error_max=0;
	theta=phy=1;ratio=0.1;
	int n_divide=211;
	double pi=3.141592654,dtheta=pi/n_divide,dphy=2*pi/n_divide,dratio=1.0/n_divide;
	ios::scientific;
	ofstream output;
	output.open("polyhedron.dat");
	for (int i = 0; i <_nume; i++)
	{
		polyhedron_ptr=&_grid_old->_cell_list[i]._cell_polyhedron;
		polyhedron_ptr->plot_this_polyehdron(output);
		volume_solver1.set_target_polyhedron(polyhedron_ptr);
		volume_solver2.set_target_polyhedron(polyhedron_ptr);
		for (int j = 0; j < n_divide; j++)
		{
			for (int k = 0; k < n_divide; k++)
			{
				theta=j*dtheta;phy=k*dphy;
				n.value(sin(theta)*cos(phy),sin(theta)*sin(phy),cos(theta));
				polyhedron_ptr->calculate_altitude_range(n,d_min,d_max);
				for (int l = 1; l < n_divide-1; l++)
				{
					ratio=l*dratio;
					d=d_min+ratio*(d_max-d_min);
					result1=volume_solver1.calculate_fraction_below(n,d);
					result2=volume_solver2.calculate_fraction_below(n,d);
					fraction_error=abs(result1.x-result2.x);
					derivative_error=abs(result1.y-result2.y);
					if (fraction_error>1e-5)
					{
						output<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
						output<<result1.x<<" "<<result2.x<<endl;
						output<<fraction_error<<endl;
						cout<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
						cout<<result1.x<<" "<<result2.x<<endl;
						cout<<fraction_error<<endl;
						output<<endl;
					}
					fraction_error_max=maxval(fraction_error_max,fraction_error);
					derivative_error_max=maxval(derivative_error_max,derivative_error);
				}
			}
		}
	}
	cout<<"The maximal fraction error is "<<fraction_error_max<<endl;
	cout<<"The maximal derivative error is "<<derivative_error_max<<endl;
	return;
}
void Tbody_pure_ALE::test_for_volume_solver()
{
	vec3D p1(0,0,0),p2(0,1,0),p3(0,0,1),p4(1,0,0);
	Tpolyhedron poly;
	poly.add_vertex(p1);poly.add_vertex(p2);poly.add_vertex(p3);poly.add_vertex(p4);
	poly.add_piece(0,1,3);poly.add_piece(1,2,3);poly.add_piece(0,2,1);poly.add_piece(0,3,2);
	double d,fraction=0.32;
	double pi=3.141592653;
	vec3D n(sqrt(2.0)/2,-sqrt(2.0)/2,0);
	//vec3D n(-sqrt(3.0)/3,-sqrt(3.0)/3,sqrt(3.0)/3);
	//vec3D n(0,0,-1);
	Tpolyhedron* poly_below;
	Tpolyhedron* poly_above;
	poly_below=new Tpolyhedron;
	poly_above=new Tpolyhedron;
	Tvolume_equation_solver volume_solver;
	volume_solver.set_target_polyhedron(&poly);	
	d=volume_solver.calculate_plane_constant_advanced(n,fraction);
	cout<<"First Volume error is "<<volume_solver.calculate_fraction_below(n,d).x-fraction<<endl;
	poly.cut_polyhedron_by_plane(n,d,poly_below,poly_above);
	volume_solver.set_target_polyhedron(poly_below);
	double theta,phy;
	int i_max,j_max;
	double fraction_error,max_fration_error=0;
	cout<<"Numer of Newton fails "<<volume_solver.G_num_Newton_fails()<<endl;
	for (int i = 0; i < 11; i++)
	{
		cout<<"i="<<i<<endl;
		for (int j = 0; j < 11; j++)
		{
			cout<<"j="<<j<<endl;
			theta=0.1*pi*i;phy=0.2*pi*j;
			n.value(sin(theta)*cos(phy),sin(theta)*sin(phy),cos(theta));
			d=volume_solver.calculate_plane_constant_advanced(n,fraction);
			fraction_error=abs(volume_solver.calculate_fraction_below(n,d).x-fraction);
			if (fraction_error>max_fration_error)
			{
				i_max=i;j_max=j;
				max_fration_error=fraction_error;
			}
			//cout<<"Second Volume error is "<<volume_solver.calculate_fraction_below(n,d).x-fraction<<endl;
		}
	}
	cout<<"Maximal fraction error is "<<max_fration_error<<" at "<<i_max<<" "<<j_max<<endl;
	cout<<"Number of calling volume solver: "<<volume_solver.G_num_solving()<<endl;
	cout<<"Numer of Newton fails "<<volume_solver.G_num_Newton_fails()<<endl;
	return;
}
void Tbody_pure_ALE::test_for_MoF_derivative()
{
	TMoF_solver MoF_solver(2,0);
	set_velocity(2);
	update_position(0.1);
	Tpolyhedron* polyhedron_ptr;
	double fraction;
	vec3D centroid;
	fraction=0.3;
	centroid.value(0.1,0.2,0.3);
	vec2D error;
	double error_max_theta=0,error_max_phy=0;
	int i_max_error_theta,i_max_error_phy;
	for (int i = 0; i <_nume; i++)
	{
		polyhedron_ptr=&_grid_old->_cell_list[i]._cell_polyhedron;
		error=MoF_solver.derivative_comparison(polyhedron_ptr,fraction,centroid);
		if (error.x>1e-4 || error.y>1e-4)
		{
			cout<<"Error cell is "<<i<<endl;
			cout<<"Error in derivative to theta is "<<error.x<<endl;
			cout<<"Error in derivative to phy is   "<<error.y<<endl;
		}
		if (error.x>error_max_theta)
		{
			error_max_theta=error.x;
			i_max_error_theta=i;
		}
		if (error.y>error_max_phy)
		{
			error_max_phy=error.y;
			i_max_error_phy=i;
		}
	}
	cout<<"Maximal error for derivative to theta is "<<error_max_theta<<" at cell "<<i_max_error_theta<<endl;
	cout<<"Maximal error for derivative to phy   is "<<error_max_phy<<" at cell "<<i_max_error_phy<<endl;
	return;
}
void Tbody_pure_ALE::test_for_MoF_patch_test()
{
	TMoF_solver MoF_solver(2,0);
	set_velocity(2);
	update_position(0.1);
	Tpolyhedron* polyhedron_ptr;
	double fraction;
	fraction=0.3;
	vec2D error;
	double error_max_theta=0,error_max_phy=0;
	int i_max_error_theta,i_max_error_phy;
	for (int i = 867; i <_nume; i++)
	{
		cout<<i<<endl;
		polyhedron_ptr=&_grid_old->_cell_list[i]._cell_polyhedron;
		error=MoF_solver.patch_test(polyhedron_ptr,fraction);
		//if (error.x>1e-4 || error.y>1e-4)
		//{
		//	cout<<"Error cell is "<<i<<endl;
		//	cout<<"Error in derivative to theta is "<<error.x<<endl;
		//	cout<<"Error in derivative to phy is   "<<error.y<<endl;
		//}
		if (error.x>error_max_theta)
		{
			error_max_theta=error.x;
			i_max_error_theta=i;
		}
		if (error.y>error_max_phy)
		{
			error_max_phy=error.y;
			i_max_error_phy=i;
		}
	}
	cout<<"Maximal error for derivative to theta is "<<error_max_theta<<" at cell "<<i_max_error_theta<<endl;
	cout<<"Maximal error for derivative to phy   is "<<error_max_phy<<" at cell "<<i_max_error_phy<<endl;
	cout<<"Number of inaccurate volume equation solution is "<<MoF_solver.G_volume_solver()->G_num_Newton_fails()<<endl;
	return;
}
void Tbody_pure_ALE::test_for_MoF_patch_test_3_material()
{
	double pi=3.141592654;
	TMoF_solver MoF_solver(3,0);
	set_velocity(2);
	update_position(0.1);
	Tpolyhedron* polyhedron_ptr;
	double fraction[2];
	fraction[0]=0.32;fraction[1]=0.44;
	vec2D angle[2];
	int n_divide=9;
	double centroid_error_max=0;
	int max_i1,max_i2,max_i3,max_i4,max_i;
	ofstream output;
	output.open("inaccurate_d.txt");
	for (int i1 = 0; i1 < n_divide+1; i1++)
	{
		for (int i2= 0; i2 < n_divide+1; i2++)
		{
			for (int i3 = 0; i3 < n_divide+1; i3++)
			{
				for (int i4 = 0; i4 < n_divide+1; i4++)
				{
					cout<<"Cycle: "<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" ";
					angle[0].value(i1*pi/n_divide,i2*2*pi/n_divide);angle[1].value(i3*pi/n_divide,i4*2*pi/n_divide);
					double max_error_i=0;
					int num_inaccurate_before=MoF_solver.G_volume_solver()->G_num_Newton_fails();
					for (int i = 0; i <_nume; i++)
					{
						//cout<<i<<endl;
						polyhedron_ptr=&_grid_old->_cell_list[i]._cell_polyhedron;
						int error_before=MoF_solver.G_volume_solver()->G_num_Newton_fails();
						double centroid_error=MoF_solver.three_material_patch_test(polyhedron_ptr,fraction,angle);
						if (centroid_error>max_error_i)
						{
							max_error_i=centroid_error;
						}
						if (centroid_error>centroid_error_max)
						{
							max_i1=i1;max_i2=i2;max_i3=i3;max_i4=i4;max_i=i;
							centroid_error_max=centroid_error;
						}
						int error_after=MoF_solver.G_volume_solver()->G_num_Newton_fails();
						if (error_before != error_after)
						{
							output<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<i<<endl;
						}
					}
					int num_inaccurate_after=MoF_solver.G_volume_solver()->G_num_Newton_fails();
					cout<<"with maximal error "<<max_error_i<<" and "<<num_inaccurate_after-num_inaccurate_before<<" inaccurate d solution"<<endl;
				}
			}	
		}
	}
	cout<<"The maximal centroid error is "<<centroid_error_max<<" at "<<max_i1<<" "<<max_i2<<" "<<max_i3<<" "<<max_i4<<" "<<max_i<<endl;
	cout<<"Number of calling volume solver: "<<MoF_solver.G_volume_solver()->G_num_solving()<<endl;
	cout<<"Number of inaccurate volume equation solution is "<<MoF_solver.G_volume_solver()->G_num_Newton_fails()<<endl;
	return;
}
