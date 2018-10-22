#include "Tbody_ALE.h"
#include "Tcell_pure_fluid.h"
#include "Tcell_mixed_fluid.h"
#include "Tcell_intersect.h"
#include "Tgrid_smooth.h"
#include "Tgrid_smooth_continum.h"
#include <omp.h>
#include <Windows.h>
#include <iomanip>
#include <cmath>
#include <iostream>
#include "public_function.h"
using namespace std;
Tbody_ALE::Tbody_ALE(int i,int num_material,string body_name):Tbody_Lagrangian(i,num_material,body_name)
{
	_grid_new=_grid_old=NULL;
};
bool Tbody_ALE::corrector_step()
{
	bool negative=false,remapping=false;
	//negative=volume_check("Correct");	
	//if (!negative)
	//{
		Tbody_Lagrangian::corrector_step();
		remapping=is_remapping();
	//}
	if (remapping || negative)
	{
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
void Tbody_ALE::S_remapping_scheme(Sremapping_scheme &remapping_scheme)
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
void Tbody_ALE::S_remesh_scheme(string name)
{
	if (name=="Origin")
	{
		_grid_smooth=new Tgrid_smooth_origin(this);
	}
	else if (name=="Continum_analogy")
	{
		_grid_smooth=new Tgrid_smooth_continum(this);
	}
	else if (name=="OST")
	{
		_grid_smooth=new Tgrid_smooth_OST(this);
	}
	else if (name=="CST")
	{
		vec2D center(0,0);
		_grid_smooth=new Tgrid_smooth_CST(this,center);
	}
	else if (name=="PST")
	{
		vec3D center(0,0,0);
		_grid_smooth=new Tgrid_smooth_PST(this,center);
	}
	else
	{
		cout<<"Error:Invalid remesh scheme name"<<endl;
		system("Pause");
		exit(0);
	}
	return;
}
bool Tbody_ALE::is_remapping()
{
	if (_remapping_shceme.name=="None" || _current_time<_remapping_start_time)
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
	else if (_remapping_shceme.name=="Adaption")
	{
		double ratio_angle,ratio_length;
		for (int i = 0; i < _nume; i++)
		{
			_cellptr_list[i]->calculate_distorted_ratio(ratio_angle,ratio_length);
			if (ratio_angle>=_remapping_shceme.tolerance_angle || ratio_length>=_remapping_shceme.tolerance_length)
			{
				cout<<ratio_angle<<" "<<ratio_length<<endl;
				return true;
			}
		}
		return false;
	}
	else if (_remapping_shceme.name=="DtReduction")
	{
		if (_dt/_remapping_shceme.initialDt<=_remapping_shceme.DtReductionRatio)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	return false;
}
void Tbody_ALE::calculate_gradient()
{
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i]->calculate_gradient_in_cell();
		_grid_old->_cell_list[i]->limit_gradient_in_cell();
	}
}
void Tbody_ALE::remapping_variables()
{
	//----------------------------------------prepare for remapping---------------------------------------------
	vec3D coor_min,coor_max;
	double cell_edge_max;
	ofstream output;
	//output.open("temp.txt");
	//system("Pause");
	//output_mesh("Final");
	cout << 1 << endl;
	_grid_smooth->generate_new_grid();
	//Smooth the old grid to generate the new grid
	cout << 2 << endl;
	reset_new_grid_geometry();
	//Reset the geometry information of the new grid
	cout << 3 << endl;
	update_old_grid_geometry();
	//Update the geometry information of the old grid
	cout << 4 << endl;
	old_grid_reconstruction();
	cout << 5 << endl;
	//Reconstruct the variable gradient and material surface on old grid
	calculate_old_grid_size(coor_min,coor_max,cell_edge_max);
	cout << 6 << endl;
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
	int cell_interval=maxval(_nume/1000,1);
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
	_iremapping<<"------Remapping result---------"<<endl;
	_iremapping<<"The average volume error in this remapping phase is: "<<av_volume_error/_nume<<endl;
	_iremapping<<"The maximal volume error in this remapping phase is: "<<max_volume_error<<endl;
	_iremapping<<"The number of intersection in this remapping phase is: "<<num_intersection<<endl;
	_iremapping<<endl;
	if (remapping_failed)
	{
		system("Pause");
	}
	//Reset the nodal velocity and apply boundary conditions on new grid
	//Calculate the nodal variables
	for (int i=0;i<_nump;i++)
	{
		vec3D velocity_origin=_grid_old->_node_list[i]->_velocity;
		_grid_new->_node_list[i]->reset_velocity_after_remapping(velocity_origin);
		_grid_new->_node_list[i]->calculate_nodal_variables();
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
	return;
}
void Tbody_ALE::reset_new_grid_geometry()
{
	bool negative;
	for (int i = 0; i < _nume; i++)
	{
		_grid_new->_cell_list[i]->calculate_coordinate_range();
		//Reset the cell coordinate range of new grid
		_grid_new->_cell_list[i]->calculate_cell_volume("Full",negative);
		if (negative)
		{
			cout<<"Error: Negative volume in calling reset_new_grid_geometry()"<<endl;
			cout<<"body id= "<<_id<<" "<<"cell id= "<<_grid_new->_cell_list[i]->_id<<endl;
		}
		//Reset cell volume and corner volume of new grid
		_grid_new->_cell_list[i]->update_cell_polyhedron();
		//Reset the cell polyhedron of new grid
		//_grid_new->_cell_list[i]->_cell_centroid=_grid_new->_cell_list[i]._cell_polyhedron.G_centroid();
		//Reset the cell centroid of new grid
	}
	return;
}
void Tbody_ALE::update_old_grid_geometry()
{
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i]->calculate_coordinate_range();
		//Update the cell coordinate range of old grid
		_grid_old->_cell_list[i]->update_cell_polyhedron();
		//Update the cell polyhedron of old grid
		//_grid_old->_cell_list[i]->_cell_centroid=_grid_old->_cell_list[i]._cell_polyhedron.G_centroid();
		_grid_old->_cell_list[i]->_cell_volume=_grid_old->_cell_list[i]->_cell_polyhedron.G_volume();
		//Update the cell centroid of old grid
	}
	return;
}
void Tbody_ALE::old_grid_reconstruction()
{
	//Reconstruct material surface on old grid
	surface_reconstrcution();
	//Reconstruct the gradient of variables
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i]->calculate_gradient_in_cell();
	    _grid_old->_cell_list[i]->limit_gradient_in_cell();
	}
	return;
}
void Tbody_ALE::calculate_old_grid_size(vec3D &coor_min,vec3D &coor_max,double &cell_edge_max)
{
	cell_edge_max=0;
	double l_ele;
	for (int i=0;i<_nume;i++)
	{
		l_ele=_grid_old->_cell_list[i]->calculate_maximal_edge();
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
			if (_grid_old->_node_list[i]->_position.access(j)<coor_min.access(j)) coor_min.value(j,_grid_old->_node_list[i]->_position.access(j));
			if (_grid_old->_node_list[i]->_position.access(j)>coor_max.access(j)) coor_max.value(j,_grid_old->_node_list[i]->_position.access(j));
		}
	}
	return;
}
void Tbody_ALE::overlapping_with_material_grid(Tmaterial_grid* material_grid_ptr,int mid)
{
	//----------------------------------------prepare for remapping---------------------------------------------
	vec3D coor_min,coor_max;
	double cell_edge_max;
	update_old_grid_geometry();
	calculate_old_grid_size(coor_min,coor_max,cell_edge_max);
	//Prepare for the bucket generation
	Tbucket_searching buckets(coor_min,coor_max,cell_edge_max);
	//Put the cells in material grid into the bucket
	bool access;
	for (int i=0;i<_nume;i++)
	{
		int ix,iy,iz;
		//vec3D centroid=material_grid_ptr->_material_tet_list[i].G_centroid();
		vec3D centroid=_grid_old->_cell_list[i]->_cell_polyhedron.G_centroid();
		buckets.search_bucket_ijk(centroid,ix,iy,iz,access);
		if (access)
		{
			cout<<"Error: The material polyhedron centroid is out of body range"<<endl;
			system("Pause");
			exit(0);
		}
		buckets.add_element(i,ix,iy,iz);
	}	
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i]->calculate_coordinate_range();
	}
	//Allocate the material information for grid old and grid new
	for (int i = 0; i < _nume; i++)
	{
		Smaterial material=material_grid_ptr->_material;
		SEOS EOS=material_grid_ptr->_EOS;
		_grid_old->_cell_list[i]->_gausspoint[mid]=generate_material(material,EOS,1); 
		_grid_new->_cell_list[i]->_gausspoint[mid]=generate_material(material,EOS,1); 
		if (material_grid_ptr->_is_peusdo)
		{
			double density_old,internal_energy_old,pressure_old,soundspeed_old;
			density_old=_grid_old->_cell_list[i]->_gausspoint[0]->G_density();
			internal_energy_old=_grid_old->_cell_list[i]->_gausspoint[0]->G_internal_energy();
			pressure_old=_grid_old->_cell_list[i]->_gausspoint[0]->G_pressure();
			soundspeed_old=_grid_old->_cell_list[i]->_gausspoint[0]->G_soundspeed();
			_grid_old->_cell_list[i]->_gausspoint[mid]->S_density(density_old);
			_grid_old->_cell_list[i]->_gausspoint[mid]->S_internal_energy(internal_energy_old);
			_grid_old->_cell_list[i]->_gausspoint[mid]->S_pressure(pressure_old);
			_grid_old->_cell_list[i]->_gausspoint[mid]->S_soundspeed(soundspeed_old);
			_grid_new->_cell_list[i]->_gausspoint[mid]->S_density(density_old);
			_grid_new->_cell_list[i]->_gausspoint[mid]->S_internal_energy(internal_energy_old);
			_grid_new->_cell_list[i]->_gausspoint[mid]->S_pressure(pressure_old);
			_grid_new->_cell_list[i]->_gausspoint[mid]->S_soundspeed(soundspeed_old);			
		}
	}
	//------------------------------remapping begin-------------------------------------------
	Ttetrahedron* material_tet;Tcell_fluid_base* old_cell;
	int nx_min,ny_min,nz_min,nx_max,ny_max,nz_max;
	//double error=0,error_max=0;int n_error=0;
	vec3D vertex_list[4];
	vector<int> error_id;
	Tcell_intersect intersection;
	cout<<"overlapping beginning!"<<endl;
	cout<<"Completing ratio: ";
	int num_material_polyhedron=material_grid_ptr->_nume;
	int cell_interval=num_material_polyhedron/20;
	bool first=true;
	//#pragma omp parallel for schedule(dynamic) private(intersection,nx_min,ny_min,nz_min,nx_max,ny_max,nz_max,old_cell,material_tet,vertex_list)
	for (int i=0;i<num_material_polyhedron;i++)
	{
		//cout<<i<<endl;
		vec3D coor_min_poly,coor_max_poly;
		material_tet=&material_grid_ptr->_material_tet_list[i];
		material_tet->calculate_vertex_range(coor_min_poly,coor_max_poly);
		if (i%cell_interval==0)
		{
			if (first)
			{
				cout<<setw(4)<<100*i/num_material_polyhedron<<"%";
				first=false;
			}
			else
			{
				cout<<'\b'<<'\b'<<'\b'<<'\b'<<'\b';
				cout<<setw(4)<<100*i/num_material_polyhedron<<"%";
			}
		}
		for (int j = 0; j < 4; j++)
		{
			vertex_list[j]=material_tet->G_vertex(j);
		}
		buckets.search_bucket_range(vertex_list,4,nx_min,ny_min,nz_min,nx_max,ny_max,nz_max,access);
		if (access)
		{
			cout<<"Error: The material polyehdron is out of the range of old grid"<<endl;
			system("Pause");
			exit(0);
		}
		intersection.S_material_poly(material_tet);
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
						vec3D coor_min_old=old_cell->_coor_min,coor_max_old=old_cell->_coor_max;
						if (!comparison_cell_coordinate(old_cell->_coor_min,old_cell->_coor_max,coor_min_poly,coor_max_poly))
						{
							intersection.S_old_cell(old_cell);
							intersection.overlapping(mid);
						}
					}
				}
			}
		}
	}
	cout<<endl;
	//------------------Calculate the initial state of material mid and output the result------------
	Tcell_fluid_base* cell_ptr=NULL;
	double material_total_volume=0;
	vec3D material_total_centroid;
	for (int i = 0; i < _nume; i++)
	{
		cell_ptr=_grid_old->_cell_list[i];
		double volume_in_cell=cell_ptr->G_material_fraction(mid);//*cell_ptr->G_cell_polyhedron()->G_volume();
		material_total_volume=material_total_volume+volume_in_cell;
		material_total_centroid=material_total_centroid+cell_ptr->G_material_centroid(mid);
	}	
	material_total_centroid=material_total_centroid/material_total_volume;
	double volume_error=abs(material_total_volume-material_grid_ptr->_total_volume);
	double centroid_error=sqrt((material_total_centroid-material_grid_ptr->_total_centroid).self_multuply());
	if (volume_error>1e-12)
	{
		cout<<"Error in overlapping for material "<<mid<<" in body "<<_id<<endl;
		cout<<"Reference volume is "<<material_grid_ptr->_total_volume<<endl;
		cout<<"Calculated volume is "<<material_total_volume<<endl;
		cout<<"Difference is "<<volume_error<<endl;
	}
	cout<<"------Overlapping result of material "<<mid<<" in body "<<_id<<"---------"<<endl;
	cout<<"Volume of material "<<mid<<" is "<<material_total_volume<<" with error "<<volume_error<<endl;
	cout<<"Centroid of material "<<mid<<" is "<<"("<<material_total_centroid.x<<","<<material_total_centroid.y<<","<<material_total_centroid.z<<")"<<endl;
	for (int i=0;i<_nume;i++)
	{	
		//Set the initial state of material mid in old grid
		_grid_old->_cell_list[i]->reset_material_initial_state(mid);
	}
	return;
}
void Tbody_ALE::input_body(ifstream &input)
{
	Skeyword keyword;
	read_in_keyword_file(input,keyword);
	//Read the body global information from keyword
	_nume=keyword.cell_8_list.size();
	_nump=keyword.node_list.size();
	_endtime=keyword.time_control.endtime;
	_remapping_start_time=-1;
	_CFL=keyword.time_control.CFL;
	_hourglass_option=keyword.hourglass_option;
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
		_grid1._node_list[i]=new Tnode_fluid(id,x,y,z);
		_grid2._node_list[i]=new Tnode_fluid(id,x,y,z);
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	if (keyword.part_list.size()==1)           //Input multi-material by other file
	{
		_grid1._cell_list.resize(_nume);
		_grid2._cell_list.resize(_nume);
		_cellptr_list.resize(_nume);
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
				node_ptr1[j]=_grid1._node_list[IENj-1];
				node_ptr2[j]=_grid2._node_list[IENj-1];
			}
			//Read cell properties (nGauss and thickness)
			part_id=keyword.cell_8_list[i].part_id;
			//Read cell material
			material_id=keyword.part_list[part_id-1].material_id;
			EOS_id=keyword.part_list[part_id-1].EOS_id;
			if (EOS_id==0 || keyword.EOS_list.size()==0)
			{
				cout<<"Error:There must be an EOS in body_ALE"<<endl;
				system("Pause");
				exit(0);
			}
			TMAT_base* mat1;
			TMAT_base* mat2;
			Smaterial this_mat=keyword.material_list[material_id-1];
			SEOS this_EOS=keyword.EOS_list[EOS_id-1];
			mat1=generate_material(this_mat,this_EOS,1);
			mat2=generate_material(this_mat,this_EOS,1);
			//Create a new cell
			if (_num_material==1)
			{
				_grid1._cell_list[i]=new Tcell_pure_fluid(id,node_ptr1,mat1);
				_grid2._cell_list[i]=new Tcell_pure_fluid(id,node_ptr2,mat2);
			}
			else
			{
				_grid1._cell_list[i]=new Tcell_mixed_fluid(id,_num_material,node_ptr1,mat1);
				_grid2._cell_list[i]=new Tcell_mixed_fluid(id,_num_material,node_ptr2,mat2);
			}
		}
	}
	else           //Read the multi-material from k file directly
	{
		_grid1._cell_list.resize(_nume);
		_grid2._cell_list.resize(_nume);
		_cellptr_list.resize(_nume);
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
				node_ptr1[j]=_grid1._node_list[IENj-1];
				node_ptr2[j]=_grid2._node_list[IENj-1];
			}
			//Allocate the materials in this problem
			TMAT_base** mat1;
			TMAT_base** mat2;
			mat1=new TMAT_base*[_num_material];
			mat2=new TMAT_base*[_num_material];
			for (int j = 0; j < keyword.part_list.size(); j++)
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
			}
			//Create a new cell
			if (_num_material == 1)
			{
				part_id = keyword.cell_8_list[i].part_id;
				material_id = keyword.part_list[part_id - 1].material_id;
				if (material_id==2)
				{
					double sss = 1;
				}
				_grid1._cell_list[i] = new Tcell_pure_fluid(id, node_ptr1, mat1[material_id-1]);
				_grid2._cell_list[i] = new Tcell_pure_fluid(id, node_ptr2, mat2[material_id-1]);
			}
			else
			{
				_grid1._cell_list[i] = new Tcell_mixed_fluid(id, _num_material, node_ptr1, mat1);
				_grid2._cell_list[i] = new Tcell_mixed_fluid(id, _num_material, node_ptr2, mat2);
				//Initialize the volume fraction and material centroid
				part_id=keyword.cell_8_list[i].part_id;
				material_id=keyword.part_list[part_id-1].material_id;
				_grid1._cell_list[i]->S_material_fraction(1,material_id-1);
				_grid1._cell_list[i]->S_material_centroid(_grid1._cell_list[i]->G_cell_polyhedron()->G_centroid(),material_id-1);
				_grid2._cell_list[i]->S_material_fraction(1,material_id-1);
				_grid2._cell_list[i]->S_material_centroid(_grid2._cell_list[i]->G_cell_polyhedron()->G_centroid(),material_id-1);
			}
		}
	}
	//-------------------------------------------------------
	//Read boundary condition
	//-------------------------------------------------------
	int node_group_id;   //Boundary condition (load) is applied on node_group_id

	int num_bd_group;    //Number of boundary condition group
	int num_node;        //Number of nodes in the node group
	num_bd_group=keyword.boundary_list.size();
	int bc_type;
	for (int i = 0; i < num_bd_group; i++)
	{
		node_group_id=keyword.boundary_list[i].id;
		num_node=keyword.node_group_list[node_group_id-1].node_id.size();
		bc_type=keyword.boundary_list[i].type;
		if (bc_type==0)
		{
			for (int j = 0; j < 3; j++)
			{
				if (keyword.boundary_list[i].pos[j]==1)
				{
					for (int k = 0; k < num_node; k++)
					{
						_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[j]=1;
						_grid2._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[j]=1;
					}
				}
			}
		}
		else if (bc_type==1)
		{
			for (int j = 0; j < 3; j++)
			{
				if (keyword.boundary_list[i].dofvel[j]==1)
				{
					for (int k = 0; k < num_node; k++)
					{
						_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[j]=2;
						_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[j]=keyword.boundary_list[i].vel[j];
						_grid2._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[j]=2;
						_grid2._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[j]=keyword.boundary_list[i].vel[j];
					}
				}
			}
		}
		else if (bc_type==2 || bc_type==3)
		{
			vec2D center(keyword.boundary_list[i].x,keyword.boundary_list[i].y);
			double v=keyword.boundary_list[i].v;
			for (int k = 0; k < num_node; k++)
			{
				for (int j = 0; j < 3; j++)
				{
					_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[j]=3;
					_grid2._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_type_position[j]=3;
				}
				//vec3D coor=_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_position;
				//vec2D r(coor.x-center.x,coor.y-center.y);
				//r=r/r.get_length();
				//_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[0]=v*r.x;
				//_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[1]=v*r.y;
				//_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[2]=0;
				//_grid2._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[0]=v*r.x;
				//_grid2._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[1]=v*r.y;
				//_grid2._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->_bc_velocity[2]=0;
			}
		}
	}

	//---------------------------------------------------------------------------------
	//Perturbance the node for implosion problem
	//---------------------------------------------------------------------------------
	if (keyword.is_implosion)
	{
		SNodePerturbation node_per=keyword.im_node_perturbation;
		if (node_per.dimension==2 ||node_per.dimension==3)
		{
			//Perturbance the node
			node_group_id=node_per.id;
			vec3D center(node_per.x,node_per.y,0);
			num_node=keyword.node_group_list[node_group_id-1].node_id.size();
			for (int k = 0; k < num_node; k++)
			{
				if (keyword.node_group_list[node_group_id-1].node_id[k]==870)
				{
					double a=0;
				}
				_grid1._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]->perturbance_node_for_implosion(center,node_per.dimension,node_per.ratio,8);
			}
			//Move the node
			//double alpha=calcualte_alpha(keyword.im.d0,keyword.im.N,keyword.im.r_in,keyword.im.r_out);
			//for (int i = 0; i < _nump; i++)
			//{
			//	_grid1._node_list[i]->move_node_for_implosion(center,node_per.dimension,keyword.im.N,keyword.im.d0,keyword.im.r_in,keyword.im.r_out,alpha);
			//}
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

	//_grid_old point to grid1 and _grid_new point to grid2 initially
	_grid_old=&_grid1;_grid_new=&_grid2;
	bool negative;
	for (int i = 0; i < _nume; i++)
	{
		_grid1._cell_list[i]->calculate_cell_volume("Full",negative);
	}
	//update_old_grid_geometry();

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
double Tbody_ALE::calculate_time_step()
{
	double time_step=Tbody_Lagrangian::calculate_time_step();
	if (_remapping_shceme.name=="DtReduction")
	{
		if (_current_step==0)
		{
			_remapping_shceme.initialDt=time_step;
		}
	}
	return time_step;
}
void Tbody_ALE::input_immersed_material()
{
	ifstream input_file,input_material_grid;
	ofstream output;
	string a;
	input_file.open("extra_information.dat");
	output.open("material_grid.dat");
	if (input_file.fail())
	{
		cout<<"Error: No extra information input"<<endl;
		system("Pause");
		exit(0);
	}
	while (true)
	{
		input_file>>a;
		if (a=="*MULTI_MATERIAL_LIST")
		{
			while (!exam_keyword(input_file))
			{
				int body_id,material_id;
				bool is_plot,is_pseudo;
				string file_name;
				input_file>>body_id>>material_id>>file_name>>is_plot>>is_pseudo;
				if (_id==body_id)
				{
					if (material_id==0)
					{
						cout<<"Error: material"<<material_id<<" in *MULTI_MATERIAL must lagrer than zero "<<endl;
						system("Pause");
						exit(0);
					}
					if (material_id>=_num_material)
					{
						cout<<"Error: material "<<material_id<<" in *MULTI_MATERIAL excess in body "<<_id<<endl;
						cout<<"Number of materials in body "<<body_id<<" is "<<_num_material<<endl;
						system("Pause");
						exit(0);
					}
					input_material_grid.open(file_name);
					if (input_material_grid.fail())
					{
						cout<<"Error:File "<<file_name<<" is not existed"<<endl;
						system("Pause");
						exit(0);
					}
					Tmaterial_grid material_grid(material_id,is_pseudo);
					material_grid.input_material_grid(input_material_grid);
					overlapping_with_material_grid(&material_grid,material_id);
					input_material_grid.close();
					if (is_plot)
					{
						//Plot the material grid is needed
						material_grid.plot_material_grid(output);
					}
				}
			}
		}
		else if (a=="*END")
		{
			break;
		}
	}
	reset_final_initial_condition();
	return;
}
//===============================================================================
//Functions for test
//===============================================================================

void Tbody_ALE::set_velocity(int type)
{
	if (type==1)
	{
		for (int i = 0; i < _nump; i++)
		{
			vec3D pos=_grid_old->_node_list[i]->_position;
			//_grid_old->_node_list[i]._velocity.x=0.25-(pos.x-0.5)*(pos.x-0.5);
			//_grid_old->_node_list[i]._velocity.y=0.25-(pos.y-0.5)*(pos.y-0.5);
			_grid_old->_node_list[i]->_velocity.z=pos.x+2*pos.y+3*pos.z;
		}
	}
	else if (type==2)
	{
		vec3D temp1,temp2;
		double l_max,l_min;
		calculate_body_size(temp1,temp2,l_max,l_min);
		for (int i = 0; i < _nump; i++)
		{
			vec3D pos=_grid_old->_node_list[i]->_position;
			_grid_old->_node_list[i]->_velocity.x=random(-0.5,0.5)*l_max;
			_grid_old->_node_list[i]->_velocity.y=random(-0.5,0.5)*l_max;
			_grid_old->_node_list[i]->_velocity.z=random(-0.5,0.5)*l_max;
			if (abs((pos.x-1)*pos.x)<1e-10)
			{
				_grid_old->_node_list[i]->_velocity.x=0;
			}
			if (abs((pos.y-1)*pos.y)<1e-10)
			{
				_grid_old->_node_list[i]->_velocity.y=0;
			}
			if (abs((pos.z-1)*pos.z)<1e-10)
			{
				_grid_old->_node_list[i]->_velocity.z=0;
			}
		}
	}
	return;
}
void Tbody_ALE::update_position(double dt)
{
	for (int i = 0; i < _nump; i++)
	{
		vec3D pos=_grid_old->_node_list[i]->_position;
		vec3D vel=_grid_old->_node_list[i]->_velocity;
		_grid_old->_node_list[i]->_position=pos+vel*dt;
	}
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i]->update_cell_polyhedron();
		//_grid_old->_cell_list[i]->_cell_centroid=_grid_old->_cell_list[i]->_cell_polyhedron.G_centroid();
	}
	return;
}
void Tbody_ALE::set_variables()
{
	for (int i = 0; i < _nume; i++)
	{
		vec3D cell_centroid=_grid_old->_cell_list[i]->_cell_polyhedron.G_centroid();
		vec3D g_ro(5,6,7),g_roe(5,6,7);
		double ro=1+g_ro*cell_centroid;
		double roe=1+g_roe*cell_centroid;
		_grid_old->_cell_list[i]->_gausspoint[0]->S_density(ro);
		_grid_old->_cell_list[i]->_gausspoint[0]->S_internal_energy(roe/ro);
	}
	return;
}
void Tbody_ALE::set_variables(int mid)
{
	int index=mid+1;
	for (int i = 0; i < _nume; i++)
	{
		if (_cellptr_list[i]->G_material_fraction(mid)>0)
		{
			vec3D material_centroid=_cellptr_list[i]->G_material_centroid(mid);
			vec3D g_ro(5,6,7),g_roe(5,6,7);
			g_ro=g_ro*index;g_roe=g_roe*index;
			double ro=1+g_ro*material_centroid;
			double roe=1+g_roe*material_centroid;
			_cellptr_list[i]->_gausspoint[mid]->S_density(ro);
			_cellptr_list[i]->_gausspoint[mid]->S_internal_energy(roe/ro);
		}
		
	}
}
void Tbody_ALE::iteration_test()
{
	for (int i = 0; i < _nume; i++)
	{
		vec3D coordinate_std(random(-1,1),random(-1,1),random(-1,1));
		vec3D coordinate_physical=_grid_old->_cell_list[i]->shape_function_interpolate(coordinate_std,"Position");
		vec3D result=_grid_old->_cell_list[i]->calculate_std_coordinate(coordinate_physical);
		double error=sqrt((result-coordinate_std).self_multuply());
		if (error>1e-10)
		{
			cout<<i<<" "<<error<<endl;
		}
	}
	return;
}
void Tbody_ALE::set_displacemant_at_conode(double dl)
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
void Tbody_ALE::plot_grid_new()
{
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]=_grid_new->_cell_list[i];
	}
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]=_grid_new->_node_list[i];
	}
	output_mesh("Final");
	return;
}
void Tbody_ALE::test_for_fraction_below()
{
	Tpolyhedron* poly0;
	Tpolyhedron* poly1;
	Tpolyhedron* poly2;
	Tpolyhedron* poly3;
	Tpolyhedron* poly4;
	poly0=&_cellptr_list[0]->_cell_polyhedron;
	poly0->calculate_piece_topology();
	double pi=3.141592654;
	double theta=0.2*pi,phy=0.2*pi;
	vec3D n;
	ofstream output;
	output.open("polyhedron.dat");
	poly0->plot_this_polyehdron(output);
	n.value(sin(theta)*cos(phy),sin(theta)*sin(phy),cos(theta));
	double d_min,d_max,d;
	poly0->calculate_altitude_range(n,d_min,d_max);
	d=0.5*(d_min+d_max);
	poly1=new Tpolyhedron;poly2=new Tpolyhedron;
	poly0->cut_polyhedron_by_plane(n,d,poly1,poly2);
	double dl=0.04;
	for (int i = 0; i < 8; i++)
	{
		poly1->move_polyhedron(n,-i*dl);
		poly2->move_polyhedron(n,i*dl);
		poly1->plot_this_polyehdron(output,1);
		poly2->plot_this_polyehdron(output,2);
	}
	theta=0.2*pi,phy=0.7*pi;
	n.value(sin(theta)*cos(phy),sin(theta)*sin(phy),cos(theta));
	poly1->calculate_altitude_range(n,d_min,d_max);
	d=0.5*(d_min+d_max);
	poly3=new Tpolyhedron;poly4=new Tpolyhedron;
	poly1->calculate_piece_topology();
	poly1->cut_polyhedron_by_plane(n,d,poly3,poly4);
	for (int i = 0; i < 8; i++)
	{
		poly3->move_polyhedron(n,-i*dl);
		poly4->move_polyhedron(n,i*dl);
		poly3->plot_this_polyehdron(output,1);
		poly4->plot_this_polyehdron(output,3);
	}
	//ofstream output;
	//output.open("polyhedron.dat");
	//poly0->plot_this_polyehdron(output);
	//poly1->plot_this_polyehdron(output,1);
	//poly2->plot_this_polyehdron(output,2);
	//poly3->plot_this_polyehdron(output,1);
	//poly4->plot_this_polyehdron(output,3);
	//Tvolume_equation_solver volume_solver1;
	//Tvolume_equation_solver volume_solver2;
	//set_velocity(2);
	//update_position(0.1);
	//Tpolyhedron* polyhedron_ptr;
	//double theta,phy,d,d_min,d_max,ratio;
	//vec3D n;
	//vec3D result1,result2;
	//double fraction_error,derivative_error;
	//double fraction_error_max=0,derivative_error_max=0;
	//theta=phy=1;ratio=0.1;
	//int n_divide=211;
	//double pi=3.141592654,dtheta=pi/n_divide,dphy=2*pi/n_divide,dratio=1.0/n_divide;
	//ios::scientific;
	//ofstream output;
	//output.open("polyhedron.dat");
	//for (int i = 0; i <_nume; i++)
	//{
	//	polyhedron_ptr=&_grid_old->_cell_list[i]->_cell_polyhedron;
	//	polyhedron_ptr->plot_this_polyehdron(output);
	//	volume_solver1.set_target_polyhedron(polyhedron_ptr);
	//	volume_solver2.set_target_polyhedron(polyhedron_ptr);
	//	for (int j = 0; j < n_divide; j++)
	//	{
	//		for (int k = 0; k < n_divide; k++)
	//		{
	//			theta=j*dtheta;phy=k*dphy;
	//			n.value(sin(theta)*cos(phy),sin(theta)*sin(phy),cos(theta));
	//			polyhedron_ptr->calculate_altitude_range(n,d_min,d_max);
	//			for (int l = 1; l < n_divide-1; l++)
	//			{
	//				ratio=l*dratio;
	//				d=d_min+ratio*(d_max-d_min);
	//				result1=volume_solver1.calculate_fraction_below(n,d);
	//				result2=volume_solver2.calculate_fraction_below(n,d);
	//				fraction_error=abs(result1.x-result2.x);
	//				derivative_error=abs(result1.y-result2.y);
	//				if (fraction_error>1e-5)
	//				{
	//					output<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
	//					output<<result1.x<<" "<<result2.x<<endl;
	//					output<<fraction_error<<endl;
	//					cout<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
	//					cout<<result1.x<<" "<<result2.x<<endl;
	//					cout<<fraction_error<<endl;
	//					output<<endl;
	//				}
	//				fraction_error_max=maxval(fraction_error_max,fraction_error);
	//				derivative_error_max=maxval(derivative_error_max,derivative_error);
	//			}
	//		}
	//	}
	//}
	//cout<<"The maximal fraction error is "<<fraction_error_max<<endl;
	//cout<<"The maximal derivative error is "<<derivative_error_max<<endl;
	return;
}
void Tbody_ALE::test_for_volume_solver()
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
void Tbody_ALE::test_for_MoF_derivative()
{
	TMoF_solver MoF_solver(2,0,3);
	set_velocity(2);
	update_position(0.1);
	Tpolyhedron* polyhedron_ptr;
	double fraction;
	vec3D centroid;
	fraction=0.32;
	centroid.value(0.2,0.2,0.2);
	vec2D error;
	double error_max_theta=0,error_max_phy=0;
	int i_max_error_theta,i_max_error_phy;
	for (int i = 0; i <_nume; i++)
	{
		polyhedron_ptr=&_grid_old->_cell_list[i]->_cell_polyhedron;
		error=MoF_solver.derivative_comparison(polyhedron_ptr,fraction,centroid);
		if (abs(error.x)>error_max_theta)
		{
			error_max_theta=abs(error.x);
			i_max_error_theta=i;
		}
		if (abs(error.y)>error_max_phy)
		{
			error_max_phy=abs(error.y);
			i_max_error_phy=i;
		}
	}
	cout<<"Maximal error for derivative to theta is "<<error_max_theta<<" at cell "<<i_max_error_theta<<endl;
	cout<<"Maximal error for derivative to phy   is "<<error_max_phy<<" at cell "<<i_max_error_phy<<endl;
	return;
}
void Tbody_ALE::test_for_MoF_patch_test()
{
	TMoF_solver MoF_solver(2,0,3);
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
		polyhedron_ptr=&_grid_old->_cell_list[i]->_cell_polyhedron;
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
void Tbody_ALE::test_for_MoF_patch_test_3_material()
{
	double pi=3.141592654;
	TMoF_solver MoF_solver(3,0,3);
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
						polyhedron_ptr=&_grid_old->_cell_list[i]->_cell_polyhedron;
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
void Tbody_ALE::test_for_MoF()
{
	ofstream output;
	output.open("MoFInput.txt");
	output.precision(20);
	for (int i = 0; i < _nume; i++)
	{
		double f0=_cellptr_list[i]->G_material_fraction(0);
		vec3D centroid=_cellptr_list[i]->G_material_centroid(0);
		output<<f0<<" "<<centroid.x<<" "<<centroid.y<<" "<<centroid.z<<endl;
	}
	surface_reconstrcution();
	plot_reconstructed_material();
	return;
}
void Tbody_ALE::set_new_grid(int type,double dt)
{
	if (type==1)
		{
			for (int i = 0; i < _nump; i++)
			{
				vec3D pos=_grid_new->_node_list[i]->_position;
				_grid_new->_node_list[i]->_velocity.x=0.25-(pos.x-0.5)*(pos.x-0.5);
				_grid_new->_node_list[i]->_velocity.y=0.25-(pos.y-0.5)*(pos.y-0.5);
				_grid_new->_node_list[i]->_velocity.z=pos.x+2*pos.y+3*pos.z;
		}
	}
	else if (type==2)
	{
		vec3D temp1,temp2;
		double l_max,l_min;
		calculate_body_size(temp1,temp2,l_max,l_min);
		for (int i = 0; i < _nump; i++)
		{
			vec3D pos=_grid_new->_node_list[i]->_position;
			_grid_new->_node_list[i]->_velocity.x=random(-0.5,0.5)*l_max;
			_grid_new->_node_list[i]->_velocity.y=random(-0.5,0.5)*l_max;
			_grid_new->_node_list[i]->_velocity.z=random(-0.5,0.5)*l_max;
			if (abs((pos.x-1)*pos.x)<1e-10)
			{
				_grid_new->_node_list[i]->_velocity.x=0;
			}
			if (abs((pos.y-1)*pos.y)<1e-10)
			{
				_grid_new->_node_list[i]->_velocity.y=0;
			}
			if (abs((pos.z-1)*pos.z)<1e-10)
			{
				_grid_new->_node_list[i]->_velocity.z=0;
			}
		}
	}
	for (int i = 0; i < _nump; i++)
	{
		vec3D pos=_grid_new->_node_list[i]->_position;
		vec3D vel=_grid_new->_node_list[i]->_velocity;
		_grid_new->_node_list[i]->_position=pos+vel*dt;
	}
	for (int i = 0; i < _nume; i++)
	{
		_grid_new->_cell_list[i]->update_cell_polyhedron();
	}
	return;
}
void Tbody_ALE::test_for_MM_remapping()
{
	surface_reconstrcution();
	plot_reconstructed_material();
	for (int i = 0; i < _num_material; i++)
	{
		set_variables(i);
	}
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i]->calculate_gradient_in_cell();
	}
	double error_g_ro[2],error_g_roe[2];
	for (int i = 0; i < 2; i++)
	{
		error_g_ro[i]=error_g_roe[i]=0;
	}
	vec3D g[2];
	g[0].value(5,6,7);g[1].value(10,12,14);
	for (int i = 0; i <  _nume; i++)
	{
		Tcell_fluid_base* cell_ptr=_grid_old->_cell_list[i];
		for (int j = 0; j < _num_material; j++)
		{
			if (cell_ptr->G_material_fraction(j)>0)
			{
				vec3D g_ro=cell_ptr->G_g_ro(j);
				error_g_ro[j]=error_g_ro[j]+(cell_ptr->G_g_ro(j)-g[j]).self_multuply();
				error_g_roe[j]=error_g_roe[j]+(cell_ptr->G_g_roe(j)-g[j]).self_multuply();
			}
		}
	}
	cout<<"initial error"<<endl;
	for (int i = 0; i < _num_material; i++)
	{
		cout<<"Gradient error of ro in material "<<i<<" is "<<error_g_ro[i]<<endl;
		cout<<"Gradient error of roe in material "<<i<<" is "<<error_g_roe[i]<<endl;
	}
	system("Pause");
	set_new_grid(2,0.1);
	remapping_variables();
	surface_reconstrcution();
	plot_reconstructed_material();
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i]->calculate_gradient_in_cell();
	}
	for (int i = 0; i < 2; i++)
	{
		error_g_ro[i]=error_g_roe[i]=0;
	}
	for (int i = 0; i <  _nume; i++)
	{
		Tcell_fluid_base* cell_ptr=_grid_old->_cell_list[i];
		for (int j = 0; j < _num_material; j++)
		{
			if (cell_ptr->G_material_fraction(j)>0)
			{
				vec3D g_ro=cell_ptr->G_g_ro(j);
				error_g_ro[j]=error_g_ro[j]+(cell_ptr->G_g_ro(j)-g[j]).self_multuply();
				error_g_roe[j]=error_g_roe[j]+(cell_ptr->G_g_roe(j)-g[j]).self_multuply();
			}
		}
	}
	cout<<"first remapping error"<<endl;
	for (int i = 0; i < _num_material; i++)
	{
		cout<<"Gradient error of ro in material "<<i<<" is "<<error_g_ro[i]<<endl;
		cout<<"Gradient error of roe in material "<<i<<" is "<<error_g_roe[i]<<endl;
	}
	system("Pause");
	set_new_grid(2,0.1);
	remapping_variables();
	surface_reconstrcution();
	plot_reconstructed_material();
	for (int i = 0; i < _nume; i++)
	{
		_grid_old->_cell_list[i]->calculate_gradient_in_cell();
	}
	for (int i = 0; i < 2; i++)
	{
		error_g_ro[i]=error_g_roe[i]=0;
	}
	for (int i = 0; i <  _nume; i++)
	{
		Tcell_fluid_base* cell_ptr=_grid_old->_cell_list[i];
		for (int j = 0; j < _num_material; j++)
		{
			if (cell_ptr->G_material_fraction(j)>0)
			{
				vec3D g_ro=cell_ptr->G_g_ro(j);
				error_g_ro[j]=error_g_ro[j]+(cell_ptr->G_g_ro(j)-g[j]).self_multuply();
				error_g_roe[j]=error_g_roe[j]+(cell_ptr->G_g_roe(j)-g[j]).self_multuply();
			}
		}
	}
	cout<<"second remapping error"<<endl;
	for (int i = 0; i < _num_material; i++)
	{
		cout<<"Gradient error of ro in material "<<i<<" is "<<error_g_ro[i]<<endl;
		cout<<"Gradient error of roe in material "<<i<<" is "<<error_g_roe[i]<<endl;
	}
	return;
}
