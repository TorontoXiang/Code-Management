#include "Tbody_pure_Lagrangian.h"
#include "public_function.h"
#include <iomanip>
using namespace std;
void Tbody_pure_Lagrangian::calculate_force_in_body()
{
	//Calcualte the corner force for every cell
	calculate_corner_force(_dt);
	//Calculate the final nodal force
	assemble_nodal_force();
	return;
}
double Tbody_pure_Lagrangian::calculate_time_step()
{
	double min_time_step=100000;
	for (int i = 0; i < _nume; i++)
	{
		min_time_step=minval(min_time_step,_cellptr_list[i]->calculate_time_step());
	}
	if (_current_step==0)
	{
		return _CFL*min_time_step;
	}
	else
	{
		return minval(_CFL*min_time_step,1.05*_dt);
	}
}
void Tbody_pure_Lagrangian::predictor_step()
{
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->predict_position(_dt);
	}
	calculate_force_in_body();
	return;
}
void Tbody_pure_Lagrangian::calculate_final_acceleration()
{
	for (int i = 0; i < _nump; i++)
	{
		//Calculate the final nodal acceleration
		_nodeptr_list[i]->calculate_acceleration();
	}
	return;
}
bool Tbody_pure_Lagrangian::corrector_step()
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
	return false;
}
void Tbody_pure_Lagrangian::calculate_nodal_inertance()
{
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->assemble_corner_inertance();
	}
	return;
}
void Tbody_pure_Lagrangian::calculate_corner_force(double dt)
{
	Tviscosity_base* Vis_ptr=_Arti_Vis;
	#pragma omp parallel for firstprivate(Vis_ptr)
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->calculate_corner_force(Vis_ptr,_dt);
	}
}
void Tbody_pure_Lagrangian::assemble_nodal_force()
{
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->assemble_corner_force();
	}
	return;
}
void Tbody_pure_Lagrangian::set_grid_topology()
{
	//Calculate cells connected to a node
	for (int i=0;i<_nume;i++)
	{
		_cellptr_list[i]->add_to_node();
	}
	//Calculate _adjacent[] and _around
	for (int i=0;i<_nume;i++)
	{
		_cellptr_list[i]->set_cell_connection();
	}
	return;
}
double Tbody_pure_Lagrangian::calculate_total_internal_energy()
{
	double total_internal_energy=0;
	for (int i = 0; i < _nume; i++)
	{
		total_internal_energy=total_internal_energy+_cellptr_list[i]->_cell_mass*_cellptr_list[i]->_gausspoint[0]->G_internal_energy();
	}
	return total_internal_energy;
}
double Tbody_pure_Lagrangian::calculate_maximal_pressure()
{
	double max_pressure=0;
	for (int i = 0; i < _nume; i++)
	{
		double cell_pressuer=_cellptr_list[i]->_gausspoint[0]->G_pressure();
		max_pressure=maxval(max_pressure,cell_pressuer);
	}
	return max_pressure;
}
void Tbody_pure_Lagrangian::calculate_body_size(vec3D &coor_min,vec3D &coor_max,double &cell_edge_max,double &cell_edge_min)
{
	cell_edge_max=0,cell_edge_min=1000000;
	double l_ele;
	vec3D eps(1e-5,1e-5,1e-5);
	for (int i=0;i<_nume;i++)
	{
		l_ele=_cellptr_list[i]->calculate_characteristic_length();
		if (l_ele>cell_edge_max)
		{
			cell_edge_max=l_ele;
		}
		if (l_ele<cell_edge_min)
		{
			cell_edge_min=l_ele;
		}
	}
	coor_min.value(1e10,1e10,1e10);coor_max.value(-1e10,-1e10,-1e10);
	for (int i=0;i<_nump;i++)
	{
		for (int j = 1; j < 4; j++)
		{
			if (_nodeptr_list[i]->_position.access(j)<coor_min.access(j)) coor_min.value(j,_nodeptr_list[i]->_position.access(j));
			if (_nodeptr_list[i]->_position.access(j)>coor_max.access(j)) coor_max.value(j,_nodeptr_list[i]->_position.access(j));
		}
	}
	//To avoid a node is exactly located on the boundary
	coor_min=coor_min-eps;coor_max=coor_max+eps;
	return;
}
//--------------------------------------------------------------------------
//Output functions
//--------------------------------------------------------------------------
void Tbody_pure_Lagrangian::output_tecplot(ofstream &output,double ratio)
{
	output<<"TITLE = \"Mesh for Brick fluid Element\""<<endl;
	output<<"VARIABLES = \"X\",\"Y\",\"Z\",\"density\",\"pressure\",\"vx\",\"vy\",\"vz\",\"mass\""<<endl;
	output<<"ZONE F=FEPOINT,N="<<_nump<<","<<"E="<<_nume<<","<<"ET=BRICK"<<endl;
	vec3D dis;
	for (int i = 0; i < _nump; i++)
	{
		ratio=1;
		dis=_nodeptr_list[i]->calculate_displacement_amplify(ratio);
		output<<dis.x<<" "<<dis.y<<" "<<dis.z<<" "<<_nodeptr_list[i]->_density<<" "<<_nodeptr_list[i]->_pressure
			  <<" "<<_nodeptr_list[i]->_velocity.x<<" "<<_nodeptr_list[i]->_velocity.y<<" "<<_nodeptr_list[i]->_velocity.z
			  <<" "<<_nodeptr_list[i]->_mass<<endl;
	}
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			output<<_cellptr_list[i]->_node_ptr[j]->_id<<" ";
		}
		output<<endl;
	}
	return;
}
void Tbody_pure_Lagrangian::output_curve()
{
	Tbody_base::output_curve();
	int num_curve=_curve_out_list.size();
	int type,id;
	Tnode_fluid* node_ptr;
	Tcell_pure_Lagrangian* cell_ptr;
	for (int i = 0; i < num_curve; i++)
	{
		type=_curve_out_list[i]->_type;
		id=_curve_out_list[i]->_id-1;
		_curve_out_list[i]->output.precision(10);
		if (type==0)
		{
			node_ptr=_nodeptr_list[id];
			vec3D dis=node_ptr->calculate_displacement();
			vec3D vel=node_ptr->_velocity;
			_curve_out_list[i]->output<<setw(20)<<_current_time<<" ";
			_curve_out_list[i]->output<<setw(20)<<vel.x<<setw(20)<<vel.y<<setw(20)<<vel.z<<endl;
		}
		else if (type==1)
		{
			cell_ptr=_cellptr_list[id];
			_curve_out_list[i]->output<<setw(20)<<_current_time<<" ";
			_curve_out_list[i]->output<<setw(20)<<cell_ptr->_gausspoint[0]->G_pressure()<<endl;
		}
	}
	return;
}
void Tbody_pure_Lagrangian::output_cell_monitor()
{
	int num_monitor=_cell_monitor_list.size();
	int cell_id;
	for (int i = 0; i < num_monitor; i++)
	{
		cell_id=_cell_monitor_list[i]->_cell_id-1;
		//_cell_monitor_list[i]->output.precision(10);
		//_cellptr_list[cell_id]->output_cell_monitor(_cell_monitor_list[i]->output);
	}
	return;
}
void Tbody_pure_Lagrangian::input_body(ifstream& input)
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
	_grid._node_list.resize(_nump);
	_nodeptr_list.resize(_nump);
	for (int i = 0; i < _nump; i++)
	{
		id=keyword.node_list[i].id;
		x=keyword.node_list[i].x;
		y=keyword.node_list[i].y;
		z=keyword.node_list[i].z;
		Tnode_fluid new_node(id,x,y,z);
		_grid._node_list[i]=new_node;
		_nodeptr_list[i]=&_grid._node_list[i];
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	_grid._cell_list.resize(_nume);
	_cellptr_list.resize(_nume);
	int nGauss;
	int part_id,material_id,EOS_id;
	Tnode_fluid* node_ptr[8];
	for (int i = 0; i < _nume; i++)
	{
		//Read cell id and IEN
		id=keyword.cell_8_list[i].cell_id;
		for (int j = 0; j < 8; j++)
		{
			int IENj=keyword.cell_8_list[i].IEN[j];
			node_ptr[j]=_nodeptr_list[IENj-1];
		}
		//Read cell properties (nGauss and thickness)
		part_id=keyword.cell_8_list[i].part_id;
		nGauss=1;
		//Read cell material
		material_id=keyword.part_list[part_id-1].material_id;
		EOS_id=keyword.part_list[part_id-1].EOS_id;
		if (EOS_id==0)
		{
			cout<<"Error:There must be an EOS in body_pure_Lagrangian"<<endl;
			system("Pause");
			exit(0);
		}
		TMAT_base** mat;
		mat=new TMAT_base*[nGauss];
		for (int k = 0; k < nGauss; k++)
		{
			Smaterial this_mat=keyword.material_list[material_id-1];
			SEOS this_EOS=keyword.EOS_list[EOS_id-1];
			mat[k]=generate_material(this_mat,this_EOS,1);
		}
		//Create a new cell
		Tcell_pure_Lagrangian new_cell(id,nGauss,node_ptr,mat);
		_grid._cell_list[i]=new_cell;
		_cellptr_list[i]=&_grid._cell_list[i];
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
					_grid._node_list[keyword.node_group_list[node_group_id-1].node_id[k]-1]._bc_type_position[j]=1;
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
	
	set_grid_topology();
	//Get the cell connection

	//Initialize the nodal variables
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->calculate_nodal_variables();
	}

	return;
}