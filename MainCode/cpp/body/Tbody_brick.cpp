#include "Tbody_brick.h"
#include "public_function.h"
#include <string>
#include <iomanip>
#include "TMAT_base.h"
void Tbody_brick::calculate_force_in_body()
{
	//Apply the external force
	apply_external_force();
	//Calcualte the corner force for every cell
	calculate_corner_force();
	//Calculate the final nodal force
	assemble_nodal_force();
	return;
}
void Tbody_brick::assemble_nodal_force()
{
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->assemble_corner_force();
	}
	return;
}
bool Tbody_brick::corrector_step()
{
	for (int i = 0; i < _nume; i++)
	{
		//_cellptr_list[i]->update_strain_energy(_dt);
	}
	for (int i = 0; i < _nump; i++)
	{
		//Calculate the velocity at the next time step
		_nodeptr_list[i]->update_vel(_dt);
		//Calculate the position at the next time point
		_nodeptr_list[i]->update_pos(_dt);
	}
	//Update the material state and calculate the body force for shell element
	calculate_force_in_body();
	return false;
}
void Tbody_brick::calculate_corner_force()
{
	for (int i = 0; i < _nume; i++)
	{
		if (i == 200)
		{
			double s = 1;
		}
		_cellptr_list[i]->calculate_corner_force(_Arti_Vis, _dt);
	}
	return;
}
void Tbody_brick::calculate_final_acceleration()
{
	for (int i = 0; i < _nump; i++)
	{
		//Calculate the final nodal acceleration
		_nodeptr_list[i]->calculate_acceleration();
	}
	return;
}
void Tbody_brick::apply_external_force()
{
	int num_load_group = _external_load_list.size();
	int num_node;
	int direction;
	double magnitude;
	//Apply concentrate load
	for (int i = 0; i < num_load_group; i++)
	{
		direction = _external_load_list[i].direction;
		magnitude = _external_load_list[i].magnitude;
		num_node = _external_load_list[i].num_node;
		for (int j = 0; j < num_node; j++)
		{
			int id = _external_load_list[i].node_id[j];
			_nodeptr_list[id - 1]->apply_external_force(direction, magnitude);
		}
	}
	//Apply gravity
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->apply_gravity(_gravity);
	}
}
double Tbody_brick::calculate_time_step()
{
	double min_time_step = 100000;
	for (int i = 0; i < _nume; i++)
	{
		min_time_step = minval(min_time_step, _cellptr_list[i]->calculate_time_step());
	}
	return _CFL * min_time_step;
}
void Tbody_brick::start_up_step()
{
	calculate_force_in_body();
	return;
}
void Tbody_brick::calculate_nodal_inertance()
{
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->assemble_corner_inertance();
	}
	return;
}
double Tbody_brick::calculate_total_internal_energy()
{
	double total_internal_energy = 0;
	for (int i = 0; i < _nume; i++)
	{
		total_internal_energy = total_internal_energy + _cellptr_list[i]->_strain_energy;
	}
	return total_internal_energy;
}
void Tbody_brick::calculate_body_size(vec3D &coor_min, vec3D &coor_max, double &cell_edge_max, double &cell_edge_min)
{
	cell_edge_max = 0, cell_edge_min = 1000000;
	double l_ele;
	for (int i = 0; i < _nume; i++)
	{
		l_ele = _cellptr_list[i]->calculate_characteristic_length();
		if (l_ele > cell_edge_max)
		{
			cell_edge_max = l_ele;
		}
		if (l_ele < cell_edge_min)
		{
			cell_edge_min = l_ele;
		}
	}
	coor_min.value(1e10, 1e10, 1e10); coor_max.value(-1e10, -1e10, -1e10);
	for (int i = 0; i < _nump; i++)
	{
		for (int j = 1; j < 4; j++)
		{
			if (_nodeptr_list[i]->_position.access(j) < coor_min.access(j)) coor_min.value(j, _nodeptr_list[i]->_position.access(j));
			if (_nodeptr_list[i]->_position.access(j) > coor_max.access(j)) coor_max.value(j, _nodeptr_list[i]->_position.access(j));
		}
	}
	return;
}
//------------------------------------------------------------------------
//Output functuions
//------------------------------------------------------------------------
void Tbody_brick::output_curve()
{
	Tbody_base::output_curve();
	int num_curve = _curve_out_list.size();
	int type, id;
	Tnode* node_ptr;
	Tcell_brick* cell_ptr;
	for (int i = 0; i < num_curve; i++)
	{
		type = _curve_out_list[i]->_type;
		id = _curve_out_list[i]->_id - 1;
		_curve_out_list[i]->output.precision(10);
		if (type == 0)
		{
			node_ptr = _nodeptr_list[id];
			vec3D dis = node_ptr->calculate_displacement();
			_curve_out_list[i]->output << setw(20) << _current_time << " ";
			_curve_out_list[i]->output << setw(20) << dis.x << setw(20) << dis.y << setw(20) << dis.z << endl;
		}
		else if (type == 1)
		{
			cell_ptr = _cellptr_list[id];
			double epeff_max = 0;
			for (int j = 0; j < 4; j++)
			{
				if (cell_ptr->_gausspoint[j]->G_epeff() > epeff_max)
				{
					epeff_max = cell_ptr->_gausspoint[j]->G_epeff();
				}
			}
			_curve_out_list[i]->output << setw(20) << _current_time << " ";
			_curve_out_list[i]->output << setw(20) << epeff_max << endl;
		}
	}
	return;
}
double Tbody_brick::calculate_maximal_effective_plastic_strain()
{
	double epeff_max = 0;
	Tcell_brick* cell_ptr;
	for (int i = 0; i < _nume; i++)
	{
		cell_ptr = _cellptr_list[i];
		double epeffi = cell_ptr->_gausspoint[0]->G_epeff();
		if (epeffi > epeff_max)
		{
			epeff_max = epeffi;
		}
	}
	return epeff_max;
}
void Tbody_brick::output_tecplot(ofstream& output, double ratio)
{
	output << "TITLE = \"Mesh for Shell Element\"" << endl;
	output << "VARIABLES = \"X\",\"Y\",\"Z\",\"EQS\"" << endl;
	output << "ZONE T= \"Time="<<_current_time<<"\", F=FEPOINT,N=" << _nump << "," << "E=" << _nume << "," << "ET=BRICK" << endl;
	vec3D dis;
	Snode_stress *node_stress;
	node_stress = new Snode_stress[_nump];
	for (int i = 0; i < _nump; i++)
	{
		node_stress[i].volume = node_stress[i].eqs = 0;
	}
	for (int i = 0; i < _nume; i++)
	{
		double volume = _cellptr_list[i]->_volume*0.125;
		double eqs = _cellptr_list[i]->_gausspoint[0]->calculate_equivalent_stress();
		for (int j = 0; j < 8; j++)
		{
			int node_id = _cellptr_list[i]->_node_ptr[j]->_id - 1;
			node_stress[node_id].eqs = node_stress[node_id].eqs + eqs * volume;
			node_stress[node_id].volume = node_stress[node_id].volume + volume;
		}
	}
	for (int i = 0; i < _nump; i++)
	{
		node_stress[i].eqs = node_stress[i].eqs / node_stress[i].volume;
	}
	for (int i = 0; i < _nump; i++)
	{
		dis = _nodeptr_list[i]->calculate_displacement_amplify(ratio);
		output << dis.x << " " << dis.y << " " << dis.z << " " << node_stress[i].eqs << endl;
	}
	for (int i = 0; i < _nume; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			output << _cellptr_list[i]->_node_ptr[j]->_id << " ";
		}
		output << endl;
	}
	return;
}
//------------------------------------------------------------------------------------
//Input function
//------------------------------------------------------------------------------------
void Tbody_brick::input_body(ifstream& input)
{
	Skeyword keyword;
	read_in_keyword_file(input, keyword);
	//Read the body global information from keyword
	_nume = keyword.cell_8_list.size();
	_nump = keyword.node_list.size();
	_endtime = keyword.time_control.endtime;
	_CFL = keyword.time_control.CFL;
	//----------------------------------------------------
	//Read node list
	//----------------------------------------------------
	int id;
	double x, y, z;
	_grid._node_list.resize(_nump);
	_nodeptr_list.resize(_nump);
	for (int i = 0; i < _nump; i++)
	{
		id = keyword.node_list[i].id;
		x = keyword.node_list[i].x;
		y = keyword.node_list[i].y;
		z = keyword.node_list[i].z;
		Tnode new_node(id, x, y, z);
		_grid._node_list[i] = new_node;
		_nodeptr_list[i] = &_grid._node_list[i];
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	_grid._cell_list.resize(_nume);
	_cellptr_list.resize(_nume);
	int nGauss=1;
	int part_id, material_id, EOS_id;
	Tnode* node_ptr[8];
	for (int i = 0; i < _nume; i++)
	{
		//Read cell id and IEN
		id = keyword.cell_8_list[i].cell_id;
		for (int j = 0; j < 8; j++)
		{
			int IENj = keyword.cell_8_list[i].IEN[j];
			node_ptr[j] = _nodeptr_list[IENj - 1];
		}
		//Read cell's material properties
		part_id = keyword.cell_8_list[i].part_id;
		material_id = keyword.part_list[part_id - 1].material_id;
		EOS_id = keyword.part_list[part_id - 1].EOS_id;
		Smaterial this_mat = keyword.material_list[material_id - 1];
		SEOS this_EOS;
		if (EOS_id!=0)
		{
			this_EOS = keyword.EOS_list[EOS_id - 1];
		}
		TMAT_base** mat;
		mat = new TMAT_base*[nGauss];
		for (int k = 0; k < nGauss; k++)
		{
			mat[k] = generate_material(this_mat, this_EOS, 0);
		}
		//Create a new cell
		Tcell_brick new_cell(id, nGauss, node_ptr, mat);
		_grid._cell_list[i] = new_cell;
		_cellptr_list[i] = &_grid._cell_list[i];
	}
	//-------------------------------------------------------
	//Read boundary condition
	//-------------------------------------------------------
	int node_group_id;   //Boundary condition (load) is applied on node_group_id

	int num_bd_group;    //Number of boundary condition group
	int num_node;        //Number of nodes in the node group
	num_bd_group = keyword.boundary_list.size();
	for (int i = 0; i < num_bd_group; i++)
	{
		node_group_id = keyword.boundary_list[i].id;
		num_node = keyword.node_group_list[node_group_id - 1].node_id.size();
		for (int j = 0; j < 3; j++)
		{
			if (keyword.boundary_list[i].pos[j] == 1)
			{
				for (int k = 0; k < num_node; k++)
				{
					_grid._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]._bc_type_position[j] = 1;
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
	num_load_group = keyword.load_list.size();
	for (int i = 0; i < num_load_group; i++)
	{
		node_group_id = keyword.load_list[i].id;
		new_load.node_id = keyword.node_group_list[node_group_id - 1].node_id;
		new_load.num_node = new_load.node_id.size();
		direction = keyword.load_list[i].dof;
		new_load.direction = direction;
		curve_id = keyword.load_list[i].curve_id;
		magnitude = keyword.curve_list[curve_id - 1].v2[0];      //Only for constant load
		new_load.magnitude = magnitude;
		_external_load_list.push_back(new_load);
	}
	//------------------------------------------------------------------------------
	//Read initial velocity
	//------------------------------------------------------------------------------
	if (keyword.is_initial_vel)
	{
		for (int i = 0; i < _nump; i++)
		{
			_nodeptr_list[i]->_velocity = keyword.initial_vel;
		}
	}
	//---------------------------------------------------------
	//Taylor bar test
	//---------------------------------------------------------
	for (int i = 0; i < _nump; i++)
	{
		if (abs(_nodeptr_list[i]->_position.z) < 1e-10)
		{
			_nodeptr_list[i]->_bc_type_position[2] = 1;
		}
	}
	return;

}