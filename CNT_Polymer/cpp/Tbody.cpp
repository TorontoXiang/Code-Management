#include "Tbody.h"
#include "public_function.h"
#include <string>
#include <iomanip>
#include "TMAT_base.h"
#include "readin.h"
#include "Shape_function.h"
void Tbody::Snode_stress::calculate_average()
{
	sxx = sxx / node_volume; syy = syy / node_volume; szz = szz / node_volume;
	sxy = sxy / node_volume; sxz = sxz / node_volume; syz = syz / node_volume;
	eqs = eqs / node_volume;
	return;
}
void Tbody::calculate_force_in_body()
{
	//Apply the external force
	apply_external_force();
	//Calcualte the corner force for every cell
	calculate_corner_force();
	//Calculate the final nodal force
	assemble_nodal_force();
	//Calculate the froce in the CNT
	calculate_CNT_force();
	//Apply the froce from the CNT
	for (int i = 0; i < _nump_CNT; i++)
	{
		contribute_CNT_force(_nodeptr_list_CNT[i]);
	}
	return;
}
void Tbody::assemble_nodal_force()
{
	for (int i = 0; i < _nume_Polymer; i++)
	{
		_cellptr_list_Polymer[i]->assemble_corner_force();
	}
	double ratio = 1;
	for (int i = 0; i < _nump_Polymer; i++)
	{
		_nodeptr_list_Polymer[i]->_force = _nodeptr_list_Polymer[i]->_force - _nodeptr_list_Polymer[i]->_velocity*ratio;
	}
	return;
}
void Tbody::update_grid()
{
	double dt_vel = 0.5*(_dt + _dt_pre);
	for (int i = 0; i < _nump_Polymer; i++)
	{
		//Calculate the velocity at the next time step
		_nodeptr_list_Polymer[i]->update_vel(dt_vel);
		//Calculate the position at the next time point
		_nodeptr_list_Polymer[i]->update_pos(_dt);
	}
	//Update the position of CNT
	for (int i = 0; i < _nump_CNT; i++)
	{
		update_CNT_position(_nodeptr_list_CNT[i]);
	}
	return;
}
void Tbody::calculate_corner_force()
{
	for (int i = 0; i < _nume_Polymer; i++)
	{
		_cellptr_list_Polymer[i]->calculate_corner_force(_dt_pre);
	}
	return;
}
void Tbody::calculate_final_acceleration()
{
	for (int i = 0; i < _nump_Polymer; i++)
	{
		_nodeptr_list_Polymer[i]->calculate_acceleration();
	}
	return;
}
void Tbody::apply_external_force()
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
			_nodeptr_list_Polymer[id - 1]->apply_external_force(direction, magnitude);
		}
	}
	//Apply gravity
	for (int i = 0; i < _nump_Polymer; i++)
	{
		_nodeptr_list_Polymer[i]->apply_gravity(_gravity);
	}
}
void Tbody::calculate_time_step()
{
	double min_time_step = 100000;
	for (int i = 0; i < _nume_Polymer; i++)
	{
		min_time_step = minval(min_time_step, _cellptr_list_Polymer[i]->calculate_time_step());
	}
	_dt_pre = _dt;
	_dt = _CFL * min_time_step;
	if (_num_step!=0)
	{
		_dt=minval(_CFL * min_time_step,1.05*_dt_pre);
	}
	if (!_is_static)
	{
		_dt = minval(_dt, _motiontime - _current_time);
	}
	_dt = minval(_dt, _endtime - _current_time);
	return;
}
void Tbody::start_up_step()
{
	calculate_force_in_body();
	return;
}
void Tbody::calculate_nodal_inertance()
{
	for (int i = 0; i < _nume_Polymer; i++)
	{
		_cellptr_list_Polymer[i]->assemble_corner_inertance();
	}
	for (int i = 0; i < _nume_CNT; i++)
	{
		contribute_CNT_mass(_cellptr_list_CNT[i]);
	}
	return;
}
double Tbody::calculate_total_internal_energy()
{
	double total_internal_energy = 0;
	for (int i = 0; i < _nume_Polymer; i++)
	{
		total_internal_energy = total_internal_energy + _cellptr_list_Polymer[i]->G_strain_energy();
	}
	return total_internal_energy;
}
//------------------------------------------------------------------------
//Output functuions
//------------------------------------------------------------------------
void Tbody::output_tecplot(ofstream& output, double ratio)
{
	output << "TITLE = \"Mesh for Shell Element\"" << endl;
	output << "VARIABLES = \"X\",\"Y\",\"Z\",\"DX\",\"DY\",\"DZ\",\"VX\",\"VY\",\"VZ\",\"SXX\",\"SYY\",\"SZZ\",\"SXY\",\"SXZ\",\"SYZ\",\"EQS\"" << endl;
	output << "ZONE F=FEPOINT,N=" << _nump_Polymer << "," << "E=" << _nume_Polymer << "," << "ET=BRICK" << endl;
	vec3D pos, dis, vel;
	Snode_stress* node_stress;
	node_stress = new Snode_stress[_nump_Polymer];
	double eqs;
	//Calculate the nodal stress
	for (int i = 0; i < _nume_Polymer; i++)
	{
		Tcell_brick* cell_ptr=_cellptr_list_Polymer[i];
		double sig[6], volume = cell_ptr->_volume*0.125;
		cell_ptr->_mat_ptr->G_stress(sig);
		double sxx = sig[0], syy = sig[1], szz = sig[2], sxy = sig[3], sxz = sig[4], syz = sig[5];
		double eqs = cell_ptr->_mat_ptr->calculate_equivalent_stress();
		for (int j = 0; j < 8; j++)
		{
			int node_id = cell_ptr->_node_ptr[j]->_id;
			node_stress[node_id - 1].node_volume = node_stress[node_id - 1].node_volume + volume;
			node_stress[node_id - 1].sxx = node_stress[node_id - 1].sxx + sxx*volume;
			node_stress[node_id - 1].syy = node_stress[node_id - 1].syy + syy * volume;
			node_stress[node_id - 1].szz = node_stress[node_id - 1].szz + szz * volume;
			node_stress[node_id - 1].sxy = node_stress[node_id - 1].sxy + sxy * volume;
			node_stress[node_id - 1].sxz = node_stress[node_id - 1].sxz + sxz * volume;
			node_stress[node_id - 1].syz = node_stress[node_id - 1].syz + syz * volume;
			node_stress[node_id - 1].eqs = node_stress[node_id - 1].eqs + eqs * volume;
		}
	}
	for (int i = 0; i < _nump_Polymer; i++)
	{
		node_stress[i].calculate_average();
	}
	for (int i = 0; i < _nump_Polymer; i++)
	{
		pos = _nodeptr_list_Polymer[i]->calculate_displacement_amplify(ratio);
		vel = _nodeptr_list_Polymer[i]->_velocity;
		dis = _nodeptr_list_Polymer[i]->_position - _nodeptr_list_Polymer[i]->_position0;
		output << pos.x << " " << pos.y << " " << pos.z << " " << dis.x << " " << dis.y << " " << dis.z<<" ";
		output << vel.x << " " << vel.y << " " << vel.z << " ";
		output << node_stress[i].sxx << " " << node_stress[i].syy << " " << node_stress[i].szz<<" ";
		output << node_stress[i].sxy << " " << node_stress[i].sxz << " " << node_stress[i].syz<<" ";
		output << node_stress[i].eqs << endl;
	}
	for (int i = 0; i < _nume_Polymer; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			output << _cellptr_list_Polymer[i]->_node_ptr[j]->_id << " ";
		}
		output << endl;
	}
	return;
}
//------------------------------------------------------------------------------------
//Input function
//------------------------------------------------------------------------------------
void Tbody::input_Polymer(ifstream& input)
{
	Skeyword keyword;
	read_in_keyword_file(input, keyword);
	if (!keyword.is_regular_grid)
	{
		input_irregular_grid(keyword);
	}
	else
	{
		input_regular_grid(keyword);
	}
	return;
}
void Tbody::input_irregular_grid(Skeyword& keyword)
{
	//Read the body global information from keyword
	_nume_Polymer = keyword.cell_8_list.size();
	_nump_Polymer = keyword.node_list.size();
	_endtime = keyword.time_control.endtime;
	_CFL = keyword.time_control.CFL;
	_num_step = 0;
	_is_static = false;
	_motiontime = keyword.time_control.motiontime;
	//----------------------------------------------------
	//Read node list
	//----------------------------------------------------
	int id;
	double x, y, z;
	_grid_polymer._node_list.resize(_nump_Polymer);
	_nodeptr_list_Polymer.resize(_nump_Polymer);
	for (int i = 0; i < _nump_Polymer; i++)
	{
		id = keyword.node_list[i].id;
		x = keyword.node_list[i].x;
		y = keyword.node_list[i].y;
		z = keyword.node_list[i].z;
		Tnode new_node(id, x, y, z);
		_grid_polymer._node_list[i] = new_node;
		_nodeptr_list_Polymer[i] = &_grid_polymer._node_list[i];
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	_grid_polymer._cell_list.resize(_nume_Polymer);
	_cellptr_list_Polymer.resize(_nume_Polymer);
	int part_id, material_id;
	Tnode* node_ptr[8];
	for (int i = 0; i < _nume_Polymer; i++)
	{
		//Read cell id and IEN
		id = keyword.cell_8_list[i].cell_id;
		for (int j = 0; j < 8; j++)
		{
			int IENj = keyword.cell_8_list[i].IEN[j];
			node_ptr[j] = _nodeptr_list_Polymer[IENj - 1];
		}
		//Read cell's material properties
		part_id = keyword.cell_8_list[i].part_id;
		material_id = keyword.part_list[part_id - 1].material_id;
		Smaterial this_mat = keyword.material_list[material_id - 1];
		TMAT_base* mat;
		mat = generate_material(this_mat, 0);
		//Create a new cell
		Tcell_brick new_cell(id, node_ptr, mat);
		_grid_polymer._cell_list[i] = new_cell;
		_cellptr_list_Polymer[i] = &_grid_polymer._cell_list[i];
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
					_grid_polymer._node_list[keyword.node_group_list[node_group_id - 1].node_id[k] - 1]._bc_type_position[j] = 1;
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
		for (int i = 0; i < _nump_Polymer; i++)
		{
			_nodeptr_list_Polymer[i]->_velocity = keyword.initial_vel;
		}
	}
	for (int i = 0; i < _nump_Polymer; i++)
	{
		if (abs(_nodeptr_list_Polymer[i]->_position.z) < 1e-10)
		{
			_nodeptr_list_Polymer[i]->_bc_type_position[2] = 1;
		}
	}
	return;
}
void Tbody::input_regular_grid(Skeyword& keyword)
{
	//Read the body global information from keyword
	_endtime = keyword.time_control.endtime;
	_CFL = keyword.time_control.CFL;
	_num_step = 0;
	_is_static = false;
	_motiontime = keyword.time_control.motiontime;
	//Read the information about the regular grid
	_x_min = keyword.regular_grid.x_min; _x_max = keyword.regular_grid.x_max;
	_nx = keyword.regular_grid.nx; _ny = keyword.regular_grid.ny; _nz = keyword.regular_grid.nz;
	_nump_Polymer = (_nx + 1)*(_ny + 1)*(_nz + 1);
	_nume_Polymer = _nx * _ny*_nz;
	_grid_polymer._node_list.resize(_nump_Polymer);
	_nodeptr_list_Polymer.resize(_nump_Polymer);
	_grid_polymer._cell_list.resize(_nume_Polymer);
	_cellptr_list_Polymer.resize(_nume_Polymer);
	_interval.x = (_x_max - _x_min).x / _nx;
	_interval.y = (_x_max - _x_min).y / _ny;
	_interval.z = (_x_max - _x_min).z / _nz;
	//Generate the node coordinate
	vec3D dx,coor;
	int index;
	for (int i = 0; i < _nx + 1; i++)
	{
		for (int j = 0; j < _ny + 1; j++)
		{
			for (int k = 0; k < _nz + 1; k++)
			{
				index = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*j + k;
				dx.value(i*_interval.x, j*_interval.y, k*_interval.z);
				coor = _x_min + dx;
				Tnode new_node(index+1, coor.x, coor.y, coor.z);
				_grid_polymer._node_list[index] = new_node;
				_nodeptr_list_Polymer[index] = &_grid_polymer._node_list[index];
			}
		}
	}
	//Generate the cell
	int index_node[8];
	Tnode* node_ptr[8];
	for (int i = 0; i < _nx; i++)
	{
		for (int j = 0; j < _ny; j++)
		{
			for (int k = 0; k < _nz; k++)
			{
				index = _ny * _nz*i + _nz * j + k; //The index convertion for cell
				index_node[0] = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*j + k; index_node[1] = (_ny + 1)*(_nz + 1)*(i + 1) + (_nz + 1)*j + k;
				index_node[2] = (_ny + 1)*(_nz + 1)*(i + 1) + (_nz + 1)*(j + 1) + k; index_node[3] = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*(j + 1) + k;
				index_node[4] = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*j + k + 1; index_node[5] = (_ny + 1)*(_nz + 1)*(i + 1) + (_nz + 1)*j + k + 1;
				index_node[6] = (_ny + 1)*(_nz + 1)*(i + 1) + (_nz + 1)*(j + 1) + k + 1; index_node[7] = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*(j + 1) + k + 1;
				for (int n = 0; n < 8; n++)
				{
					node_ptr[n] = _nodeptr_list_Polymer[index_node[n]];
				}
				Smaterial this_mat = keyword.material_list[0];
				TMAT_base* mat;
				mat = generate_material(this_mat, 0);
				Tcell_brick new_cell(index, node_ptr, mat);
				_grid_polymer._cell_list[index] = new_cell;
				_cellptr_list_Polymer[index] = &_grid_polymer._cell_list[index];
			}
		}
	}
	int num_bc = keyword.regular_bc_list.size();
	for (int i = 0; i < num_bc; i++)
	{
		Sregular_grid_bc bc = keyword.regular_bc_list[i];
		apply_regular_bc(bc);
	}
	return;
}
void Tbody::input_CNT(ifstream& input)
{
	Skeyword keyword;
	read_in_keyword_file(input, keyword);
	//Read the body global information from keyword
	_nume_CNT = keyword.cell_2_list.size();
	_nump_CNT = keyword.node_list.size();
	//----------------------------------------------------
	//Read node list
	//----------------------------------------------------
	int id;
	double x, y, z;
	_grid_CNT._node_list.resize(_nump_CNT);
	_nodeptr_list_CNT.resize(_nump_CNT);
	for (int i = 0; i < _nump_CNT; i++)
	{
		id = keyword.node_list[i].id;
		x = keyword.node_list[i].x;
		y = keyword.node_list[i].y;
		z = keyword.node_list[i].z;
		Tnode_CNT new_node(id, x, y, z);
		_grid_CNT._node_list[i] = new_node;
		_nodeptr_list_CNT[i] = &_grid_CNT._node_list[i];
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	_grid_CNT._cell_list.resize(_nume_CNT);
	_cellptr_list_CNT.resize(_nume_CNT);
	int part_id, material_id, section_id;
	Tnode_CNT* node_ptr[2];
	for (int i = 0; i < _nume_CNT; i++)
	{
		//Read cell id and IEN
		id = keyword.cell_2_list[i].cell_id;
		for (int j = 0; j < 2; j++)
		{
			int IENj = keyword.cell_2_list[i].IEN[j];
			node_ptr[j] = _nodeptr_list_CNT[IENj - 1];
		}
		//Read cell's material properties
		part_id = keyword.cell_2_list[i].part_id;
		material_id = keyword.part_list[part_id - 1].material_id;
		section_id = keyword.part_list[part_id - 1].section_id;
		double area = keyword.section_list[section_id - 1].area;
		Smaterial this_mat = keyword.material_list[material_id - 1];
		TMAT_base* mat;
		mat = generate_material(this_mat, 0);
		//Create a new cell
		Tcell_CNT new_cell(id, node_ptr, area, mat);
		_grid_CNT._cell_list[i] = new_cell;
		_cellptr_list_CNT[i] = &_grid_CNT._cell_list[i];
	}
	calculate_CNT_location();
	return;
}
void Tbody::create_curve(int type, int id)
{
	Scurve_output* new_curve;
	new_curve = new Scurve_output(type, id);
	_curve_out_list.push_back(new_curve);
	return;
}
void Tbody::create_Tecplot_mesh(double amplify_ratio, int total_frame)
{
	Smesh_output* new_mesh = NULL;
	new_mesh = new Smesh_output(amplify_ratio, total_frame, _endtime);
	_mesh_out = new_mesh;
	return;
}
void Tbody::input_simulation_control()
{
	ifstream input;
	string a;
	input.open("SimulationControl.dat");
	if (input.fail())
	{
		cout << "Error: No extra information input" << endl;
		system("Pause");
		exit(0);
	}
	while (true)
	{
		next_keyword(input);
		input >> a;
		if (a == "*CURVE_NODE")
		{
			int node_id;
			while (!exam_keyword(input))
			{
				input >> node_id;
				if (node_id >= _nump_Polymer)
				{
					cout << "Error: Node ID" << node_id << " excess" << endl;
					cout << "Number of node in body is " << _nump_Polymer << endl;
					system("Pause");
					exit(0);
				}
				create_curve(0, node_id);
			}
		}
		else if (a == "*CURVE_CELL")
		{
			int cell_id;
			while (!exam_keyword(input))
			{
				input >> cell_id;
				if (cell_id >= _nume_Polymer)
				{
					cout << "Error: Cell ID " << cell_id << " excess" << endl;
					cout << "Number of cell in body is " << _nume_Polymer << endl;
					system("Pause");
					exit(0);
				}
				create_curve(1, cell_id);
			}
		}
		else if (a == "*TECPLOT_MESH_OUTPUT")
		{
			int total_frame;
			double amplify_ratio;
			while (!exam_keyword(input))
			{
				input >> amplify_ratio >> total_frame;
				create_Tecplot_mesh(amplify_ratio, total_frame);
			}
		}
		else if (a=="*END")
		{
			input.close();
			return;
		}
	}
	return;
}
void Tbody::Output_Tecplot_mesh(string type)
{
	if (_mesh_out == NULL)
	{
		return;
	}
	if (type == "Initial")
	{
		cout << "Output grid at initial"<< endl;
		output_tecplot(_mesh_out->output_initial, _mesh_out->_amplify_ratio);
	}
	else if (type == "Process")
	{

		if (_mesh_out->is_output(_current_time, _dt_pre))
		{
			cout << "Output grid at time=" << setw(8) << _current_time<< endl;
			output_tecplot(_mesh_out->output_process, _mesh_out->_amplify_ratio);
		}
	}
	else if (type == "Final")
	{
		cout << "Output the final grid"<< endl;
		output_tecplot(_mesh_out->output_final, _mesh_out->_amplify_ratio);
	}
	else
	{
		cout << "Error: Invalid type in output mesh" << endl;
	}
	return;
}
void Tbody::Output_Curve()
{
	int num_curve = _curve_out_list.size();
	int type, id;
	Tnode* node_ptr;
	//Tcell_brick* cell_ptr;
	for (int i = 0; i < num_curve; i++)
	{
		type = _curve_out_list[i]->_type;
		id = _curve_out_list[i]->_id - 1;
		_curve_out_list[i]->output.precision(10);
		if (type == 0)
		{
			node_ptr = _nodeptr_list_Polymer[id];
			vec3D dis = node_ptr->calculate_displacement();
			_curve_out_list[i]->output << setw(20) << _current_time << " ";
			_curve_out_list[i]->output << setw(20) << dis.x << setw(20) << dis.y << setw(20) << dis.z << endl;
		}
	}
	return;
}
void Tbody::apply_regular_bc(Sregular_grid_bc& bc)
{
	int nx_min = 0, ny_min = 0, nz_min = 0;
	int nx_max = _nx + 1, ny_max = _ny + 1, nz_max = _nz + 1;
	if (bc.FaceName=="x_min")
	{
		nx_min = 0; nx_max = 1;
	}
	else if (bc.FaceName == "x_max")
	{
		nx_min = _nx; nx_max = _nx + 1;
	}
	else if (bc.FaceName == "y_min")
	{
		ny_min = 0; ny_max = 1;
	}
	else if (bc.FaceName == "y_max")
	{
		ny_min = _ny; ny_max = _ny + 1;
	}
	else if (bc.FaceName == "z_min")
	{
		nz_min = 0; nz_max = 1;
	}
	else if (bc.FaceName == "z_max")
	{
		nz_min = _nz; nz_max = _nz + 1;
	}
	else
	{
		cout << "Invalid boundary contition for regular grid" << endl;
		system("Pause");
		exit(0);
	}
	for (int i = nx_min; i < nx_max; i++)
	{
		for (int j = ny_min; j < ny_max; j++)
		{
			for (int k = nz_min; k < nz_max; k++)
			{
				int index = (_ny + 1)*(_nz + 1)*i + (_nz + 1)*j + k;
				Tnode* nodeptr = _nodeptr_list_Polymer[index];
				if (bc.bc_type==0)
				{
					nodeptr->_bc_type_position[bc.direction] = 1;
				}
				else if (bc.bc_type==1)
				{
					nodeptr->_bc_type_position[bc.direction] = 2;
					nodeptr->_bc_velocity[bc.direction] = bc.velocity;
				}
			}
		}
	}
	return;
}
void Tbody::change_state()
{
	if (abs(_current_time-_motiontime)<1e-15)
	{
		for (int i = 0; i < _nump_Polymer; i++)
		{
			Tnode* node_ptr = _nodeptr_list_Polymer[i];
			for (int j = 0; j < 3; j++)
			{
				if (node_ptr->_bc_type_position[j]==2)
				{
					node_ptr->_bc_velocity[j] = 0;
				}
			}
		}
		_is_static = true;
	}
	return;
}
void Tbody::calculate_average_stress(double(&stress)[6])
{
	double sxx = 0, syy = 0, szz = 0, sxz = 0, sxy = 0, syz = 0;
	double total_volume = 0, cell_volume;
	Tcell_brick* cell_ptr = NULL;
	double sigma[6];
	for (int i = 0; i < _nume_Polymer; i++)
	{
		cell_ptr = _cellptr_list_Polymer[i];
		cell_volume = cell_ptr->_volume;
		total_volume = total_volume + cell_volume;
		cell_ptr->_mat_ptr->G_stress(sigma);
		sxx = sxx + sigma[0] * cell_volume; syy = syy + sigma[1] * cell_volume;
		szz = szz + sigma[2] * cell_volume; sxy = sxy + sigma[3] * cell_volume;
		sxz = sxz + sigma[4] * cell_volume; syz = syz + sigma[5] * cell_volume;
	}
	sxx = sxx / total_volume; syy = syy / total_volume; szz = szz / total_volume;
	sxy = sxy / total_volume; sxz = sxz / total_volume; syz = syz / total_volume;
	stress[0] = sxx; stress[1] = syy; stress[2] = szz;
	stress[3] = sxy; stress[4] = sxz; stress[5] = syz;
	cout << total_volume << endl;
	return;
}
void Tbody::calculate_CNT_location()
{
	vec3D coor;
	for (int i = 0; i < _nump_CNT; i++)
	{
		CNT_location(_nodeptr_list_CNT[i]);
	}
	for (int i = 0; i < _nume_CNT; i++)
	{
		CNT_location(_cellptr_list_CNT[i]);
	}
	return;
}
void Tbody::CNT_location(Tnode_CNT* node_ptr)
{
	vec3D coor = node_ptr->_position;
	vec3D delta = coor - _x_min;
	double ix, iy, iz;
	int nx = (int)floor(delta.x / _interval.x);
	int ny = (int)floor(delta.y / _interval.y);
	int nz = (int)floor(delta.z / _interval.z);
	nx = minval(nx, _nx - 1); ny = minval(ny, _ny - 1); nz = minval(nz, _nz - 1);
	node_ptr->_location_id = _ny * _nz*nx + _nz * ny + nz;
	double dx = coor.x - nx * _interval.x;
	ix = node_ptr->_ix = -1 + 2 * dx / _interval.x;
	double dy = coor.y - ny * _interval.y;
	iy = node_ptr->_iy = -1 + 2 * dy / _interval.y;
	double dz = coor.z - nz * _interval.z;
	iz = node_ptr->_iz = -1 + 2 * dz / _interval.z;
	//if (abs(ix)>1 || abs(iy)>1 || abs(iz)>1)
	//{
	//	cout << "Error: The standard coordiante excess!" << endl;
	//	system("Pause");
	//	exit(0);
	//}
	return;
}
void Tbody::CNT_location(Tcell_CNT* cell_ptr)
{
	vec3D coor = ((cell_ptr->_node_ptr[0]->_position) + (cell_ptr->_node_ptr[1]->_position))*0.5;
	vec3D delta = coor - _x_min;
	int nx = (int)floor(delta.x / _interval.x);
	int ny = (int)floor(delta.y / _interval.y);
	int nz = (int)floor(delta.z / _interval.z);
	nx = minval(nx, _nx - 1); ny = minval(ny, _ny - 1); nz = minval(nz, _nz - 1);
	cell_ptr->_location_id = _ny * _nz*nx + _nz * ny + nz;
	double dx = coor.x - nx * _interval.x;
	double ix, iy, iz;
	ix = cell_ptr->_ix = -1 + 2 * dx / _interval.x;
	double dy = coor.y - ny * _interval.y;
	iy = cell_ptr->_iy = -1 + 2 * dy / _interval.y;
	double dz = coor.z - nz * _interval.z;
	iz = cell_ptr->_iz = -1 + 2 * dz / _interval.z;
	//if (abs(ix) > 1 || abs(iy) > 1 || abs(iz) > 1)
	//{
	//	cout << "Error: The standard coordiante excess!" << endl;
	//	system("Pause");
	//	exit(0);
	//}
	return;
}
void Tbody::contribute_CNT_mass(Tcell_CNT* cell_ptr)
{
	double N[8], ix = cell_ptr->_ix, iy = cell_ptr->_iy, iz = cell_ptr->_iz;
	int Polymer_cell_id = cell_ptr->_location_id;
	double mass_CNT = cell_ptr->_l0*cell_ptr->_area*cell_ptr->_mat_ptr->G_density();
	Tcell_brick* cellptr_polymer = _cellptr_list_Polymer[Polymer_cell_id];
	for (int i = 0; i < 8; i++)
	{
		N[i] = N_brick_8(i, ix, iy, iz);
		cellptr_polymer->_node_ptr[i]->_mass = cellptr_polymer->_node_ptr[i]->_mass + N[i] * mass_CNT;
	}
	return;
}
void Tbody::contribute_CNT_force(Tnode_CNT* node_ptr)
{
	double N[8], ix = node_ptr->_ix, iy = node_ptr->_iy, iz = node_ptr->_iz;
	int Polymer_cell_id = node_ptr->_location_id;
	vec3D force_CNT2Polymer = node_ptr->_force;
	Tcell_brick* cellptr_polymer = _cellptr_list_Polymer[Polymer_cell_id];
	for (int i = 0; i < 8; i++)
	{
		N[i] = N_brick_8(i, ix, iy, iz);
		cellptr_polymer->_node_ptr[i]->_force = cellptr_polymer->_node_ptr[i]->_force + force_CNT2Polymer*N[i];
	}
	//Initialize the CNT force after the contribution
	node_ptr->_force.value(0, 0, 0);
	return;
}
void Tbody::update_CNT_position(Tnode_CNT* node_ptr)
{
	double N[8], ix = node_ptr->_ix, iy = node_ptr->_iy, iz = node_ptr->_iz;
	int Polymer_cell_id = node_ptr->_location_id;
	Tcell_brick* cellptr_polymer = _cellptr_list_Polymer[Polymer_cell_id];
	vec3D temp;
	for (int i = 0; i < 8; i++)
	{
		N[i] = N_brick_8(i, ix, iy, iz);
		temp = temp + cellptr_polymer->_node_ptr[i]->_position * N[i];
	}
	node_ptr->_position = temp;
	return;
}
void Tbody::calculate_CNT_force()
{
	for (int i = 0; i < _nume_CNT; i++)
	{
		_cellptr_list_CNT[i]->calculate_corner_force();
		_cellptr_list_CNT[i]->assemble_corner_force();
	}
	return;
}
void Tbody::output_CNT_stress()
{
	ofstream output;
	output.open("CNT_Result.dat");
	Tcell_CNT* cellptr;
	output << "Cell force result" << endl;
	for (int i = 0; i < _nume_CNT; i++)
	{
		cellptr = _cellptr_list_CNT[i];
		cellptr->calculate_corner_force();
		output << i<<" "<<cellptr->_corner_force[0].x << " "<<cellptr->_corner_force[0].y << " " << cellptr->_corner_force[0].z << endl;
		_cellptr_list_CNT[i]->assemble_corner_force();
	}
	output << "Node force result" << endl;
	Tnode_CNT* nodeptr;
	for (int i = 0; i < _nump_CNT; i++)
	{
		nodeptr = _nodeptr_list_CNT[i];
		output << i << " " << nodeptr->_force.x << " " << nodeptr->_force.y << " " << nodeptr->_force.z << endl;
	}
	output << "Node coordinate" << endl;
	for (int i = 0; i < _nump_CNT; i++)
	{
		nodeptr = _nodeptr_list_CNT[i];
		output << "k,"<< i << "," << nodeptr->_position.x << "," << nodeptr->_position.y << "," << nodeptr->_position.z << endl;
	}
}	
void Tbody::output_Polymer_stress()
{
	Tcell_brick* cellptr;
	ofstream output;
	output.open("Polymer_Result.dat");
	double sig[6];
	for (int i = 0; i < _nume_Polymer; i++)
	{
		cellptr = _cellptr_list_Polymer[i];
		cellptr->_mat_ptr->G_stress(sig);
		output << i << " " << sig[0] << " " << sig[1] << " " << sig[2] << " " << sig[3] << " " << sig[4] << " " << sig[5] << endl;
	}
	return;
}