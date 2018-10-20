#include "Tbody.h"
#include "public_function.h"
#include <string>
#include <iomanip>
#include "TMAT_base.h"
#include "readin.h"
void Tbody::calculate_force_in_body()
{
	//Apply the external force
	apply_external_force();
	//Calcualte the corner force for every cell
	calculate_corner_force();
	//Calculate the final nodal force
	assemble_nodal_force();
	return;
}
void Tbody::assemble_nodal_force()
{
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->assemble_corner_force();
	}
	return;
}
void Tbody::update_grid()
{
	double dt_vel = 0.5*(_dt + _dt_pre);
	for (int i = 0; i < _nump; i++)
	{
		//Calculate the velocity at the next time step
		_nodeptr_list[i]->update_vel(dt_vel);
		//Calculate the position at the next time point
		_nodeptr_list[i]->update_pos(_dt);
	}
	return;
}
void Tbody::calculate_corner_force()
{
	for (int i = 0; i < _nume; i++)
	{
		if (i==200)
		{
			double s = 1;
		}
		_cellptr_list[i]->calculate_corner_force(_dt_pre);
	}
	return;
}
void Tbody::calculate_final_acceleration()
{
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->calculate_acceleration();
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
			_nodeptr_list[id - 1]->apply_external_force(direction, magnitude);
		}
	}
	//Apply gravity
	for (int i = 0; i < _nump; i++)
	{
		_nodeptr_list[i]->apply_gravity(_gravity);
	}
}
void Tbody::calculate_time_step()
{
	double min_time_step = 100000;
	for (int i = 0; i < _nume; i++)
	{
		min_time_step = minval(min_time_step, _cellptr_list[i]->calculate_time_step());
	}
	_dt_pre = _dt;
	_dt = _CFL * min_time_step;
	if (_num_step!=0)
	{
		_dt=minval(_CFL * min_time_step,1.05*_dt_pre);
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
	for (int i = 0; i < _nume; i++)
	{
		_cellptr_list[i]->assemble_corner_inertance();
	}
	return;
}
double Tbody::calculate_total_internal_energy()
{
	double total_internal_energy = 0;
	for (int i = 0; i < _nume; i++)
	{
		total_internal_energy = total_internal_energy + _cellptr_list[i]->G_strain_energy();
	}
	return total_internal_energy;
}
//------------------------------------------------------------------------
//Output functuions
//------------------------------------------------------------------------
void Tbody::output_tecplot(ofstream& output, double ratio)
{
	output << "TITLE = \"Mesh for Shell Element\"" << endl;
	output << "VARIABLES = \"X\",\"Y\",\"Z\"" << endl;
	output << "ZONE F=FEPOINT,N=" << _nump << "," << "E=" << _nume << "," << "ET=BRICK" << endl;
	vec3D dis;
	for (int i = 0; i < _nump; i++)
	{
		dis = _nodeptr_list[i]->calculate_displacement_amplify(ratio);
		output << dis.x << " " << dis.y << " " << dis.z << endl;
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
void Tbody::input_body(ifstream& input)
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
	_nume = keyword.cell_8_list.size();
	_nump = keyword.node_list.size();
	_endtime = keyword.time_control.endtime;
	_CFL = keyword.time_control.CFL;
	_num_step = 0;
	//----------------------------------------------------
	//Read node list
	//----------------------------------------------------
	int id;
	double x, y, z;
	_grid_polymer._node_list.resize(_nump);
	_nodeptr_list.resize(_nump);
	for (int i = 0; i < _nump; i++)
	{
		id = keyword.node_list[i].id;
		x = keyword.node_list[i].x;
		y = keyword.node_list[i].y;
		z = keyword.node_list[i].z;
		Tnode new_node(id, x, y, z);
		_grid_polymer._node_list[i] = new_node;
		_nodeptr_list[i] = &_grid_polymer._node_list[i];
	}
	//-----------------------------------------------------
	//Read cell list
	//-----------------------------------------------------
	_grid_polymer._cell_list.resize(_nume);
	_cellptr_list.resize(_nume);
	int part_id, material_id;
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
		Smaterial this_mat = keyword.material_list[material_id - 1];
		TMAT_base* mat;
		mat = generate_material(this_mat, 0);
		//Create a new cell
		Tcell_brick new_cell(id, node_ptr, mat);
		_grid_polymer._cell_list[i] = new_cell;
		_cellptr_list[i] = &_grid_polymer._cell_list[i];
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
		for (int i = 0; i < _nump; i++)
		{
			_nodeptr_list[i]->_velocity = keyword.initial_vel;
		}
	}
	for (int i = 0; i < _nump; i++)
	{
		if (abs(_nodeptr_list[i]->_position.z) < 1e-10)
		{
			_nodeptr_list[i]->_bc_type_position[2] = 1;
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
	//Read the information about the regular grid
	_x_min = keyword.regular_grid.x_min; _x_max = keyword.regular_grid.x_max;
	_nx = keyword.regular_grid.nx; _ny = keyword.regular_grid.ny; _nz = keyword.regular_grid.nz;
	_nump = (_nx + 1)*(_ny + 1)*(_nz + 1);
	_nume = _nx * _ny*_nz;
	_grid_polymer._node_list.resize(_nump);
	_nodeptr_list.resize(_nump);
	_grid_polymer._cell_list.resize(_nume);
	_cellptr_list.resize(_nume);
	_interval.x = (_x_max - _x_min).x / _nx;
	_interval.y = (_x_max - _x_min).x / _ny;
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
				Tnode new_node(index, coor.x, coor.y, coor.z);
				_grid_polymer._node_list[index] = new_node;
				_nodeptr_list[index] = &_grid_polymer._node_list[index];
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
					node_ptr[n] = _nodeptr_list[index_node[n]];
				}
				Smaterial this_mat = keyword.material_list[0];
				TMAT_base* mat;
				mat = generate_material(this_mat, 0);
				Tcell_brick new_cell(index, node_ptr, mat);
				_grid_polymer._cell_list[index] = new_cell;
				_cellptr_list[index] = &_grid_polymer._cell_list[index];
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
				if (node_id >= _nump)
				{
					cout << "Error: Node ID" << node_id << " excess" << endl;
					cout << "Number of node in body is " << _nump << endl;
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
				if (cell_id >= _nume)
				{
					cout << "Error: Cell ID " << cell_id << " excess" << endl;
					cout << "Number of cell in body is " << _nume << endl;
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
				Tnode* nodeptr = _nodeptr_list[index];
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