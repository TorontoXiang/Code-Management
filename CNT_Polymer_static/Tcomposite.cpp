#include "Tcomposite.h"
void Tcomposite::Input()
{
	_grid_Polymer.Input_Polymer();
	_grid_Polymer.Calculate_stiffness_matrix();
	cout << "Input CNT" << endl;
	input_CNT_list();
}
void Tcomposite::Calculate_Fp()
{
	_grid_Polymer.initialize_Fp();
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < _num_CNT; i++)
	{
		_CNT_list[i].calculate_CNT_boundary_p(&_grid_Polymer);
		_CNT_list[i].calculate_load_from_constraint();
		_CNT_list[i].Solving_equilibrium_equation();
		_CNT_list[i].calculate_reacting_force();
	}
	for (int i = 0; i < _num_CNT; i++)
	{
		_grid_Polymer.assemble_Fp(&_CNT_list[i]);
	}
}
void Tcomposite::CG_iteration()
{
	CG_iteration_initialization();
	double diff = _grid_Polymer.calculate_diff();
	if (diff < 1e-12)
	{
		return;
	}
	int num_iteration = 0;
	cout << "The 0th iteration step with error = " << diff << endl;
	while (diff > 1e-12)
	{
		Calculate_Fp();
		_grid_Polymer.Calculate_Kp();
		_grid_Polymer.CG_update();
		num_iteration = num_iteration + 1;
		diff = _grid_Polymer.calculate_diff();
		if (num_iteration % 1000 == 0)
		{
			cout << "The " << num_iteration << " iteration step with error = " << diff << endl;
		}
	}
	cout << num_iteration << " steps are used to reach error = " << diff << endl;
	return;
}
void Tcomposite::output_result()
{
	_grid_Polymer.Output_Tecplot("Polymer.dat");
	ofstream output;
	output.open("CNT.dat");
	for (int i = 0; i < _num_CNT; i++)
	{
		_CNT_list[i].calculate_CNT_boundary_displacement(&_grid_Polymer);
		_CNT_list[i].calculate_load_from_constraint();
		_CNT_list[i].Solving_equilibrium_equation();
		_CNT_list[i].output_tecplot(output);
	}
	//Output the external load
	ofstream output_force;
	output_force.open("Racting_force.k");
	Calculate_reacting_force();
	vec3D external_force1, external_force2;
	external_force1 = _grid_Polymer.calculate_external_load("x_min");
	external_force2 = _grid_Polymer.calculate_external_load("x_max");
	output_force << "External force at x_min is: " << external_force1.x << " " << external_force1.y << " " << external_force1.z << endl;
	output_force << "External force at x_max is: " << external_force2.x << " " << external_force2.y << " " << external_force2.z << endl;
	external_force1 = _grid_Polymer.calculate_external_load("y_min");
	external_force2 = _grid_Polymer.calculate_external_load("y_max");
	output_force << "External force at y_min is: " << external_force1.x << " " << external_force1.y << " " << external_force1.z << endl;
	output_force << "External force at y_max is: " << external_force2.x << " " << external_force2.y << " " << external_force2.z << endl;
	external_force1 = _grid_Polymer.calculate_external_load("z_min");
	external_force2 = _grid_Polymer.calculate_external_load("z_max");
	output_force << "External force at z_min is: " << external_force1.x << " " << external_force1.y << " " << external_force1.z << endl;
	output_force << "External force at z_max is: " << external_force2.x << " " << external_force2.y << " " << external_force2.z << endl;
}
void Tcomposite::CG_iteration_initialization()
{
	Calculate_Fp();
	_grid_Polymer.Set_F0();
	_grid_Polymer.calculate_load_from_constraint();
	_grid_Polymer.CG_initialization();
}
void Tcomposite::input_CNT_list()
{
	ifstream input;
	input.open("CNT_grid.k");
	string a;
	string temp;
	vector<Skeyword> keyword_list;
	while (true)
	{
		next_keyword(input);
		input >> a;
		if (a == "*NEW_CNT")
		{
			Skeyword new_keyword;
			while (true)
			{
				input >> a;
				if (a == "*NODE")
				{
					while (!exam_keyword(input))
					{
						Snode_temp new_node;
						input >> new_node.id >> new_node.x >> new_node.y >> new_node.z;
						new_keyword.node_list.push_back(new_node);
					}
				}
				else if (a == "*MAT_ELASTIC")
				{
					Smaterial new_material;
					double temp;
					next_data(input);
					input >> new_material.material_id;
					input >> temp;
					input >> new_material.Youngs >> new_material.Possion;
					//Add new_material into material_list and keep in order
					unsigned int material_id = new_material.material_id;
					if (material_id > new_keyword.material_list.size())
					{
						new_keyword.material_list.resize(material_id);
					}
					new_keyword.material_list[material_id - 1] = new_material;
				}
				else if (a == "*ELEMENT_SOLID")
				{
					while (!exam_keyword(input))
					{
						Scell_8 new_cell;
						input >> new_cell.cell_id >> new_cell.mid;
						input >> new_cell.IEN[0] >> new_cell.IEN[1] >> new_cell.IEN[2] >> new_cell.IEN[3];
						input >> new_cell.IEN[4] >> new_cell.IEN[5] >> new_cell.IEN[6] >> new_cell.IEN[7];
						new_keyword.cell_8_list.push_back(new_cell);
					}
				}
				else if (a=="*END_CNT")
				{
					break;
				}
			}
			keyword_list.push_back(new_keyword);
		}
		else if (a == "*END")
		{
			break;
		}
	}
	_num_CNT = keyword_list.size();
	_CNT_list.resize(_num_CNT);
	int num_CNT_element = 0,num_CNT_node = 0;
	for (int i = 0; i < _num_CNT; i++)
	{
		_CNT_list[i].Input_CNT(keyword_list[i]);
		_CNT_list[i].Create_MKL_solver();
		_CNT_list[i].calculate_CNT_location(&_grid_Polymer);
		num_CNT_element = num_CNT_element + _CNT_list[i].G_nume();
		num_CNT_node = num_CNT_node + _CNT_list[i].G_nump();
	}
	cout << "Number of CNT element is: " << num_CNT_element << endl;
	cout << "Number of CNT node is: " << num_CNT_node << endl;
	return;
}
void Tcomposite::Calculate_reacting_force()
{
	_grid_Polymer.calculate_reacting_force();
	for (int i = 0; i < _num_CNT; i++)
	{
		_CNT_list[i].calculate_CNT_boundary_displacement(&_grid_Polymer);
		_CNT_list[i].calculate_load_from_constraint();
		_CNT_list[i].Solving_equilibrium_equation();
		_CNT_list[i].calculate_reacting_force();
		_grid_Polymer.assemble_reacting_froce_from_CNT(&_CNT_list[i]);
	}
}