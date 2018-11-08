#include "Tcomposite.h"
void Tcomposite::Input()
{
	Skeyword keyword;
	ifstream input;
	input.open("one_CNT.k");
	read_in_keyword_file(input, keyword);
	_grid_CNT.Input_CNT(keyword);
	_grid_CNT.Create_MKL_solver();

	_grid_Polymer.Input_Polymer();
	_grid_Polymer.Calculate_stiffness_matrix();
	_grid_CNT.calculate_CNT_location(&_grid_Polymer);
}
void Tcomposite::Calculate_Fp()
{
	_grid_Polymer.initialize_Fp();
	for (int i = 0; i < _num_CNT; i++)
	{
		_CNT_list[i].calculate_CNT_boundary_p(&_grid_Polymer);
		_CNT_list[i].calculate_load_from_constraint();
		_CNT_list[i].Solving_equilibrium_equation();
		_CNT_list[i].calculate_reacting_force();
		_grid_Polymer.assemble_Fp(&_CNT_list[i]);
	}

}
void Tcomposite::CG_iteration()
{
	CG_iteration_initialization();
	double diff = _grid_Polymer.calculate_diff();
	if (diff<1e-10)
	{
		return;
	}
	int num_iteration = 0;
	while (diff>1e-15)
	{
		Calculate_Fp();
		_grid_Polymer.Calculate_Kp();
		_grid_Polymer.CG_update();
		num_iteration = num_iteration + 1;
		diff = _grid_Polymer.calculate_diff();
		cout << num_iteration<<" "<<diff << endl;
	}
	return;
}
void Tcomposite::output_result()
{
	_grid_Polymer.Output_Tecplot("Polymer.dat");
	_grid_CNT.calculate_CNT_boundary_displacement(&_grid_Polymer);
	_grid_CNT.calculate_load_from_constraint();
	_grid_CNT.Solving_equilibrium_equation();
	_grid_CNT.Output_Tecplot("CNT.dat");
}
void Tcomposite::CG_iteration_initialization()
{
	Calculate_Fp();
	_grid_Polymer.Set_F0();
	_grid_Polymer.calculate_load_from_constraint();
	_grid_Polymer.CG_initialization();
}