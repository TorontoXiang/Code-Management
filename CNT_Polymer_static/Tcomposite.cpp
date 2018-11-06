#include "Tcomposite.h"
void Tcomposite::Input()
{
	Skeyword keyword;
	ifstream input;
	input.open("single_cell_move.k");
	read_in_keyword_file(input, keyword);
	_grid_CNT.Input_CNT(keyword);
	_grid_CNT.Create_MKL_solver();
	_grid_Polymer.Input_Polymer();
	_grid_Polymer.Create_MKL_solver();
	_grid_Polymer.calculate_load_from_constraint();
	_grid_CNT.calculate_CNT_location(&_grid_Polymer);
}
void Tcomposite::Calculate_F0_b()
{
	//Using zero Polymer displacement to calculate the constant term of _Fp
	_grid_CNT.calculate_CNT_boundary_displacement(&_grid_Polymer);
	_grid_CNT.calculate_load_from_constraint();
	_grid_CNT.Solving_equilibrium_equation();
	_grid_CNT.calculate_reacting_force();
	_grid_Polymer.initialize_load_total(1);
	_grid_Polymer.calculate_load_from_CNT(&_grid_CNT,1);
	_grid_Polymer.Calculate_b();
}
void Tcomposite::Calculate_Fp()
{
	_grid_CNT.calculate_CNT_boundary_p(&_grid_Polymer);
	_grid_CNT.calculate_load_from_constraint();
	_grid_CNT.Solving_equilibrium_equation();
	_grid_CNT.calculate_reacting_force();
	_grid_Polymer.initialize_load_total(2);
	_grid_Polymer.calculate_load_from_CNT(&_grid_CNT, 2);
}
void Tcomposite::Calculate_Fd()
{
	_grid_CNT.calculate_CNT_boundary_displacement(&_grid_Polymer);
	_grid_CNT.calculate_load_from_constraint();
	_grid_CNT.Solving_equilibrium_equation();
	_grid_CNT.calculate_reacting_force();
	_grid_Polymer.initialize_load_total(3);
	_grid_Polymer.calculate_load_from_CNT(&_grid_CNT, 3);
}
void Tcomposite::CG_iteration()
{
	Calculate_F0_b();
	_grid_Polymer.Solving_equilibrium_equation();
	_grid_Polymer.CG_initialize_p();
	_grid_Polymer.Calculate_Kp();
	Calculate_Fp();
	_grid_Polymer.CG_initialize_r();
	double diff = _grid_Polymer.calculate_diff();
	if (diff<1e-10)
	{
		return;
	}
	//_grid_Polymer.Calculate_Kp();
	//Calculate_Fp();
	//_grid_Polymer.Calculate_Ap();
	//_grid_Polymer.CG_update();

	//_grid_Polymer.Calculate_Kd();
	//Calculate_Fd();
	//_grid_Polymer.Calculate_Ad();
	//_grid_Polymer.calculate_rd();
	//diff = _grid_Polymer.calculate_diff();
	//cout<< diff << endl;
	int num_iteration = 0;
	while (diff>1e-15)
	{
		_grid_Polymer.Calculate_Kp();
		Calculate_Fp();
		_grid_Polymer.Calculate_Ap();
		_grid_Polymer.CG_update();
		num_iteration = num_iteration + 1;
		diff = _grid_Polymer.calculate_diff();
		cout << num_iteration<<" "<<diff << endl;
	}
	//_grid_Polymer.CG_initialize_p();
	//_grid_Polymer.Calculate_Kp();
	//Calculate_Fp();
	//_grid_Polymer.CG_initialize_r();
	//diff = _grid_Polymer.calculate_diff();
	//cout << "Final error: " << diff << endl;
	return;
}
void Tcomposite::output_result()
{
	ofstream output;
	output.open("dis.dat");
	_grid_Polymer.output_dis(output);
	_grid_Polymer.Output_Tecplot("Polymer.dat");
	_grid_CNT.calculate_CNT_boundary_displacement(&_grid_Polymer);
	_grid_CNT.calculate_load_from_constraint();
	_grid_CNT.Solving_equilibrium_equation();
	_grid_CNT.Output_Tecplot("CNT.dat");
}