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
	//_grid_Polymer.calculate_load_from_constraint();
	_grid_CNT.calculate_CNT_location(&_grid_Polymer);
	//system("Pause");
}
//void Tcomposite::Calculate_F0_b()
//{
//	//Using zero Polymer displacement to calculate the constant term of _Fp
//	_grid_CNT.calculate_CNT_boundary_displacement(&_grid_Polymer);
//	_grid_CNT.calculate_load_from_constraint();
//	_grid_CNT.Solving_equilibrium_equation();
//	_grid_CNT.calculate_reacting_force();
//	_grid_Polymer.initialize_load_total(1);
//	_grid_Polymer.calculate_load_from_CNT(&_grid_CNT,1);
//	_grid_Polymer.Calculate_b();
//}
void Tcomposite::Calculate_Fp()
{
	_grid_CNT.calculate_CNT_boundary_p(&_grid_Polymer);
	_grid_CNT.calculate_load_from_constraint();
	_grid_CNT.Solving_equilibrium_equation();
	_grid_CNT.calculate_reacting_force();
	_grid_Polymer.initialize_Fp();
	_grid_Polymer.assemble_Fp(&_grid_CNT);
}
//void Tcomposite::Calculate_Fd()
//{
//	_grid_CNT.calculate_CNT_boundary_displacement(&_grid_Polymer);
//	_grid_CNT.calculate_load_from_constraint();
//	_grid_CNT.Solving_equilibrium_equation();
//	_grid_CNT.calculate_reacting_force();
//	_grid_Polymer.initialize_load_total(3);
//	_grid_Polymer.calculate_load_from_CNT(&_grid_CNT, 3);
//}
void Tcomposite::CG_iteration()
{
	CG_iteration_initialization();
	//_grid_Polymer.CG_initialization();
	//Calculate_F0_b();
	//_grid_Polymer.Solving_equilibrium_equation();
	//_grid_Polymer.CG_initialize_p();
	//_grid_Polymer.Calculate_Kp();
	//Calculate_Fp();
	//_grid_Polymer.CG_initialize_r();
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

		//_grid_Polymer.Calculate_Ap();
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