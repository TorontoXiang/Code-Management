#include "TFE_grid.h"
#include <iostream>
#include "readin.h"
#include <fstream>
#include "Tcomposite.h"
using namespace std;
void fixed_point_iteration();
void CG_iteration();
int main()
{
	//Skeyword keyword;
	//ifstream input;
	//input.open("benchmark.k");
	//read_in_keyword_file(input, keyword);
	//Tgrid_CNT grid_CNT;
	//grid_CNT.Input_CNT(keyword);
	//grid_CNT.Create_MKL_solver();
	//system("Pause");
	//grid_CNT.calculate_load_from_constraint();
	//grid_CNT.Solving_equilibrium_equation();
	//grid_CNT.Output_Tecplot("benchmark_CNT.dat",1);
	//grid_CNT.Output_Tecplot("benchmark_Polymer.dat", 0);
	//fixed_point_iteration();
	CG_iteration();
	return 0;
}
void CG_iteration()
{
	Tcomposite composite;
	composite.Input();
	composite.CG_iteration();
	composite.output_result();
	system("Pause");
}
void fixed_point_iteration()
{
	//Skeyword keyword;
	//ifstream input;
	//input.open("single_cell_move.k");
	//read_in_keyword_file(input, keyword);
	//Tgrid_CNT grid_CNT;
	//Tgrid_CNT* CNT_ptr = &grid_CNT;
	//grid_CNT.Input_CNT(keyword);
	//grid_CNT.Create_MKL_solver();

	//Tgrid_Polymer grid_Polymer;
	//Tgrid_Polymer* Polymer_ptr = &grid_Polymer;
	//grid_Polymer.Input_Polymer();
	//grid_Polymer.Create_MKL_solver();
	//grid_Polymer.calculate_load_from_constraint();
	//grid_Polymer.Solving_equilibrium_equation();
	//ofstream output;
	//output.open("dis.dat");
	////grid_Polymer.output_dis(output);
	////grid_Polymer.Output_Tecplot("Test.dat");
	//int i = 0;
	//grid_CNT.calculate_CNT_location(Polymer_ptr);
	//while (i < 30)
	//{
	//	grid_CNT.calculate_CNT_boundary_displacement(Polymer_ptr);
	//	grid_CNT.calculate_load_from_constraint();
	//	grid_CNT.Solving_equilibrium_equation();
	//	grid_CNT.calculate_reacting_force();
	//	grid_Polymer.initialize_load_total(0);
	//	grid_Polymer.calculate_load_from_CNT(CNT_ptr, 0);
	//	grid_Polymer.Solving_equilibrium_equation();
	//	grid_Polymer.output_dis(output);
	//	i = i + 1;
	//	cout << i << endl;
	//}

	//grid_CNT.Output_Tecplot("test_CNT.dat");

	////grid_Polymer.initialize_load_total();
	////grid_Polymer.calculate_load_from_CNT(CNT_ptr);
	////grid_Polymer.Solving_equilibrium_equation();
	//grid_Polymer.Output_Tecplot("test_Polymer.dat");
	////vec3D external_force;
	////grid_Polymer.calculate_reacting_force();
	////external_force = grid_Polymer.calculate_external_load("x_min");
	//system("Pause");
	//return;
}