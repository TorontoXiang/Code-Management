#include "TFE_grid.h"
#include <iostream>
#include "readin.h"
#include <fstream>
using namespace std;
int main()
{
	Skeyword keyword;
	ifstream input;
	input.open("Single_cell_CNT_test.k");
	read_in_keyword_file(input, keyword);
	Tgrid_CNT grid_CNT;
	Tgrid_CNT* CNT_ptr = &grid_CNT;
	grid_CNT.Input_CNT(keyword);
	grid_CNT.Create_MKL_solver();

	Tgrid_Polymer grid_Polymer;
	Tgrid_Polymer* Polymer_ptr = &grid_Polymer;
	grid_Polymer.Input_Polymer();
	grid_Polymer.Create_MKL_solver();
	grid_Polymer.calculate_load_from_constraint();
	grid_Polymer.Solving_equilibrium_equation();
	ofstream output;
	output.open("dis.dat");
	grid_Polymer.output_dis(output);
	//grid_Polymer.Output_Tecplot("Test.dat");
	int i = 0;
	grid_CNT.calculate_CNT_location(Polymer_ptr);
	while (i<30)
	{
		grid_CNT.calculate_CNT_boundary_displacement(Polymer_ptr);
		grid_CNT.calculate_load_from_constraint();
		grid_CNT.Solving_equilibrium_equation();
		grid_CNT.calculate_reacting_force();
		grid_Polymer.initialize_load_total();
		grid_Polymer.calculate_load_from_CNT(CNT_ptr);
		grid_Polymer.Solving_equilibrium_equation();
		grid_Polymer.output_dis(output);
		i = i + 1;
	}

	grid_CNT.Output_Tecplot("test_CNT.dat");

	//grid_Polymer.initialize_load_total();
	//grid_Polymer.calculate_load_from_CNT(CNT_ptr);
	//grid_Polymer.Solving_equilibrium_equation();
	grid_Polymer.Output_Tecplot("Test.dat");
	//vec3D external_force;
	//grid_Polymer.calculate_reacting_force();
	//external_force = grid_Polymer.calculate_external_load("x_min");
	system("Pause");
	return 0;
}