#include "Tbody.h"
#include <fstream>
#include <iostream>
using namespace std;
int main()
{
	Tbody simulation;
	//Input grid
	ifstream input_Polymer, input_CNT;
	input_Polymer.open("RegularGridTest_tensile.k");
	simulation.input_Polymer(input_Polymer);
	input_CNT.open("horizon_bar_x_y_incliened_50.k");
	simulation.input_CNT(input_CNT);
	simulation.input_simulation_control();
	//Simulation
	simulation.calculate_nodal_inertance();
	simulation.Output_Tecplot_mesh("Initial");
	while (!simulation.is_finish())
	{
		simulation.calculate_time_step();
		simulation.calculate_force_in_body();
		simulation.calculate_final_acceleration();
		simulation.update_grid();
		simulation.advance_in_time();
		simulation.change_state();

		//Output results in screen and files
		cout << simulation.G_step()<<" step with t= "<<simulation.G_current_time() << " dt= "<<simulation.G_dt()<< endl;
		simulation.Output_Tecplot_mesh("Process");
		simulation.Output_Curve();
		//system("Pause");
	}
	simulation.Output_Tecplot_mesh("Final");
	simulation.output_CNT_stress();
	simulation.output_Polymer_stress();
	double sig[6];
	simulation.calculate_average_stress(sig);
	cout << "sxx = " << sig[0] << " syy = " << sig[1] << " szz = " << sig[2] << endl;
	cout << "sxy = " << sig[3] << " sxz = " << sig[4] << " syz = " << sig[5] << endl;
	system("Pause");
	return 0;
}