#include "Tbody.h"
#include <fstream>
#include <iostream>
using namespace std;
int main()
{
	Tbody simulation;
	ifstream input;
	input.open("RegularGridTest_tensile.k");
	simulation.input_body(input);
	simulation.input_simulation_control();
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
	double sig[6];
	simulation.calculate_average_stress(sig);
	cout << "sxx = " << sig[0] << " syy = " << sig[1] << " szz = " << sig[2] << endl;
	cout << "sxy = " << sig[3] << " sxz = " << sig[4] << " syz = " << sig[5] << endl;
	system("Pause");
	return 0;
}