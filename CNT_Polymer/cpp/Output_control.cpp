#include "Output_control.h"
#include <iomanip>
using namespace std;
Smesh_output::Smesh_output(double amplify_ratio, int total_frame, double endtime) :
	_amplify_ratio(amplify_ratio), _total_frame(total_frame)
{
	char a_bodyid[32];
	string s_bodyid, file_name;
	file_name = "Tecplot_mesh_initial.dat";
	output_initial.open(file_name);
	file_name = "Tecplot_mesh_process.dat";
	output_process.open(file_name);
	file_name = "Tecplot_mesh_final.dat";
	output_final.open(file_name);
	_output_point.resize(_total_frame - 1);
	double dt = endtime / _total_frame;
	for (int i = 0; i < _total_frame - 1; i++)
	{
		_output_point[i] = dt * (i + 1);
	}
}
bool Smesh_output::is_output(double current_time, double dt_pre)
{
	for (int i = 0; i < _total_frame - 1; i++)
	{
		if (current_time > _output_point[i])
		{
			if (current_time - dt_pre < _output_point[i])
			{
				_output_point[i] = -1;
				return true;
			}
		}
	}
	return false;
}
Scurve_output::Scurve_output(int type, int id): _type(type), _id(id)
{
	char a_id[32];
	string s_id, file_name;
	_itoa_s(_id, a_id, 10);
	s_id = a_id;
	if (type == 0)
	{
		file_name = "Curve_for_Node" + s_id + ".dat";
		output.open(file_name);
		output << "Time-history curve for node " << _id << endl;
		output << setw(20) << "Time" << setw(20) << "x-displacement" << setw(20) << "y-displacement" << setw(20) << "z-displacement" << endl;
	}
	else if (type == 1)
	{
		file_name = "Curve_for_Cell" + s_id + ".dat";
		output.open(file_name);
		output << "Time-history curve for cell " << _id << endl;
	}
	return;
}