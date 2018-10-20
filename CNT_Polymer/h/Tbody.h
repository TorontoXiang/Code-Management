#ifndef TBODY
#define TBODY
//Define the brick body
#include "readin.h"
#include "Tcell_brick.h"
#include "Tnode.h"
#include <vector>
#include <fstream>
#include "Output_control.h"
using namespace std;
class Tbody
{
public:
	Tbody() {};

	void calculate_force_in_body();
	//Calculate the nodal acceleration caused by the body 
	void calculate_time_step();
	//Calculate the time step 
	void input_body(ifstream& input);
	//Input the shell body information from input data
	void input_simulation_control();
	//Input the control parameters of the simulation
	void update_grid();
	//Update the DOF in the corrector step and calculate the in-body force for next time step
	void calculate_final_acceleration();
	//Calculte the final acceleration by in-body and interaction force
	bool predictor_step() { return false; };
	//Calculate the nodal force at n+1/2 step in the body
	//For shell element, the force at n+1/2 equals to force at n
	void start_up_step();
	//Calculate the acceleration in first time step
	void calculate_nodal_inertance();
	//Calculate the nodal inertance
	double calculate_total_internal_energy();
	//Calcualte the total internal energy
	bool is_finish() { return _current_time == _endtime; };
	//Determine whether the simulation is finished
	void advance_in_time() { _current_time = _current_time + _dt; _num_step = _num_step + 1; };
	//Update the current time
	void Output_Tecplot_mesh(string type);
	//Output the Tecplot mesh
	void Output_Curve();
	//Output the history curve for specific nodes and cells

	//Access functions
	double G_current_time() { return _current_time; };
	int G_step() { return _num_step; };
	double G_dt() { return _dt; };
protected:
	struct Sexternal_load
	{
		int num_node;           //The number of node in this load group
		vector<int> node_id;    //The node id in this load group
		int direction;          //The direction of this load group
		double magnitude;       //The magnitude of this load group
	};
	struct Sbrick_grid
	{
		vector<Tcell_brick> _cell_list;
		vector<Tnode> _node_list;
	};

	int _nume;              //Number of element
	int _nump;              //Number of node
	int _nx, _ny, _nz;      //Number of cell in x,y,z direction (Only for regular grid)
	vec3D _x_min, _x_max;   //Extreme coordinate of the regular grid (Only for regular grid)
	vec3D _interval;        //The coordinate increment at x,y,z direction (Only for regular grid)
	int _num_step;          //Simulation steps
	double _dt,_dt_pre;     //Current and previous time increment
	double _endtime;        //The end time
	double _current_time;   //The current time
	double _CFL;            //The CFL number
	vec3D _gravity;         //The gravity in this body
	vector<Sexternal_load> _external_load_list;     //The external load list

	Sbrick_grid _grid_polymer;
	vector<Tcell_brick*> _cellptr_list;
	vector<Tnode*> _nodeptr_list;

	vector<Scurve_output*> _curve_out_list;   //Output curve
	Smesh_output* _mesh_out;                  //Output Tecplot mesh

	void calculate_corner_force();
	//Calculate the corner of cells
	void assemble_nodal_force();
	//Assemble corner force to nodal force and calculate the acceleration 
	void apply_external_force();
	//Apply the external force
	void output_tecplot(ofstream& output, double ratio);
	//Output the tecplot file
	void create_curve(int type, int id);
	//Create a Scurve_output to output curve
	void create_Tecplot_mesh(double amplify_ratio, int total_frame);
	//Create the tecplot mesh output
	void input_irregular_grid(Skeyword& keyword);
	//Input a irregular grid
	void input_regular_grid(Skeyword& keyword);
	//Input a regular grid
	void apply_regular_bc(Sregular_grid_bc& bc);
	//Apply the boundary condition for a regular grid


};
#endif
