#ifndef TBODY_BASE
#define TBODY_BASE
#include <vector>
#include "data_structure.h"
#include <fstream>
#include "readin.h"
#include <stdlib.h>
#include <string>
#include "Tviscosity_base.h"
#include "TMoF_solver.h"
#include "Tmaterial_grid.h"
using namespace std;
//Define the base class of body
class Tbody_base
{
public:
	template<class T1,class T2>
	struct Sgrid
	{
		vector<T1> _cell_list;
		vector<T2> _node_list;
	};
public:
	Tbody_base(int i,string body_name);
//------------------------------------------------------------------------
//Virtual functions
//------------------------------------------------------------------------
	virtual double calculate_time_step(){return 0;};
	//Calculate the time step 
	virtual void input_body(ifstream& input){};
	//Input the body information from input data
	virtual bool corrector_step(){return false;};
	//Update the DOF and material state
	virtual void calculate_final_acceleration(){};
	//Calculate the final acceleration (inclduing interaction)
	virtual bool predictor_step(){return false;};
	//Calculate the nodal force at n+1/2 step in the body
	virtual void start_up_step(){};
	//Calculate the in-body force in first time step
	virtual void calculate_nodal_inertance(){};
	//Calculate the nodal inertance
	virtual void output_curve(){};
	//Output curve variables
	virtual void output_cell_monitor(){};
	//Output the cell monitor
	virtual double calculate_maximal_pressure(){return 0;};
	//Calculate the maximal pressure in the body
	virtual double calculate_maximal_effective_plastic_strain(){return 0;};
	//Calcualte the maximal effective plastic strain
	virtual double calculate_total_internal_energy(){return 0;};
	//Calcualte the total internal energy
	virtual void overlapping_with_material_grid(Tmaterial_grid* material_grid_ptr,int mid){};
	//Overlapping the ALE grid with the immersed material mid grid
	virtual void standard_test_initial_conditions(){};
	//Set the specific initial conditions for standard test
	virtual void output_k_file() {};
	//Output the k file of the explosion region
//-----------------------------------------------------------------
//Non-virtual functions
//-----------------------------------------------------------------
	void create_curve(int body_id,int type,int id);
	//Create a Scurve_output to output curve
	void create_mesh(int body_id,double amplify_ratio,int interval_step);
	//Create a Smesh_output to output mesh
	void create_cell_monitor(int body_id,int cell_id);
	//Create a Scell_monitor to output cell deformation
	void advance_in_time();
	//Advance in time
	void output_mesh(string type);
	//Output the mesh
	//initial: output the mesh of the initial phase
	//process: output the mesh in the process
	//final: output the mesh of the final phase
	void output_global_curve();
	//Output the global variables
	void exam_av_type();
	//Exam whether artificial type is valide
	string calculate_body_type_id();
	//Calculate the body type id

	//Access function
	double G_dt(){return _dt;};
	virtual int G_num_material(){return 0;};

	//Value functions
	void S_dt(double dt_global){_dt=dt_global;};
	virtual void S_remapping_scheme(Sremapping_scheme &remapping_scheme){};
	virtual void S_remesh_scheme(string name){};
	virtual void S_material_plot(int mid,string type){};
	virtual void S_detonation_time(double detonation_time){};
	virtual void S_remapping_start_time(double t){};
	//-------------------------------------------------------------------------
	//Fucntions for test
	//-------------------------------------------------------------------------
	virtual void set_velocity(int type){};
	virtual void set_displacemant_at_conode(double dl){};
	virtual void test_for_fraction_below(){};
	virtual void test_for_volume_solver(){};
	virtual void test_for_MoF_derivative();
	virtual void test_for_MoF_patch_test();
	virtual void test_for_MoF_patch_test_3_material();
	virtual void test_for_MoF(){};
	virtual void test_for_MM_remapping(){};
protected:
	struct Sexternal_load
	{
		int num_node;           //The number of node in this load group
		vector<int> node_id;    //The node id in this load group
		int direction;          //The direction of this load group
		double magnitude;       //The magnitude of this load group
	};
	struct Smesh_output
	{
		int _body_id;
		double _amplify_ratio;   //The amplify ratio in mesh plot
		int _total_frame;        //The step interval for mesh plot
		ofstream output_initial;
		ofstream output_process;
		ofstream output_final;
		ofstream distribution;
		string file_name;
		vector<double> _output_point;

		Smesh_output(int body_id,double amplify_ratio,int total_frame,double endtime);
		Smesh_output(){};
		bool is_output(double current_time,double dt);
	};
	struct Scurve_output
	{
		int _body_id;
		int _type;               //0-node curve;1-cell curve; 
		int _id;                 //Node or cell id
		ofstream output;    
		string file_name;    

		Scurve_output(int body_id,int type,int id);
	};
	struct Scell_monitor
	{
		int _body_id;
		int _cell_id;
		ofstream output;
		string file_name;

		Scell_monitor(int body_id,int cell_id);
	};
	struct Sbody_global_variable
	{
		int _body_id;
		bool _is_output;
		ofstream output;
		string file_name;
		bool _internal_energy,_max_p,_epeqv;
		Sbody_global_variable(int body_id,bool is_output);
		void create_ofstream();
	};
	string _body_name;
	int _id;
	double _dt;          //Time step in this body
	double _endtime;     //The end time
	double _current_time;//The current time
	double _CFL;         //The CFL number
	int _current_step;   //The current step
	int _nume;            //The number of cell
	int _nump;            //The number of node
	Tviscosity_base* _Arti_Vis;    //The artificial viscosity in this body
	vec3D _gravity;       //The gravity in this body
	string _test_name;    //Whether it is a standard test
	int _dimension;       //The dimension of this problem
	bool _hourglass_option;

	vector<Sexternal_load> _external_load_list;     //The external load list

	vector<Scurve_output*> _curve_out_list;   //Files for output curve
	vector<Scell_monitor*> _cell_monitor_list; //Files for cell monitor

	//Body global variables
	Sbody_global_variable* _global_output;
public:
	Smesh_output* _mesh_out;             //Files for output mesh

	virtual void output_tecplot(ofstream& output,double ratio){};
	//Output the tecplot file
	virtual void output_distribution(ofstream &output){};
	//Output the variables distributions
	virtual void calculate_force_in_body(){};
	//Calculate the nodal force in the body
	virtual void modify_body(ifstream &input){};

	friend class Tregion;
};
#endif