#ifndef TBODY_LAGRANGIAN
#define TBODY_LAGRANGIAN
//Define the Lagrangian body
#include "Tbody_base.h"
#include "Tcell_fluid_base.h"
#include "Tcell_pure_fluid.h"
#include "Tnode_fluid.h"
#include <vector>
#include <iostream>
class Tbody_Lagrangian : public Tbody_base
{
public:
	Tbody_Lagrangian(int i,int num_material,string body_name);

	virtual double calculate_time_step();
	//Calculate the time step 
	virtual void input_body(ifstream& input);
	//Input the shell body information from input data
	virtual bool corrector_step();
	//Update the DOF in the corrector step and update the material state
	virtual void calculate_final_acceleration();
	//Calculte the final acceleration by in-body and interaction force
	virtual bool predictor_step();
	//Calculate the nodal force at n+1/2 step in the body
	virtual void calculate_nodal_inertance();
	//Calculate the nodal inertance
	virtual void output_curve();
	//Output curve variables
	virtual void output_cell_monitor();
	//Output the cell monitor
	virtual double calculate_maximal_pressure();
	//Calculate the maximal pressure in the body
	virtual double calculate_total_internal_energy();
	//Calcualte the total internal energy
	void calculate_body_size(vec3D &coor_min,vec3D &coor_max,double &cell_edge_max,double &cell_edge_min);
	//Calculate the size of the body
	void S_material_plot(int mid,string type);
	//Set the _material_output
	void plot_reconstructed_material();
	//Plot the reconstructed material information
	bool volume_check(string type);
	//check whether the volume is negative before predict and update step

protected:	
	int _num_material;                                //Number of materials in this body
	vector<Tcell_fluid_base*> _cellptr_list;      
	vector<Tnode_fluid*> _nodeptr_list;

	void set_grid_topology();
	//Set the connection of cells and nodes
	void calculate_corner_force(double dt);
	//Calculate the corner of cells
	void assemble_nodal_force();
	//Assemble corner force to nodal force and calculate the acceleration 
	virtual void output_tecplot(ofstream& output,double ratio);
	//Output the tecplot file
	virtual void output_distribution(ofstream &output);
	//Output the distributions
	virtual void calculate_force_in_body();
	//Calculate the nodal acceleration caused by the body 
	void surface_reconstrcution();
	//Reconstrcut the surface using MoF method
	void plot_material_surface(int mid,ofstream &output);
	//Polt the surface of material mid
	void plot_material_polyhedron(int mid,ofstream &output);
	//Plot the material polyhedron of material mid
	void plot_all_material_polyhedron(ofstream &output);
	//Plot all the material polyhedron
	bool volume_check(double dt);
	//Check the volume in the correct step
private:
	Sgrid<Tcell_pure_fluid,Tnode_fluid> _grid_pure_fluid;
	Sgrid<Tcell_pure_fluid,Tnode_fluid> _grid_mix_fluid;
protected:
	struct Smaterial_plot
	{
		int _mid;
		string _output_type;
		ofstream _output;
		string _file_name;

		Smaterial_plot(int body_id,int mid,string output_type);
	};
	vector<Smaterial_plot*> _material_output;

	void reset_final_initial_condition();
	//Reset the initial conditions for multi-material ALE before time integration
	void standard_test_initial_conditions();
	//Set the specific initial conditions for standard test

	template<class T1,class T2>
	friend class Tinteraction_conode;
};
#endif