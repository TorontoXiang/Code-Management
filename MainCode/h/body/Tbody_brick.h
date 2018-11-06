#ifndef TBODY_BRICK
#define TBODY_BRICK
//Define the brick body

#include "Tbody_base.h"
#include "Tcell_brick.h"
#include "Tnode.h"
#include <vector>
#include <fstream>
using namespace std;
class Tbody_brick : public Tbody_base
{
	struct Snode_stress
	{
		double volume;
		double eqs;
	};
public:
	Tbody_brick(int i, string body_name) : Tbody_base(i, body_name) {}

	virtual double calculate_time_step();
	//Calculate the time step 
	virtual void input_body(ifstream& input);
	//Input the shell body information from input data
	virtual bool corrector_step();
	//Update the DOF in the corrector step and calculate the in-body force for next time step
	virtual void calculate_final_acceleration();
	//Calculte the final acceleration by in-body and interaction force
	virtual bool predictor_step() { return false; };
	//Calculate the nodal force at n+1/2 step in the body
	//For shell element, the force at n+1/2 equals to force at n
	virtual void start_up_step();
	//Calculate the acceleration in first time step
	virtual void calculate_nodal_inertance();
	//Calculate the nodal inertance
	virtual void output_curve();
	//Output curve variables
	virtual double calculate_maximal_effective_plastic_strain();
	//Calcualte the maximal effective plastic strain
	virtual double calculate_total_internal_energy();
	//Calcualte the total internal energy
	void calculate_body_size(vec3D &coor_min, vec3D &coor_max, double &cell_edge_max, double &cell_edge_min);
	//Calculate the coordinate range and the maximal/minimal cell edge of the grid
protected:
	vector<Tcell_brick*> _cellptr_list;
	vector<Tnode*> _nodeptr_list;

	void calculate_corner_force();
	//Calculate the corner of cells
	void assemble_nodal_force();
	//Assemble corner force to nodal force and calculate the acceleration 
	void apply_external_force();
	//Apply the external force
	virtual void output_tecplot(ofstream& output, double ratio);
	//Output the tecplot file
	virtual void calculate_force_in_body();
	//Calculate the nodal acceleration caused by the body 
private:
	Sgrid<Tcell_brick, Tnode> _grid;

	template<class T1, class T2>
	friend class Tinteraction_conode;

};
#endif
