#ifndef TMPM_BACKGROUND
#define TMPM_BACKGROUND
//Define the class of particle in MPM background
#include <vector>
#include "data_structure.h"
#include "Tcell_MPM.h"
#include "TMPM_particle.h"
#include "Tnode.h"
#include <fstream>
using namespace std;
class TMPM_background
{
public:
	void generate_background(vec3D& x_min,vec3D& x_max,int nx,int ny,int nz);
	//Generate a background
	void initialize();
	//Initialize the variables on the background grid
	void map_particle_to_cell(TMPM_particle* particle_ptr);
	//Distribute the particle mass and moment to the background node
	void calculate_character_length();
	//Calculate the character length of each cell
	void distribute_particle_force(TMPM_particle* particle_ptr);
	//Distribute the particle force to the background node
	void update_background_DOF(double dt);
	//Update the DOF of the background
	void update_particle_DOF(TMPM_particle* particle_ptr,double dt);
	//Update the DOF of the particle
	void update_particle_state(TMPM_particle* particle_ptr,double dt);
	//Update the stress on particle_ptr
	void assemble_nodal_force();
	//Assemble the nodal force
	void reset_nodal_velocity();
	//Reset the node velocity
	void apply_background_load();
	//Apply the load on the background grid
	void plot_background(ofstream& output);
	//Polt the background
	void plot_background_cell(ofstream& output,int cell_id);

	vector<Tnode> _node_list;           //The node list of the background grid
	vector<Tcell_MPM> _cell_list;       //The cell list of the background grid
	vec3D _x_min,_x_max;                //The range of the background
	vec3D _interval;                    //The interval of each dimension
	int _nx,_ny,_nz;                    //The division of each dimension
	int _total_node;                    //The total number of node in the background
	int _total_cell;                    //The total number of cell in the background 

	int calculate_paritcle_loaction(TMPM_particle* particle_ptr);
	//Calculate the cell index of a particle
};
#endif