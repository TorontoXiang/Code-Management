#ifndef TBODY_ALEMPM
#define TBODY_ALEMPM
//Define the body for couple ALE-MPM method
#include "Tbody_ALE.h"
#include "Tcell_fluid_base.h"
#include "Tnode_fluid.h"
#include "Tbucket_searching.h"
#include <vector>
#include "Tvolume_equation_solver.h"
#include <iostream>
#include "TMPM_particle.h"
#include "Tcell_mixed_fluid_MPM.h"
#include "Tcell_pure_fluid_MPM.h"
#include <Windows.h>
#include "Tcell_intersect.h"
#include "Tgrid_smooth.h"
#include <omp.h>
#include <iomanip>
class Tbody_ALEMPM : public Tbody_ALE
{
public:
	Tbody_ALEMPM(int i,int num_material,string body_name);
	virtual void input_body(ifstream& input);
	//Input the ALEMPM body
	virtual double calculate_time_step();
	//Calculate the time step in the ALEMPM body
	virtual bool predictor_step();
	//Calculate the corner nodal force from fluid and particles
	virtual bool corrector_step();
	//The corrector step in ALEMPM
	virtual void calculate_nodal_inertance();
	//Calculate the inertance of the node
	virtual void output_tecplot(ofstream& output,double ratio);
	//Output the tecplot file

protected:
	vec3D _x_min,_x_max;
	int _nx,_ny,_nz;
	vec3D _interval;
	int _num_particle_part;                       //The number of particle part in the body
	vector<TMPM_particle>* _particle_group_list;  //The particle group list
	ofstream output_MPM;                          //The file stream for MPM particles
	ofstream output_temp;

	void generate_back_ground(vec3D x_min,vec3D x_max,int nx,int ny, int nz);
	//Generate the back ground for MPM and the computational grid for fluid
	double calculate_time_step_particle();
	//Calculate the time step of particles
	void calculate_particle_force();
	//Calculate the corner force from particles
	int find_particle_location(TMPM_particle* particle_ptr);
	//Find which cell the particle located in
	void MPM_tec_output();
	//Output the particles in tecplot format
	void remapping_variables();
	//Remap the variable from old grid to the new grid

	friend class Tgrid_smooth_base;
	friend class Tgrid_smooth_origin;
	friend class Tgrid_smooth_continum;
};
#endif