#ifndef TBODY_EXPLOSION
#define TBODY_EXPLOSION
#include "Tbody_ALE.h"
#include <fstream>
#include "Tcell_fluid_base.h"
#include "Tnode_fluid.h"
#include "Tbucket_searching.h"
#include <vector>
#include "Tvolume_equation_solver.h"
#include "Tbody_ALEMPM.h"
#include <iostream>
class Tgrid_smooth_base;
class Tbody_explosion : public Tbody_ALEMPM
{
public:
	Tbody_explosion(int i,int num_material,string body_name);

	virtual void input_body(ifstream& input);
	//Input the explosion body
	virtual void remapping_variables();
	//Reamp the variables from old grid to new grid
	virtual bool is_remapping();
	//Whether remapping in the explosion body
	virtual void S_detonation_time(double detonation_time){_detonation_time=detonation_time;};
	virtual void S_remapping_scheme(Sremapping_scheme &remapping_scheme);
protected:
	Sfluid_grid _grid_initial;
	int _nume_initial,_nume_com;
	int _nump_initial,_nump_com;
//	vec3D _x_min,_x_max;
//	int _nx,_ny,_nz;
//	vec3D _interval;
	int _num_remapping;
	double _detonation_time;
	int _num_deactive_particle_part;

	void generate_back_ground(vec3D x_min,vec3D x_max,int nx,int ny, int nz);
	void generate_cylinder_particles(double r,double z_min,double z_max,vec2D center);
};
#endif