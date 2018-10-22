#ifndef TBODY_ALE
#define TBODY_ALE
//Define the body for ALE method
#include "Tcell_fluid_base.h"
#include "Tbody_Lagrangian.h"
#include "Tnode_fluid.h"
#include "Tbucket_searching.h"
#include <vector>
#include "Tvolume_equation_solver.h"
#include <iostream>
class Tgrid_smooth_base;
class Tgrid_smooth_OST;
class Tgrid_smooth_CST;
class Tgrid_smooth_PST;
class Tbody_ALE : public Tbody_Lagrangian
{
protected:
	struct Sfluid_grid
	{
		vector<Tcell_fluid_base*> _cell_list;
		vector<Tnode_fluid*> _node_list;
	};
public:
	Tbody_ALE(int i,int num_material,string body_name);

	virtual void input_body(ifstream& input);
	//Input the ALE body information from input data
	virtual double calculate_time_step();
	//Calculate the time step in ALE frame
	void input_immersed_material();
	//Input the multi-material information
	virtual bool corrector_step();
	//The variables may be interpolated to a new grid if the remapping phase occur
	virtual bool is_remapping();
	//Determine whether the remapping phase is needed to perform
	virtual void remapping_variables();
	//Reamp the variables from old grid to new grid
	void calculate_old_grid_size(vec3D &coor_min,vec3D &coor_max,double &cell_edge_max);
	//Calculate the coordinate range and the maximal/minimal cell edge of the old grid
	virtual void overlapping_with_material_grid(Tmaterial_grid* material_grid_ptr,int mid);
	//Overlapping the ALE grid with the immersed material mid grid


	//Access functions
	virtual int G_num_material(){return _num_material;};

	//Value functions
	virtual void S_remapping_scheme(Sremapping_scheme &remapping_scheme);
	virtual void S_remesh_scheme(string name);
	virtual void S_remapping_start_time(double t){_remapping_start_time=t;};
	//-----------------------------------------------------------------------
	//Functions for test
	//-----------------------------------------------------------------------
	virtual void set_velocity(int type);
	void update_position(double dt);
	void set_variables();
	void set_variables(int mid);
	void iteration_test();
	void plot_grid_new();
	void set_new_grid(int type,double dt);
	virtual void set_displacemant_at_conode(double dl);
	virtual void test_for_fraction_below();
	virtual void test_for_volume_solver();
	virtual void test_for_MoF_derivative();
	virtual void test_for_MoF_patch_test();
	virtual void test_for_MoF_patch_test_3_material();
	virtual void test_for_MoF();
	virtual void test_for_MM_remapping();

protected:
	Tgrid_smooth_base* _grid_smooth;
	Sremapping_scheme _remapping_shceme;

	void reset_new_grid_geometry();
	//Reset the geometry information after the new grid is generated
	void update_old_grid_geometry();
	//Prepare the old grid before remapping
	void old_grid_reconstruction();
	//Reconstruct the gradient of variables of old grid
	void calculate_gradient();
	//Calculate the gradient of density and energy
protected:
	//The two grid
	Sfluid_grid _grid1;
	Sfluid_grid _grid2;
	//The new and old grid pointer
	Sfluid_grid* _grid_old;
	Sfluid_grid* _grid_new;
	ofstream _iremapping;        //Output the remapping information
	double _remapping_start_time;//Time for remapping beginning


	friend class Tgrid_smooth_base;
	friend class Tgrid_smooth_origin;
	friend class Tgrid_smooth_continum;
	friend class Tgrid_smooth_OST;
	friend class Tgrid_smooth_CST;
	friend class Tgrid_smooth_PST;
	template<class T1, class T2>
	friend class Tinteraction_conode;
};
#endif