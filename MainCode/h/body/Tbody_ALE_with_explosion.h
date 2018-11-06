#ifndef TBODY_ALE_WITH_EXPLOSION
#define TBODY_ALE_WITH_EXPLOSION
//Define the body for ALE method
#include "Tcell_fluid_base.h"
#include "Tbody_ALE.h"
#include "Tnode_fluid.h"
#include "Tbucket_searching.h"
#include <vector>
#include "Tvolume_equation_solver.h"
#include <iostream>
class Tgrid_smooth_base;
class Tgrid_smooth_OST;
class Tgrid_smooth_CST;
class Tgrid_smooth_PST;
class Tbody_ALE_with_explosion : public Tbody_ALE
{
protected:
	struct Sfluid_grid
	{
		vector<Tcell_fluid_base*> _cell_list;
		vector<Tnode_fluid*> _node_list;
	};
public:
	Tbody_ALE_with_explosion(int i, int num_material, string body_name);

	virtual void input_body(ifstream& input);


protected:
	Sfluid_grid* input_explosion_region();
	//Input the information of the explosion region
	void overlapping_with_explosion_region();
	//Overlapping the initial grid with the compulation region
	void calculate_explosion_region_size(Sfluid_grid* grid_exp,vec3D& coor_min, vec3D& coor_max, double cell_edge_max);
	//Calculate the size of the explosion region

	friend class Tgrid_smooth_base;
	friend class Tgrid_smooth_origin;
	friend class Tgrid_smooth_continum;
	friend class Tgrid_smooth_OST;
	friend class Tgrid_smooth_CST;
	friend class Tgrid_smooth_PST;
	template<class T1, class T2>
	friend class Tinteraction_conode;
};
#endif#pragma once
