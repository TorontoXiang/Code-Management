#ifndef TGRID_SMOOTH
#define TGRID_SMOOTH
//Define the base class for grid smooth
#include "Tbody_ALE.h"
#include <string>
using namespace std;
class Tgrid_smooth_base
{ 
public:
	Tgrid_smooth_base(Tbody_ALE* the_body):_ALE_body(the_body){};

	//Generate a new grid in this_body
	virtual void generate_new_grid()=0;
protected:
	string _smooth_name;
	Tbody_ALE* _ALE_body;
};
//Smooth the grid old into its origin position
class Tgrid_smooth_origin : public Tgrid_smooth_base
{
public:
	Tgrid_smooth_origin(Tbody_ALE* the_body):Tgrid_smooth_base(the_body){_smooth_name="Origin";};

	//Generate a new grid in this_body
	virtual void generate_new_grid()
	{
		for (int i = 0; i < _ALE_body->_nump; i++)
		{
			_ALE_body->_grid_new->_node_list[i]->reset_position();
		}
		return;
	}
};
//Smooth the grid by orthogonality similarity transformation
class Tgrid_smooth_OST : public Tgrid_smooth_base
{
public:
	Tgrid_smooth_OST(Tbody_ALE* the_body);

	//Generate a new grid in this_body
	virtual void generate_new_grid();
protected:
	double _x_min0,_y_min0,_z_min0;
	double _x_max0,_y_max0,_z_max0;
};
//Sommth the grid by cylinder similarity transformation
class Tgrid_smooth_CST : public Tgrid_smooth_base
{
public:
	Tgrid_smooth_CST(Tbody_ALE* the_body,vec2D center);

	//Generate a new grid in this_body
	virtual void generate_new_grid();
protected:
	vec2D _center;
	double _r0;
};
//Smooth the grid by sphere similary transforamtion
class Tgrid_smooth_PST : public Tgrid_smooth_base
{
public:
	Tgrid_smooth_PST(Tbody_ALE* the_body,vec3D center);

	//Generate a new grid in this_body
	virtual void generate_new_grid();
protected:
	vec3D _center;
	double _r0;
};
#endif