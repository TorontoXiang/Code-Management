#ifndef  TGRID_AND_CELL
#define TGRID_AND_CELL
#include "Tpolygon.h"
#include "public.h"
#include <vector>
#include <iostream>
using namespace std;
struct Tcell
{
	Tcell(){};
	int IEN[4];
	double fraction[2];
	vec2D centroid[2];
	Tpolygon polygon;
	Tpolygon* material_polygon[2];
};

struct Tgrid
{	
	int nume,nump;
	vector<vec2D> node_list;
	vector<Tcell> cell_list;

	void plot_material_polygon(ofstream& output,int id);
	void plot_grid(ofstream& output);
};
struct Tcell_polygon
{
	vector<vec2D*> _node_list;
	Tpolygon _cell_polygon;
	double _fraction[2];
	vec2D _centroid[2];
	Tpolygon* _material_polygon[2];

	void generate_cell_polygon();
	void intersecting_with_other(Tcell* other);
};
struct Tgrid_polygon
{
	vec2D** _node_list;
	vector<Tcell_polygon> _cell_list;
	int _nx_node,_ny_node;

	void generate_node_list();
	void generate_cell_list();
	void remapping_surface(Tgrid* grid_old);

	void plot_grid(ofstream& output,int id);
	void plot_material(ofstream& output,int id);
};
#endif // ! TGRID_AND_CELL
