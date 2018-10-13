#ifndef TCELL_INTERSECT
#define TCELL_INTERSECT
//Define the cell intersection algorithm
#include "Tcell_fluid_base.h"
#include "Ttetrahedron.h"
#include "Tpolyhedron.h"
class Tcell_intersect
{
public:
	Tcell_intersect():_new_cell(NULL),_old_cell(NULL){};
	//Tcell_intersect_new(int num_material):_new_cell(NULL),_old_cell(NULL),_num_material(num_material){};

	void intersect();
	//Intersect new and old cell and give the intersection portion to the new cell (by tetrahedron intersecting with a polyhedron)
	void overlapping(int mid);
	//Overlap the old grid with the material i th polyhedron
	void create_tet_list();
	//Create the tetrahedron list of the new cell
	//Initialize the transform matrix and volume of every tetrahedrons

	//Value functions
	void S_new_cell(Tcell_fluid_base* new_cell){_new_cell=new_cell;};
	void S_old_cell(Tcell_fluid_base* old_cell){_old_cell=old_cell;};
	void S_material_poly(Ttetrahedron* tet_ptr){_material_tet=tet_ptr;};

private:
	//int _num_material;                 //The number of materials in cell
	Tcell_fluid_base* _new_cell;       //Cell in new grid
	Tcell_fluid_base* _old_cell;       //Cell in old grid
	Ttetrahedron* _material_tet;       //The material polyhedron in overlapping process

	Ttetrahedron _tet_list[48];        //Tetrahedron list in new cell
	int _ass_node[48];                 //The corresponding local node id of every tetrahedrons in tet_list
	Tpolyhedron _corner_polyhedron[8];        

	void create_corner_polyhedron();
	//Create the corner polyhedron for cell change when the velocity is discontinuous
	bool whether_dis(double tolerance);
	//Determine whether the velocity of cell_change is continuous with respect to the tolerance
};
#endif