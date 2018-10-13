#ifndef TCELL_BASE
#define TCELL_BASE
#include <vector>
#include "data_structure.h"
#include "TMAT_base.h"
using namespace std;
//Define the base class of cell
class Tcell_base
{

public:
	Tcell_base(){};
	Tcell_base(int cell_id,int nGauss);

	//void calculate_corner_force(double dt);
	//Calculate the corner force of the cell
	//void assemble_corner_force();
	////Assemble the corner force into nodal force
	//void assemble_corner_inertance();
	////Assemble the corner inertance into node

	//Access functions
	TMAT_base* G_gausspoint(int i){return _gausspoint[i];};
protected:
	int _id;       //cell id
	vector<TMAT_base*> _gausspoint;
	vector<vec3D> _corner_stress_force;   //corner force by stress

	friend class Tbody_ALEMPM;
};
#endif