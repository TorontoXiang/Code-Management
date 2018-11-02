#ifndef TCELL
#define TCELL
#include "data_structure.h"
class Tcell
{
public:
	Tcell() {};
	Tcell(int id, Snode* (&node_ptr)[8], int mid);
protected:
	int _id;
	Snode* _node_ptr[8];        //The 8 node of a cell
	int _mid;                   //The material number of the cell
	bool _is_boundary;          //Whether the cell connects to the node with boundary condition

	//Local stiffness matrix calculation
	void calculate_element_stiffness(double(&Ke)[24][24],Smat& mat);
	//Calculate the cell stiffness matrix directly from material property (For Tgrid_CNT)
	void calculate_element_stiffness_simplified_modification() {};
	//Calculate the element stiffness by simplified modification of the material property
	void calculate_element_stiffness_accurate_modification() {};
	//Calculate the element stiffness accurately by considering the immersed CNT cells
	void calculate_cell_stress(Sstress* cell_stress,double(&d_cell)[8][3],Smat& mat);
	//Calculate the stress at the 8 Gauss points
	Sstress calculate_gausspiont_stress(vec3D(&DN)[8],double(&d_cell)[8][3],Smat& mat);
	//Calculate the stress at a Gauss point

	void calculate_DN_matrix(vec3D(&DN)[8], double &detJ, double xi, double eta, double zeta);
	//Calculate the dervative of shape functions with respect to x,y,z at Gauss point
	void calcualte_Jacobi(double(&Jacobi)[3][3], vec3D(&GN)[8], double xi, double eta, double zeta);
	//Calculate the GN matrix and Jacobi matrix at Gauss point
	void accumulate_Gausspoint_stiffness(double(&Ke)[24][24], vec3D(&DN)[8], double detJ, double weight,Smat& mat);
	//Accumulate the element stiffness matrix from Gauss points

	friend class Tgrid;
	friend class Tgrid_CNT;
	friend class Tgrid_Polymer;
};
#endif

