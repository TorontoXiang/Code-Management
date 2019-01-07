#ifndef TCELL
#define TCELL
#include "data_structure.h"
struct Snode
{
	Snode() {};
	Snode(int id, double(&coor)[3])
	{
		_id = id; _pos[0] = coor[0]; _pos[1] = coor[1]; _pos[2] = coor[2];
		_bc[0] = _bc[1] = _bc[2] = 0; _bc_value[0] = _bc_value[1] = _bc_value[2] = 0;
		_is_surface = false;
	};

	int _id;
	double _pos[3];
	int _bc[3];                //Boundary condition  of a node
	double _bc_value[3];       //The value of boundary condition
	bool _is_surface;          //Whether a node is a surface node
	Slocation_info* _location; //The location information of a CNT node (Only allocate memory for the surface node in CNT grid)

	void calculate_location(double(&x_min)[3], double(&x_max)[3], double(&interval)[3], int nx_max, int ny_max, int nz_max);
	//Calculate the location of a surface node in CNT grid

};
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
	friend class Tgrid_CNT_base;
	friend class Tgrid_CNT_T;
};
#endif

