#ifndef TCELL_8H
#define TCELL_8H
#include "Tnode.h"
#include "Tcell_base.h"
#include "TMAT_base.h"
#include "data_structure.h"
//Define the 8H element for grid smooth and brick solid simulation
class Tcell_8H : public Tcell_base
{
public:
	Tcell_8H(){};
	Tcell_8H(int cell_id,int nGauss,Tnode* (&node_ptr)[8],TMAT_base** mat);
	//nGauss is the number of Gausspoint at one dimension

	void calculate_element_stiffness(double (&Ke)[24][24]);
	//Calculate the element stiffness matrix and the corresponding freedom degrees
	void modify_Youngs();
	//Modify the Youngs model by the cell distorted quantity

private:
	Tnode* _node_ptr[8];        //The connected node information
	int _nGauss;                //The number of Gausspoint at one dimension
	double _distance_to_conode; //The distance of the cell to the co-node surface

	void calculate_elastic_matrix(double &D0,double &k1,double &k2);
	//Calculate the elastic matrix of the element
	void calculate_DN_matrix(vec3D (&DN)[8],double &detJ,double xi,double eta,double zeta);
	//Calculate the dervative of shape functions with respect to x,y,z at Gauss point
	void calcualte_Jacobi(double (&Jacobi)[3][3],vec3D (&GN)[8],double xi,double eta,double zeta,int type);
	//Calculate the GN matrix and Jacobi matrix at Gauss point
	//type=0:Calculate the Jacobi matrix with respect to the initial position
	//type=1:Calculate the Jacobi matrix with respect ro the current position
	void accumulate_Gausspoint_stiffness(double (&Ke)[24][24],vec3D (&DN)[8],double detJ,double weight);
	//Accumulate the element stiffness matrix from Gauss points
	double calculate_distorted_quantity();
	//Calculate the distorted quantity of the cell


	friend class Tgrid_smooth_continum;
	friend class Tgrid_smooth_continum_new;

};
#endif