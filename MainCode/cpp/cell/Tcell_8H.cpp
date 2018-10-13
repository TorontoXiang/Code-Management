#include "Tcell_8H.h"
#include "Shape_function.h"
#include "public_function.h"
#include "Gauss_weight.h"
#include "public_function.h"
#include "Ttetrahedron.h"
Tcell_8H::Tcell_8H(int cell_id,int nGauss,Tnode* (&node_ptr)[8],TMAT_base** mat)
{
	_id=cell_id;
	_nGauss=nGauss;
	_gausspoint.resize(_nGauss*_nGauss*_nGauss);  
	//Allocate memory
	_corner_stress_force.resize(8);
	//Set cell node
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]=node_ptr[i];
	}
	//Set Gauss point material
	for (int i = 0; i < _nGauss*_nGauss*_nGauss; i++)
	{
		_gausspoint[i]=mat[i];
	}
	vec3D cell_center;
	for (int i = 0; i < 8; i++)
	{
		cell_center=cell_center+_node_ptr[i]->_position0;
	}
	cell_center=cell_center/8.0;
	//_distance_to_conode=cell_center.get_length()-1;
	_distance_to_conode=1;
}
void Tcell_8H::calculate_element_stiffness(double (&Ke)[24][24])
{
	//Initialize the element stiffness
	for (int i = 0; i < 24; i++)
	{
		for (int j = 0; j < 24; j++)
		{
			Ke[i][j]=0;
		}
	}
	//Calculate the element stiffness
	double xi,eta,zeta,detJ,weight;
	vec3D DN[8];
	for (int i = 0; i <_nGauss; i++)
	{
		for (int j = 0; j < _nGauss; j++)
		{
			for (int k = 0; k < _nGauss; k++)
			{
				xi=Gauss_coordinate[_nGauss-1][i],eta=Gauss_coordinate[_nGauss-1][j],zeta=Gauss_coordinate[_nGauss-1][k];
				weight=Gauss_weight[_nGauss-1][i]*Gauss_weight[_nGauss-1][j]*Gauss_weight[_nGauss-1][k];
				calculate_DN_matrix(DN,detJ,xi,eta,zeta);
				accumulate_Gausspoint_stiffness(Ke,DN,detJ,weight);
			}
		}
	}
	return;
}
void Tcell_8H::accumulate_Gausspoint_stiffness(double (&Ke)[24][24],vec3D (&DN)[8],double detJ,double weight)
{
	double D0,k1,k2;
	calculate_elastic_matrix(D0,k1,k2);
	double BiDBj[3][3];
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			BiDBj[0][0]=DN[i].x*DN[j].x+k2*DN[i].y*DN[j].y+k2*DN[i].z*DN[j].z;
			BiDBj[1][0]=k1*DN[i].y*DN[j].x+k2*DN[i].x*DN[j].y;
			BiDBj[2][0]=k1*DN[i].z*DN[j].x+k2*DN[i].x*DN[j].z;
			BiDBj[0][1]=k1*DN[i].x*DN[j].y+k2*DN[i].y*DN[j].x;
			BiDBj[1][1]=DN[i].y*DN[j].y+k2*DN[i].x*DN[j].x+k2*DN[i].z*DN[j].z;
			BiDBj[2][1]=k1*DN[i].z*DN[j].y+k2*DN[i].y*DN[j].z;
			BiDBj[0][2]=k1*DN[i].x*DN[j].z+k2*DN[i].z*DN[j].x;
			BiDBj[1][2]=k1*DN[i].y*DN[j].z+k2*DN[i].z*DN[j].y;
			BiDBj[2][2]=DN[i].z*DN[j].z+k2*DN[i].y*DN[j].y+k2*DN[i].x*DN[j].x;
			for (int m = 0; m < 3; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					Ke[3*i+m][3*j+n]=Ke[3*i+m][3*j+n]+D0*BiDBj[m][n]*weight*detJ;
				}
			}
		}
	}
	return;
}
void Tcell_8H::calculate_elastic_matrix(double &D0,double &k1,double &k2)
{
	double E=_gausspoint[0]->G_Youngs(),nu=_gausspoint[0]->G_Poisson();
	D0=E*(1-nu)/((1+nu)*(1-2*nu));
	k1=nu/(1-nu);
	k2=0.5*(1-2*nu)/(1-nu);
	return;
}
void Tcell_8H::calculate_DN_matrix(vec3D (&DN)[8],double &detJ,double xi,double eta,double zeta)
{
	double Jacobi[3][3];
	vec3D GN[8];
	bool degenerate;
	calcualte_Jacobi(Jacobi,GN,xi,eta,zeta,0);
	inverse(Jacobi,detJ,degenerate);
	for (int i = 0; i < 8; i++)
	{
		DN[i]=GN[i].multiply_by_matrix(Jacobi);
	}
	return;
}
void Tcell_8H::calcualte_Jacobi(double (&Jacobi)[3][3],vec3D (&GN)[8],double xi,double eta,double zeta,int type)
{
	for (int i = 0; i < 8; i++)
	{
		GN[i]=GN_brick_8(i,xi,eta,zeta);
	}
	vec3D temp[3];
	if (type==0)
	{
		//Calculate the Jacobi matrix with respected to the initial position
		for (int i=0;i<8;i++)
		{
			for (int j = 0; j < 3; j++)
			{
				temp[j]=temp[j]+GN[i]*_node_ptr[i]->_position0.access(j+1);
			}
		}
	}
	else if (type==1)
	{
		//Calculate the Jacobi matrix with respected to the current position
		for (int i=0;i<8;i++)
		{
			for (int j = 0; j < 3; j++)
			{
				temp[j]=temp[j]+GN[i]*_node_ptr[i]->_position.access(j+1);
			}
		}
	}
	else
	{
		cout<<"Error:Invalid type in calling calcualte_Jacobi() for 8H cell"<<endl;
		system("Pause");
		exit(0);
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Jacobi[i][j]=temp[j].access(i+1);
		}
	}
	return;
}
double Tcell_8H::calculate_distorted_quantity()
{
	double phy[8];
	for (int i = 0; i < 4; i++)
	{
		Ttetrahedron tet_temp1(_node_ptr[0]->_position0,_node_ptr[3]->_position0,_node_ptr[1]->_position0,_node_ptr[i+4]->_position0,1);
		Ttetrahedron tet_temp2(_node_ptr[1]->_position0,_node_ptr[2]->_position0,_node_ptr[3]->_position0,_node_ptr[i+4]->_position0,1);
		phy[2*i]=tet_temp1.G_distort_ratio();
		phy[2*i+1]=tet_temp2.G_distort_ratio();
	}
	double distort_ratio=0;
	for (int i = 0; i < 8; i++)
	{
		distort_ratio=maxval(distort_ratio,phy[i]);
	}
	return distort_ratio;
}
void Tcell_8H::modify_Youngs()
{
	double distorted_quantity=calculate_distorted_quantity();
	double Youngs=_gausspoint[0]->G_Youngs();
	//Youngs=1/sqrt(_distance_to_conode);
	Youngs=1/sqrt(_distance_to_conode);
	_gausspoint[0]->S_Youngs(Youngs);
	return;
}