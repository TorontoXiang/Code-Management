#include "Tcell.h"
#include "Shape_function.h"
#include "Gauss_weight.h"
Tcell::Tcell(int id, Snode* (&node_ptr)[8],int mid)
{
	_id = id;
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i] = node_ptr[i];
	}
	_mid = mid;
	_is_boundary = false;
	return;
}
void Tcell::calculate_element_stiffness(double(&Ke)[24][24],Smat& mat)
{
	//Initialize the element stiffness
	for (int i = 0; i < 24; i++)
	{
		for (int j = 0; j < 24; j++)
		{
			Ke[i][j] = 0;
		}
	}
	//Calculate the element stiffness
	double xi, eta, zeta, detJ, weight;
	vec3D DN[8];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				xi = Gauss_coordinate[1][i], eta = Gauss_coordinate[1][j], zeta = Gauss_coordinate[1][k];
				weight = Gauss_weight[1][i] * Gauss_weight[1][j] * Gauss_weight[1][k];
				calculate_DN_matrix(DN, detJ, xi, eta, zeta);
				accumulate_Gausspoint_stiffness(Ke, DN, detJ, weight,mat);
			}
		}
	}
	return;
}
void Tcell::calculate_DN_matrix(vec3D(&DN)[8], double &detJ, double xi, double eta, double zeta)
{
	double Jacobi[3][3];
	vec3D GN[8];
	bool degenerate;
	calcualte_Jacobi(Jacobi, GN, xi, eta, zeta);
	inverse(Jacobi, detJ, degenerate);
	for (int i = 0; i < 8; i++)
	{
		DN[i] = GN[i].multiply_by_matrix(Jacobi);
	}
	return;
}
void Tcell::calcualte_Jacobi(double(&Jacobi)[3][3], vec3D(&GN)[8], double xi, double eta, double zeta)
{
	for (int i = 0; i < 8; i++)
	{
		GN[i] = GN_brick_8(i, xi, eta, zeta);
	}
	vec3D temp[3];
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			temp[j] = temp[j] + GN[i] * _node_ptr[i]->_pos[j];
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Jacobi[i][j] = temp[j].access(i + 1);
		}
	}
	return;
}
void Tcell::accumulate_Gausspoint_stiffness(double(&Ke)[24][24], vec3D(&DN)[8], double detJ, double weight,Smat& mat)
{
	double D0, k1, k2;
	double E = mat._E, mu = mat._mu;
	D0 = E * (1 - mu) / ((1 + mu)*(1 - 2 * mu));
	k1 = mu / (1 - mu);
	k2 = 0.5*(1 - 2 * mu) / (1 - mu);
	double BiDBj[3][3];
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			BiDBj[0][0] = DN[i].x*DN[j].x + k2 * DN[i].y*DN[j].y + k2 * DN[i].z*DN[j].z;
			BiDBj[1][0] = k1 * DN[i].y*DN[j].x + k2 * DN[i].x*DN[j].y;
			BiDBj[2][0] = k1 * DN[i].z*DN[j].x + k2 * DN[i].x*DN[j].z;
			BiDBj[0][1] = k1 * DN[i].x*DN[j].y + k2 * DN[i].y*DN[j].x;
			BiDBj[1][1] = DN[i].y*DN[j].y + k2 * DN[i].x*DN[j].x + k2 * DN[i].z*DN[j].z;
			BiDBj[2][1] = k1 * DN[i].z*DN[j].y + k2 * DN[i].y*DN[j].z;
			BiDBj[0][2] = k1 * DN[i].x*DN[j].z + k2 * DN[i].z*DN[j].x;
			BiDBj[1][2] = k1 * DN[i].y*DN[j].z + k2 * DN[i].z*DN[j].y;
			BiDBj[2][2] = DN[i].z*DN[j].z + k2 * DN[i].y*DN[j].y + k2 * DN[i].x*DN[j].x;
			for (int m = 0; m < 3; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					Ke[3 * i + m][3 * j + n] = Ke[3 * i + m][3 * j + n] + D0 * BiDBj[m][n] * weight*detJ;
				}
			}
		}
	}
	return;
}
Sstress Tcell::calculate_gausspiont_stress(vec3D(&DN)[8], double(&d_cell)[8][3],Smat& mat)
{
	double D0;
	double E = mat._E, mu = mat._mu;
	D0 = E / ((1 + mu)*(1 - 2 * mu));
	double exx, eyy, ezz, eyz, exz, exy;
	exx = eyy = ezz = eyz = exz = exy = 0;
	for (int i = 0; i < 8; i++)
	{
		exx = exx + d_cell[i][0] * DN[i].x;
		eyy = eyy + d_cell[i][1] * DN[i].y;
		ezz = ezz + d_cell[i][2] * DN[i].z;
		exy = exy + d_cell[i][0] * DN[i].y + d_cell[i][1] * DN[i].x;
		exz = exz + d_cell[i][0] * DN[i].z + d_cell[i][2] * DN[i].x;
		eyz = eyz + d_cell[i][1] * DN[i].z + d_cell[i][2] * DN[i].y;
	}
	Sstress stress;
	stress.sxx = D0 * ((1 - mu)*exx + mu*eyy + mu * ezz);
	stress.syy = D0 * (mu*exx + (1 - mu)*eyy + mu * ezz);
	stress.szz = D0 * (mu*exx + mu * eyy + (1 - mu)*ezz);
	stress.syz = D0 * (1 - 2 * mu)*0.5*eyz;
	stress.sxz = D0 * (1 - 2 * mu)*0.5*exz;
	stress.sxy = D0 * (1 - 2 * mu)*0.5*exy;
	return stress;
}
void Tcell::calculate_cell_stress(Sstress* cell_stress, double(&d_cell)[8][3], Smat& mat)
{
	double xi, eta, zeta, detJ;
	vec3D DN[8];
	int index;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				index = local_node_id(i,j,k);
				xi = Gauss_coordinate[1][i], eta = Gauss_coordinate[1][j], zeta = Gauss_coordinate[1][k];
				calculate_DN_matrix(DN, detJ, xi, eta, zeta);
				cell_stress[index] = calculate_gausspiont_stress(DN, d_cell, mat);
			}
		}
	}
}