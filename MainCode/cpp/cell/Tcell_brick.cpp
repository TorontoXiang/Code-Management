#include "Tcell_brick.h"
#include "public_function.h"
#include <cmath>
#include "public_function.h"
#include "Gauss_weight.h"
#include "Shape_function.h"

Tcell_brick::Tcell_brick(int cell_id, int nGauss, Tnode* (&node_ptr)[8], TMAT_base** mat) : Tcell_base(cell_id, nGauss)
{
	//Allocate memory
	_corner_stress_force.resize(8);
	//Set cell node
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i] = node_ptr[i];
	}
	//Set Gauss point material
	for (int i = 0; i < nGauss; i++)
	{
		_gausspoint[i] = mat[i];
	}
	//Set geometry information
	_volume = calculate_cell_volume();
	_clength = calculate_characteristic_length();
	_soundspeed_max = _gausspoint[0]->G_soundspeed();
	//Initialize variables
	_strain_energy = 0;
	return;
}
double Tcell_brick::calculate_characteristic_length()
{
	int cell_face[6][4] = { {0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{0,4,7,3},{4,5,6,7} };
	vec3D p[4];
	double s_max = 0;
	for (int i = 0; i < 6; i++)
	{
		vec3D v1 = _node_ptr[cell_face[i][2]]->_position - _node_ptr[cell_face[i][0]]->_position;
		vec3D v2 = _node_ptr[cell_face[i][3]]->_position - _node_ptr[cell_face[i][1]]->_position;
		double s = (v1%v2).get_length()*0.5;
		s_max = maxval(s, s_max);
	}
	return _volume / s_max;
}
void Tcell_brick::calculate_Jacobi(double(&Jacobi)[3][3], vec3D(&GN)[8])
{
	for (int i = 0; i < 8; i++)
	{
		GN[i] = GN_brick_8(i, 0, 0, 0);
	}
	vec3D temp[3];
	//Calculate the Jacobi matrix with respected to the current position
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			temp[j] = temp[j] + GN[i] * _node_ptr[i]->_position.access(j + 1);
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
void Tcell_brick::assemble_corner_inertance()
{
	double mass = _volume *_gausspoint[0]->G_density();
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->accumulate_inertance(mass / 8.0);
	}
	return;
}
void Tcell_brick::calculate_DN_matrix(vec3D(&DN)[8], double &detJ)
{
	double Jacobi[3][3];
	vec3D GN[8];
	bool degenerate;
	calculate_Jacobi(Jacobi, GN);
	inverse(Jacobi, detJ, degenerate);
	for (int i = 0; i < 8; i++)
	{
		DN[i] = GN[i].multiply_by_matrix(Jacobi);
	}
	return;
}
void Tcell_brick::calculate_corner_force(Tviscosity_base* Arti_Vis, double dt)
{
	int num_Gauss = _gausspoint.size();
	vec3D DN[8];                                   //The derivatives of shape function with respected to x,y,z
	double det;
	calculate_DN_matrix(DN, det);
	double volume_pre = _volume;
	_volume = 8.0*det;
	//Calculate the gradient of the velocity
	double dij[3][3];
	for (int i = 0; i < 3; i++)
	{
		vec3D temp;
		for (int j = 0; j < 8; j++)
		{
			temp = temp + DN[j] * _node_ptr[j]->_velocity.access(i + 1);
		}
		dij[i][0] = temp.x; dij[i][1] = temp.y; dij[i][2] = temp.z;
	}
	//Calculate strain rate and rotation rate
	double de[6] = { 0 }, vort[3] = { 0 };         //The strain rate & rotation rate
	de[0] = dij[0][0]; de[1] = dij[1][1]; de[2] = dij[2][2];
	de[3] = dij[1][2] + dij[2][1]; de[4] = dij[0][2] + dij[2][0]; de[5] = dij[0][1] + dij[1][0];
	vort[0] = 0.5*(dij[1][2] - dij[2][1]);
	vort[1] = 0.5*(dij[0][2] - dij[2][0]);
	vort[2] = 0.5*(dij[0][1] - dij[1][0]);
	//Update the stress
	_gausspoint[0]->update_state(de, vort, dt, volume_pre);
	//Calcualte the artificial viscosity and modify the soundspeed
	double vrate = de[0] + de[1] + de[2];
	//double q = Arti_Vis->calculate_q(_gausspoint[0]->G_density(), vrate, _gausspoint[0]->G_soundspeed(), _clength);
	//double soundspeed = Arti_Vis->modify_soundspeed(vrate, _gausspoint[0]->G_soundspeed(), _clength);
	//_gausspoint[0]->S_soundspeed(soundspeed);
	//Calculate corner force
	double sig[6];
	_gausspoint[0]->G_stress(sig);
	double sxx = sig[0], syy = sig[1], szz = sig[2], sxy = sig[3], sxz = sig[4], syz = sig[5];
	for (int i = 0; i < 8; i++)
	{
		_corner_stress_force[i].x = - (sxx*DN[i].x + sxy * DN[i].y + sxz * DN[i].z)*_volume;
		_corner_stress_force[i].y = - (sxy*DN[i].x + syy * DN[i].y + syz * DN[i].z)*_volume;
		_corner_stress_force[i].z = - (sxz*DN[i].x + syz * DN[i].y + szz * DN[i].z)*_volume;
	}
	vec3D HG_force[8];
	calculate_hourglass_viscosity_force(HG_force);
	for (int i = 0; i < 8; i++)
	{
		_corner_stress_force[i] = _corner_stress_force[i] + HG_force[i];
	}
	return;
}
void Tcell_brick::calculate_hourglass_viscosity_force(vec3D(&HG_force)[8])
{
	double qh = 0.2;
	double a = qh * _gausspoint[0]->G_density() / 4.0;
	double v32 = pow(_volume, 2.0 / 3.0);
	double ah = a * v32 * _gausspoint[0]->G_soundspeed();
	double al = a * v32 * 100 * qh;
	double v2367[3], v1247[3], v0356[3], v0145[3];
	for (int i = 0; i < 3; i++)
	{
		v2367[i] = (_node_ptr[2]->_velocity - _node_ptr[3]->_velocity - _node_ptr[6]->_velocity + _node_ptr[7]->_velocity).access(i+1);
		v1247[i] = (_node_ptr[1]->_velocity - _node_ptr[2]->_velocity - _node_ptr[4]->_velocity + _node_ptr[7]->_velocity).access(i + 1);
		v0356[i] = (_node_ptr[0]->_velocity - _node_ptr[3]->_velocity - _node_ptr[5]->_velocity + _node_ptr[6]->_velocity).access(i + 1);
		v0145[i] = (_node_ptr[0]->_velocity - _node_ptr[1]->_velocity - _node_ptr[4]->_velocity + _node_ptr[5]->_velocity).access(i + 1);
	}
	double hgr[3][4];
	for (int i = 0; i < 3; i++)
	{
		hgr[i][0] = v0356[i] - v1247[i];
		hgr[i][1] = v0356[i] + v1247[i];
		hgr[i][2] = v0145[i] - v2367[i];
		hgr[i][3] = v0145[i] + v2367[i];
	}
	for (int i = 0; i < 3; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			hgr[i][k] = hgr[i][k] * (ah + abs(al*hgr[i][k]));
		}
	}
	double hap[3], ham[3], hbp[3], hbm[3];
	for (int i = 0; i < 3; i++)
	{
		hap[i] = hgr[i][0] + hgr[i][3];
		ham[i] = hgr[i][0] - hgr[i][3];
		hbp[i] = hgr[i][1] + hgr[i][2];
		hbm[i] = hgr[i][1] - hgr[i][2];
	}
	for (int i = 0; i < 3; i++)
	{
		HG_force[0].value(i + 1, -hap[i] - hbp[i]);
		HG_force[1].value(i + 1, hap[i] - hbm[i]);
		HG_force[2].value(i + 1, -hap[i] + hbp[i]);
		HG_force[3].value(i + 1, hap[i] + hbm[i]);
		HG_force[4].value(i + 1, -ham[i] + hbp[i]);
		HG_force[5].value(i + 1, ham[i] + hbm[i]);
		HG_force[6].value(i + 1, -ham[i] - hbp[i]);
		HG_force[7].value(i + 1, ham[i] - hbm[i]);
	}
	return;
}
void Tcell_brick::assemble_corner_force()
{
	for (int i = 0; i < 8; i++)
	{
		_node_ptr[i]->_force = _node_ptr[i]->_force + _corner_stress_force[i];
		_node_ptr[i]->_force0 = _node_ptr[i]->_force0 + _corner_stress_force[i];
	}
	return;
}
double Tcell_brick::calculate_time_step()
{
	_clength = calculate_characteristic_length();
	return _clength / _soundspeed_max;
}
void Tcell_brick::update_strain_energy(double dt)
{
	double de = 0;
	for (int i = 0; i < 4; i++)
	{
		de = de + _corner_stress_force[i] * _node_ptr[i]->_velocity;
		_corner_stress_force[i].value(0, 0, 0);
	}
	_strain_energy = _strain_energy + de * dt;
	return;
}
double Tcell_brick::calculate_cell_volume()
{
	double Jacobi[3][3];
	vec3D GN[8];
	double detJ;
	bool degenerate;
	calculate_Jacobi(Jacobi, GN);
	inverse(Jacobi, detJ, degenerate);
	return detJ * 8.0;
}
