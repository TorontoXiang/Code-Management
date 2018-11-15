#include "TCNT_grid.h"
void Ttruncation::generate_truncation_grid(double r, double l, int m, int n)
{
	//Generate grid on two curve portion
	double pi = 3.141592654;
	int i_max1 = 4 * n + 4, j_max1 = m + 1;
	int i_max2 = n + 1, j_max2 = n + 1;
	int nump1 = i_max1 * (j_max1 + 1), nume1 = i_max1 * j_max1;
	int nump2 = (i_max2 - 1)*(j_max2 - 1);
	int nume2 = i_max2 * j_max2;
	_nume = nume1 + nume2; _nump = nump1 + nump2;
	_node_list = new vec3D[nump1 + nump2];
	_IEN_list = new int*[nume1 + nume2];
	for (int i = 0; i < nume1 + nume2; i++)
	{
		_IEN_list[i] = new int[4];
	}
	vec3D **pos;
	pos = new vec3D*[i_max1];
	for (int i = 0; i < i_max1; i++)
	{
		pos[i] = new vec3D[j_max1 + 1];
	}
	//Generate the outer grid
	for (int i = 0; i < i_max1; i++)
	{
		double theta = -3 * pi / 4 + 2 * pi * i / i_max1;
		pos[i][0].value(r*cos(theta), r*sin(theta), 0);
		if (i <= n + 1)
		{
			pos[i][j_max1].value(-0.5*l + l * i / (n + 1), -0.5*l, 0);
		}
		else if (i > n + 1 && i <= 2 * (n + 1))
		{
			pos[i][j_max1].value(0.5*l, -0.5*l + l * (i - n - 1) / (n + 1), 0);
		}
		else if (i > 2 * (n + 1) && i <= 3 * (n + 1))
		{
			pos[i][j_max1].value(0.5*l - l * (i - 2 * n - 2) / (n + 1), 0.5*l, 0);
		}
		else
		{
			pos[i][j_max1].value(-0.5*l, 0.5*l - l * (i - 3 * n - 3) / (n + 1), 0);
		}
	}
	for (int i = 0; i < i_max1; i++)
	{
		for (int j = 1; j < j_max1; j++)
		{
			vec3D temp = (pos[i][j_max1] - pos[i][0])*j / (m + 1);
			pos[i][j] = pos[i][0] + temp;
		}
	}
	for (int i = 0; i < i_max1; i++)
	{
		for (int j = 0; j < j_max1 + 1; j++)
		{
			int index = (i_max1)*j + i;
			_node_list[index] = pos[i][j];
		}
	}
	for (int i = 0; i < i_max1; i++)
	{
		for (int j = 0; j < j_max1; j++)
		{
			int index = i_max1 * j + i;
			int n1 = i_max1 * j + i, n2 = i_max1 * j + (i + 1) % i_max1, n3 = i_max1 * (j + 1) + (i + 1) % i_max1, n4 = i_max1 * (j + 1) + i;
			_IEN_list[index][0] = n1, _IEN_list[index][1] = n2, _IEN_list[index][2] = n3, _IEN_list[index][3] = n4;
		}
	}
	//Generate the inner grid
	vec3D o(-0.5*l, -0.5*l, 0);
	for (int i = 0; i < i_max2 - 1; i++)
	{
		for (int j = 0; j < j_max2 - 1; j++)
		{
			vec3D temp((i + 1)*l / i_max2, (j + 1)*l / j_max2, 0);
			int index = (i_max2 - 1)*j + i + nump1;
			_node_list[index] = o + temp;
		}
	}
	for (int i = 0; i < i_max2; i++)
	{
		for (int j = 0; j < j_max2; j++)
		{
			int index = i_max2 * j + i + nume1;
			int n1 = calculate_index(i, j, i_max2, i_max1, j_max1, nump1);
			int n2 = calculate_index(i + 1, j, i_max2, i_max1, j_max1, nump1);
			int n3 = calculate_index(i + 1, j + 1, i_max2, i_max1, j_max1, nump1);
			int n4 = calculate_index(i, j + 1, i_max2, i_max1, j_max1, nump1);
			_IEN_list[index][0] = n1, _IEN_list[index][1] = n2, _IEN_list[index][2] = n3, _IEN_list[index][3] = n4;
		}
	}
	for (int i = 0; i < i_max1; i++)
	{
		delete pos[i];
	}
	delete pos;
	for (int i = 0; i < _nump; i++)
	{
		_node_list[i] = _r0 + _node_list[i].multiply_by_matrix(_A);
	}
	return;
}
TCNT_grid::TCNT_grid(vector<vec3D> node_list)
{
	_num_node = node_list.size();
	_truncation_list.resize(_num_node);
	for (int i = 0; i < _num_node; i++)
	{
		_truncation_list[i]._r0 = node_list[i];
	}
	//Calculate the transformation matrix at each truncation face
	vec3D z0 = _truncation_list[1]._r0 - _truncation_list[0]._r0;
	z0.normalize();
	T_matrix(z0, _truncation_list[0]._A);
	_truncation_list[0]._c_t = 1;
	for (int i = 1; i < _num_node; i++)
	{
		//Calculate the normal of ith truncation in global system
		vec3D z;
		if (i< _num_node -1)
		{
			vec3D z1= _truncation_list[i]._r0 - _truncation_list[i-1]._r0;
			vec3D z2= _truncation_list[i+1]._r0 - _truncation_list[i]._r0;
			z1.normalize(); z2.normalize();
			z =(z1 + z2)*0.5;
			z.normalize();
			_truncation_list[i]._c_t = z * z1;
		}
		else
		{
			z = _truncation_list[i]._r0 - _truncation_list[i - 1]._r0;
			z.normalize();
			_truncation_list[i]._c_t = 1;
		}
		z=z.multiply_by_matrix_transpose(_truncation_list[i - 1]._A);  //The normal of ith truncation in the coordinate system of (i-1) truncation face
		double T[3][3],temp[3][3];
		T_matrix(z, T);         //Calculate the transformation from i to i-1 truncation face
		matrix_multiply(_truncation_list[i - 1]._A, T, temp);    //Calculate the transformation from i to global system
		//Calculate the axis of the new system and recalculate the transformation matrix to avoid numerical instability
		vec3D axis_x(temp[0][0],temp[1][0],temp[2][0]), axis_y(temp[0][1], temp[1][1], temp[2][1]), axis_z(temp[0][2], temp[1][2], temp[2][2]);
		axis_z.normalize();
		axis_y = (axis_z%axis_x); axis_y.normalize();
		axis_x = (axis_y%axis_z); axis_z.normalize();
		_truncation_list[i ]._A[0][0] = axis_x.x; _truncation_list[i]._A[1][0] = axis_x.y; _truncation_list[i]._A[2][0] = axis_x.z;
		_truncation_list[i]._A[0][1] = axis_y.x; _truncation_list[i]._A[1][1] = axis_y.y; _truncation_list[i]._A[2][1] = axis_y.z;
		_truncation_list[i]._A[0][2] = axis_z.x; _truncation_list[i]._A[1][2] = axis_z.y; _truncation_list[i]._A[2][2] = axis_z.z;
	}
	return;
}

void TCNT_grid::generate_CNT_grid(double r, double l, int m, int n)
{
	for (int i = 0; i < _num_node; i++)
	{
		//double ri=r/_truncation_list[i]._c_t, li=l/ _truncation_list[i]._c_t;
		_truncation_list[i].generate_truncation_grid(r, l, m, n);
	}
}
void TCNT_grid::output_CNT_tecplot(ofstream& output)
{
	int nump = _num_node * _truncation_list[0]._nump;
	int nume = (_num_node - 1)*_truncation_list[1]._nume;
	output << "TITLE = \"Tecplot Grid\"" << endl;
	output << "VARIABLES = \"X\",\"Y\",\"Z\"" << endl;
	output << "ZONE F=FEPOINT,N=" << nump << "," << "E=" << nume << "," << "ET=BRICK" << endl;
	int nump_each = _truncation_list[0]._nump;
	int nume_each = _truncation_list[0]._nume;
	for (int i = 0; i < _num_node; i++)
	{
		vec3D pos;
		for (int j = 0; j < nump_each; j++)
		{
			pos = _truncation_list[i]._node_list[j];
			output << pos.x << " " << pos.y << " " << pos.z << endl;
		}
	}
	for (int i = 0; i < _num_node - 1; i++)
	{
		for (int j = 0; j < nume_each; j++)
		{
			int IEN[8];
			for (int k = 0; k < 4; k++)
			{
				IEN[k] = _truncation_list[i]._IEN_list[j][k] + i * nump_each;
				IEN[k + 4] = IEN[k] + nump_each;
			}
			for (int k = 0; k < 8; k++)
			{
				output << IEN[k] + 1 << " ";
			}
			output << endl;
		}
	}
}
void TCNT_grid::output_CNT_k_file(ofstream& output,double E,double mu)
{
	output.precision(10);
	output << "*NEW_CNT" << endl;
	output << "*NODE" << endl;
	int nump_each = _truncation_list[0]._nump;
	int nume_each = _truncation_list[0]._nume;
	for (int i = 0; i < _num_node; i++)
	{
		vec3D pos;
		for (int j = 0; j < nump_each; j++)
		{
			pos = _truncation_list[i]._node_list[j];
			output << j+i*nump_each+1<<" "<<pos.x << " " << pos.y << " " << pos.z << endl;
		}
	}
	output << "*ELEMENT_SOLID" << endl;
	for (int i = 0; i < _num_node - 1; i++)
	{
		for (int j = 0; j < nume_each; j++)
		{
			output << j + i * nume_each + 1 << " " << 1 << " ";
			int IEN[8];
			for (int k = 0; k < 4; k++)
			{
				IEN[k] = _truncation_list[i]._IEN_list[j][k] + i * nump_each;
				IEN[k + 4] = IEN[k] + nump_each;
			}
			for (int k = 0; k < 8; k++)
			{
				output << IEN[k] + 1 << " ";
			}
			output << endl;
		}
	}
	output << "*MAT_ELASTIC" << endl;
	output << 1 << " 1 " << E << " " << mu << endl;
	output << "*END_CNT" << endl;
}
int calculate_index(int i, int j, int i_max2, int i_max1, int j_max1, int nump1)
{
	int index;
	int i1, j1 = j_max1;
	if (i == 0)
	{
		i1 = (4 * i_max2 - j) % (4 * i_max2); index = i_max1 * j1 + i1;
	}
	else if (j == 0)
	{
		i1 = i; index = i_max1 * j1 + i1;
	}
	else if (i == i_max2)
	{
		i1 = i_max2 + j; index = i_max1 * j1 + i1;
	}
	else if (j == i_max2)
	{
		i1 = 3 * i_max2 - i; index = i_max1 * j1 + i1;
	}
	else
	{
		i = i - 1; j = j - 1; index = (i_max2 - 1)*j + i + nump1;
	}
	return index;
}