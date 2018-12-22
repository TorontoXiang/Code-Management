#include "Tsubstrate.h"
Tsubstrate::Tsubstrate(double lx, double ly, double r_NW, double r_P, double l_NW, double density, double theta_max)
{
	double pi = 3.141592653;
	_lx = lx; _ly = ly; _r_NW = r_NW; _r_P = r_P;  _l_NW = l_NW; _density = density;
	_num_NW = lx*ly*density/(pi*r_NW*r_NW);
	//cout << "Create a substrate with straight NW" << endl;
	int num_input = 0;
	while (num_input < _num_NW)
	{
		vector<vec3D> new_NW = generate_valid_straight_NW(0, 1, theta_max);
		_NW_list.push_back(new_NW);
		_root_list.push_back(new_NW[0]);
		num_input = num_input + 1;
		//cout << num_input << endl;
	}
	//cout << "Finish the substrate generation" << endl;
}
Tsubstrate::Tsubstrate(double lx, double ly, double r_NW, double r_P, double l_NW, double density, double curvature, int n_divide)
{
	double pi = 3.141592653;
	_lx = lx; _ly = ly; _r_NW = r_NW; _r_P = r_P;  _n_divide = n_divide; _l_NW = l_NW; _density = density;
	_num_NW = lx * ly*density / (pi*r_NW*r_NW);
	//cout << "Create a substrate with wavy NW" << endl;
	int num_input = 0;
	vec3D x0, dx;
	double phy;
	bool is_survive;
	while (num_input < _num_NW)
	{
		vector<vec3D> new_NW = generate_valid_wavy_NW(0, 1, curvature,n_divide);
		_NW_list.push_back(new_NW);
		_root_list.push_back(new_NW[0]);
		num_input = num_input + 1;
		//cout << num_input << endl;
	}
	//cout << "Finish the substrate generation" << endl;
}
void Tsubstrate::generate_coupled_straight_substrate(double lx, double ly, double r_NW, double r_P, double l_NW, double density, double theta_max, Tsubstrate &coupled_substrate, bool &is_successed)
{
	double pi = 3.141592653;
	int num_trial=0;
	coupled_substrate._lx = lx; coupled_substrate._ly = ly;
	coupled_substrate._l_NW = l_NW; coupled_substrate._num_NW = lx * ly*density / (pi*r_NW*r_NW);
	coupled_substrate._r_NW = r_NW; coupled_substrate._r_P = r_P;
	coupled_substrate._density = density;
	int num_input = 0;
	vec3D x0, dx;
	double z = max_z();
	bool is_survive;
	is_successed = false;
	while (num_input<coupled_substrate._num_NW)
	{
		vector<vec3D> new_NW = coupled_substrate.generate_valid_straight_NW(z, -1, theta_max);
		num_trial = num_trial + 1;
		if (num_trial>=1e6)
		{
			return;
		}
		is_survive = true;
		for (int i = 0; i < _num_NW; i++)
		{
			if (is_survive)
			{
				if (!is_curve_valid(new_NW, _NW_list[i], 2 * _r_NW))
				{
					is_survive = false;
				}
			}
		}
		if (is_survive)
		{
			coupled_substrate._NW_list.push_back(new_NW);
			coupled_substrate._root_list.push_back(new_NW[0]);
			num_input = num_input + 1;
			//cout << num_input << " " << num_trial << endl;
			num_trial = 0;
		}
	}
	is_successed = true;
}
void Tsubstrate::generate_coupled_wavy_substrate(double lx, double ly, double r_NW, double r_P, double l_NW, double density, double curvature, int n_divide, Tsubstrate &coupled_substrate, bool &is_successed)
{
	double pi = 3.141592653;
	int num_trial = 0;
	coupled_substrate._lx = lx; coupled_substrate._ly = ly;
	coupled_substrate._l_NW = l_NW; coupled_substrate._num_NW = lx * ly*density / (pi*r_NW*r_NW);
	coupled_substrate._r_NW = r_NW; coupled_substrate._r_P = r_P;
	coupled_substrate._n_divide = _n_divide;
	int num_input = 0;
	vec3D x0, dx;
	double z = max_z();
	bool is_survive;
	is_successed = false;
	while (num_input < coupled_substrate._num_NW)
	{
		vector<vec3D> new_NW = coupled_substrate.generate_valid_wavy_NW(z, -1,curvature,coupled_substrate._n_divide);
		num_trial = num_trial + 1;
		if (num_trial >= 1e6)
		{
			return;
		}
		is_survive = true;
		for (int i = 0; i < _num_NW; i++)
		{
			if (is_survive)
			{
				if (!is_curve_valid(new_NW, _NW_list[i], 2 * _r_NW))
				{
					is_survive = false;
				}
			}
		}
		if (is_survive)
		{
			coupled_substrate._NW_list.push_back(new_NW);
			coupled_substrate._root_list.push_back(new_NW[0]);
			num_input = num_input + 1;
			//cout << num_input << " " << num_trial << endl;
			num_trial = 0;
		}
	}
	is_successed = true;
}
bool Tsubstrate::is_root_valid(vec3D x0)
{
	for (int i = 0; i < _root_list.size(); i++)
	{
		if ((x0-_root_list[i]).get_length()<(_r_NW+_r_P))
		{
			return false;
		}
	}
	return true;
}
void Tsubstrate::output_substrate(ofstream &output,int value)
{
	for (int i = 0; i < _NW_list.size(); i++)
	{
		TCNT_grid CNT_grid(_NW_list[i]);
		CNT_grid.generate_CNT_grid(_r_NW, 0.6*_r_NW, 1, 2);
		CNT_grid.output_CNT_tecplot(output,value);
	}
}
double Tsubstrate::max_z()
{
	double z_max = -100000;
	for (int i = 0; i < _num_NW; i++)
	{
		for (int j = 0; j < _NW_list[i].size(); j++)
		{
			z_max = maxval(z_max, _NW_list[i][j].z);
		}
	}
	return z_max;
}
vector<vec3D> Tsubstrate::generate_valid_straight_NW(double z, int direction, double theta_max)
{
	vec3D x0,dx;
	double phy;
	double pi = 3.141592653;
	bool is_survive;
	//srand((int)time(0));
	while (true)
	{
		x0.value(_lx*(1.0*rand() / RAND_MAX-0.5), _ly*(1.0*rand() / RAND_MAX-0.5), z);
		if (is_root_valid(x0))
		{
			phy = 2 * pi*rand() / RAND_MAX;
			if (direction==1)
			{
				dx.value(_l_NW*sin(theta_max)*cos(phy), _l_NW*sin(theta_max)*sin(phy), _l_NW*cos(theta_max));
			}
			else if (direction==-1)
			{
				dx.value(_l_NW*sin(theta_max)*cos(phy), _l_NW*sin(theta_max)*sin(phy), -_l_NW*cos(theta_max));
			}
			is_survive = true;
			vec3D p1_trial = x0;
			vec3D p2_trial = x0 + dx;
			for (int i = 0; i < _NW_list.size(); i++)
			{
				if (is_survive)
				{
					if (min_dis(p1_trial, p2_trial, _NW_list[i][0], _NW_list[i][1]) < 2 * _r_NW)
					{
						is_survive = false;
					}
				}
			}
			if (is_survive)
			{
				vector<vec3D> new_NW;
				new_NW.push_back(p1_trial); new_NW.push_back(p2_trial);
				return new_NW;
			}
		}
	}
}
vector<vec3D> Tsubstrate::generate_valid_wavy_NW(double z, int direction, double curvature, int n_divide)
{
	vec3D x0;
	bool is_survive;
	//srand((int)time(0));
	while (true)
	{
		x0.value(_lx*(1.0*rand() / RAND_MAX-0.5), _ly*(1.0*rand() / RAND_MAX-0.5), z);
		if (is_root_valid(x0))
		{
			vector<vec3D> new_NW = Generate_wavy_CNT(x0, curvature, n_divide, _l_NW, direction);
			is_survive = true;
			for (int i = 0; i < _NW_list.size(); i++)
			{
				if (is_survive)
				{
					if (!is_curve_valid(new_NW, _NW_list[i], 2 * _r_NW))
					{
						is_survive = false;
					}
				}
			}
			if (is_survive)
			{
				return new_NW;
			}
		}
	}
}