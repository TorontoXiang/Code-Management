#include "CNT_net_generation.h"
#include "readin.h"
#include <iomanip>
#include <istream>
#include <string>
using namespace std;
Slink_node::Slink_node(vec3D coor, int State)
{
	pos = coor;
	state = State;
	pre_node = NULL;
	next_node = NULL;
}
void Create_straight_CNT_net(int num_CNT, double l_CNT, double lx, double ly, double lz)
{
	double r = 0.335, d_vdw = 0.34;
	vec3D x_min[3][3][3];
	//Create the minimal coordinate of the 27 boxes
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				x_min[i][j][k].value((i - 1)*lx, (j - 1)*ly, (k - 1)*lz);
			}
		}
	}
	int num_input = 0;
	vector<vec3D> p1_list, p2_list;
	vec3D dx, dl;
	double phy, theta;
	double pi = 3.141592654;
	srand((int)time(0));
	bool is_survive;
	while (num_input < num_CNT)
	{
		dx.value(lx*rand() / RAND_MAX, ly*rand() / RAND_MAX, lz*rand() / RAND_MAX);
		phy = 2 * pi*rand() / RAND_MAX;
		double c_theta = 1.0*rand() / RAND_MAX;
		theta = acos(c_theta);
		dl.value(l_CNT*sin(theta)*cos(phy), l_CNT*sin(theta)*sin(phy), l_CNT*c_theta);
		is_survive = true;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					if (is_survive)
					{
						vec3D p1_trial = x_min[i][j][k] + dx;
						vec3D p2_trial = x_min[i][j][k] + dx + dl;
						for (int l = 0; l < 27 * num_input; l++)
						{
							if (min_dis(p1_trial, p2_trial, p1_list[l], p2_list[l]) <= 2 * r + d_vdw)
							{
								is_survive = false;
							}
						}
					}
				}
			}
		}
		if (is_survive)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						p1_list.push_back(x_min[i][j][k] + dx);
						p2_list.push_back(x_min[i][j][k] + dx + dl);
					}
				}
			}
			num_input = num_input + 1;
			cout << num_input << endl;
		}
	}
	vec3D x1, x2;
	vector<vec3D> CNT_list1, CNT_list2;
	vec3D* ptr1 = NULL;
	vec3D* ptr2 = NULL;
	bool is_in = false;
	for (int i = 0; i < p1_list.size(); i++)
	{
		is_in = true;
		//x direction cutting
		x1 = p1_list[i], x2 = p2_list[i];
		double ratio;
		if (x1.x < 0 && x2.x < 0)
		{
			is_in = false;
		}
		else if (x1.x < 0 && x2.x >= 0)
		{
			ratio = -x2.x / (x1.x - x2.x);
			x1 = x2 + (x1 - x2)*ratio;
		}
		else if (x1.x >= 0 && x2.x < 0)
		{
			ratio = -x2.x / (x1.x - x2.x);
			x2 = x2 + (x1 - x2)*ratio;
		}
		if (x1.x > lx && x2.x > lx)
		{
			is_in = false;
		}
		else if (x1.x <= lx && x2.x > lx)
		{
			ratio = (lx - x2.x) / (x1.x - x2.x);
			x2 = x2 + (x1 - x2)*ratio;
		}
		else if (x1.x > lx && x2.x <= lx)
		{
			ratio = (lx - x2.x) / (x1.x - x2.x);
			x1 = x2 + (x1 - x2)*ratio;
		}
		//y direction cutting
		if (x1.y < 0 && x2.y < 0)
		{
			is_in = false;
		}
		else if (x1.y < 0 && x2.y >= 0)
		{
			ratio = -x2.y / (x1.y - x2.y);
			x1 = x2 + (x1 - x2)*ratio;
		}
		else if (x1.y >= 0 && x2.y < 0)
		{
			ratio = -x2.y / (x1.y - x2.y);
			x2 = x2 + (x1 - x2)*ratio;
		}
		if (x1.y > ly && x2.y > ly)
		{
			is_in = false;
		}
		else if (x1.y <= lx && x2.y > lx)
		{
			ratio = (ly - x2.y) / (x1.y - x2.y);
			x2 = x2 + (x1 - x2)*ratio;
		}
		else if (x1.y > lx && x2.y <= ly)
		{
			ratio = (ly - x2.y) / (x1.y - x2.y);
			x1 = x2 + (x1 - x2)*ratio;
		}
		//z direction cutting
		if (x1.z < 0 && x2.z < 0)
		{
			is_in = false;
		}
		else if (x1.z < 0 && x2.z >= 0)
		{
			ratio = -x2.z / (x1.z - x2.z);
			x1 = x2 + (x1 - x2)*ratio;
		}
		else if (x1.z >= 0 && x2.z < 0)
		{
			ratio = -x2.z / (x1.z - x2.z);
			x2 = x2 + (x1 - x2)*ratio;
		}
		if (x1.z > lz && x2.z > lz)
		{
			is_in = false;
		}
		else if (x1.z <= lz && x2.z > lz)
		{
			ratio = (lz - x2.z) / (x1.z - x2.z);
			x2 = x2 + (x1 - x2)*ratio;
		}
		else if (x1.z > lz && x2.z <= lz)
		{
			ratio = (lz - x2.z) / (x1.z - x2.z);
			x1 = x2 + (x1 - x2)*ratio;
		}
		if (is_in)
		{
			CNT_list1.push_back(x1), CNT_list2.push_back(x2);
		}
	}
	double total_length = 0;
	for (int i = 0; i < CNT_list1.size(); i++)
	{
		total_length = total_length + (CNT_list1[i] - CNT_list2[i]).get_length();
	}
	cout << "The total length is: " << total_length << endl;

	//Output the CNT grid
	int m = 1, n = 2;
	ofstream output_k, output_tec,output_bc;
	output_k.open("CNT_grid.k");
	output_tec.open("CNT_grid.dat");
	output_bc.open("CNT_bc.dat");
	output_bc << CNT_list1.size() << " " << m << " " << n << endl;
	double size_l = 0.67 * 3;
	vector<int> index;
	for (int i = 0; i < CNT_list1.size(); i++)
	{
		double l = (CNT_list1[i] - CNT_list2[i]).get_length();
		int n_divided = int(l / size_l) + 2;
		Generate_straight_CNT(CNT_list1[i], CNT_list2[i], n_divided, 0.335, 0.335*0.6, m, n, output_k, output_tec, 860, 0.17);
		identify_straight_CNT_at_boundary(CNT_list1[i], CNT_list2[i], index, lx, ly, lz);
		for (int i = 0; i < 12; i++)
		{
			output_bc << index[i] << " ";
		}
		output_bc << endl;
	}
	output_k << "*END" << endl;
}
void TCNT_net_generator::Create_straight_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info,int output_control)
{
	vec3D x_min[3][3][3];
	//Create the minimal coordinate of the 27 boxes
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				x_min[i][j][k].value((i - 1)*_lx, (j - 1)*_ly, (k - 1)*_lz);
			}
		}
	}
	double volume_CNT = 0.25*_d_CNT*_d_CNT*pi*_l_CNT;
	int num_target = int(_lx * _ly*_lz*_volume_fraction / volume_CNT);
	int num_input = 0;
	vector<vec3D> p1_list, p2_list;
	vec3D dx, dl;
	double phy, theta;
	double pi = 3.141592654;
	bool is_survive;
	srand((int)time(0));
	while (num_input < num_target)
	{
		dx.value(_lx*rand() / RAND_MAX, _ly*rand() / RAND_MAX, _lz*rand() / RAND_MAX);
		phy = 2 * pi*rand() / RAND_MAX;
		double c_theta = 1.0*rand() / RAND_MAX;
		theta = acos(c_theta);
		dl.value(_l_CNT*sin(theta)*cos(phy), _l_CNT*sin(theta)*sin(phy), _l_CNT*c_theta);
		is_survive = true;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					if (is_survive)
					{
						vec3D p1_trial = x_min[i][j][k] + dx;
						vec3D p2_trial = x_min[i][j][k] + dx + dl;
						for (int l = 0; l < 27 * num_input; l++)
						{
							if (min_dis(p1_trial, p2_trial, p1_list[l], p2_list[l]) <= _d_CNT + _d_VDW)
							{
								is_survive = false;
							}
						}
					}
				}
			}
		}
		if (is_survive)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						p1_list.push_back(x_min[i][j][k] + dx);
						p2_list.push_back(x_min[i][j][k] + dx + dl);
					}
				}
			}
			num_input = num_input + 1;
			//cout << num_input << endl;
			if (output_control==0)
			{
				if (num_input==1)
				{
					cout << setw(4) << 100 * num_input / num_target << "%";
				}
				else
				{
					cout << '\b' << '\b' << '\b' << '\b' << '\b';
					cout << setw(4) << 100 * num_input / num_target << "%";
				}
			}
		}
	}
	if (output_control==0)
	{
		cout << endl;
	}
	vec3D x1, x2;
	vector<vec3D> CNT_list1, CNT_list2;
	vec3D* ptr1 = NULL;
	vec3D* ptr2 = NULL;
	bool is_in = false;
	for (int i = 0; i < p1_list.size(); i++)
	{
		x1 = p1_list[i], x2 = p2_list[i];
		Clipping_straight_CNT(x1, x2, is_in);
		if (is_in)
		{
			CNT_list1.push_back(x1), CNT_list2.push_back(x2);
		}
	}
	double size_l = _l_CNT/_n_divided;
	for (int i = 0; i < CNT_list1.size(); i++)
	{
		double l = (CNT_list1[i] - CNT_list2[i]).get_length();
		int N_segment = int(l / size_l) + 2;
		vector<vec3D> node_list;
		vector<int> index;
		for (int j = 0; j < N_segment+1; j++)
		{
			vec3D pos = CNT_list1[i] + (CNT_list2[i] - CNT_list1[i])*j / N_segment;
			node_list.push_back(pos);
		}
		identify_straight_CNT_at_boundary(CNT_list1[i], CNT_list2[i], index, _lx, _ly, _lz);
		CNT_net.push_back(node_list);
		bc_Info.push_back(index);
	}
	return;
}
void TCNT_net_generator::Create_straight_CNT_net_advanced(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info, int output_control)
{
	vec3D x_min[3][3][3];
	//Create the minimal coordinate of the 27 boxes
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				x_min[i][j][k].value((i - 1)*_lx, (j - 1)*_ly, (k - 1)*_lz);
			}
		}
	}
	double volume_CNT = 0.25*_d_CNT*_d_CNT*pi*_l_CNT;
	int num_target = int(_lx * _ly*_lz*_volume_fraction / volume_CNT);
	int num_input = 0;
	vector<vec3D> p1_list, p2_list;
	vec3D dx, dl;
	double phy, theta;
	double pi = 3.141592654;
	bool is_survive;
	srand((int)time(0));
	while (num_input < num_target)
	{
		dx.value(_lx*rand() / RAND_MAX, _ly*rand() / RAND_MAX, _lz*rand() / RAND_MAX);
		phy = 2 * pi*rand() / RAND_MAX;
		double c_theta = 1.0*rand() / RAND_MAX;
		theta = acos(c_theta);
		dl.value(_l_CNT*sin(theta)*cos(phy), _l_CNT*sin(theta)*sin(phy), _l_CNT*c_theta);
		is_survive = true;
		vector<vec3D> p1_trial_list;
		vector<vec3D> p2_trial_list;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					if (is_survive)
					{
						vec3D p1_trial = x_min[i][j][k] + dx;
						vec3D p2_trial = x_min[i][j][k] + dx + dl;
						bool is_in=true;
						Clipping_straight_CNT(p1_trial, p2_trial, is_in);
						if (is_in)
						{
							for (int l = 0; l < p1_list.size(); l++)
							{
								if (min_dis(p1_trial, p2_trial, p1_list[l], p2_list[l]) <= _d_CNT + _d_VDW)
								{
									is_survive = false;
								}
							}
							if (is_survive)
							{
								p1_trial_list.push_back(p1_trial);
								p2_trial_list.push_back(p2_trial);
							}
						}
					}
				}
			}
		}
		if (is_survive)
		{
			for (int i = 0; i < p1_trial_list.size(); i++)
			{
				p1_list.push_back(p1_trial_list[i]);
				p2_list.push_back(p2_trial_list[i]);
			}
			num_input = num_input + 1;
			if (output_control == 0)
			{
				if (num_input == 1)
				{
					cout << setw(4) << 100 * num_input / num_target << "%";
				}
				else
				{
					cout << '\b' << '\b' << '\b' << '\b' << '\b';
					cout << setw(4) << 100 * num_input / num_target << "%";
				}
			}
		}
	}
	if (output_control == 0)
	{
		cout << endl;
	}
	//vec3D x1, x2;
	//vector<vec3D> CNT_list1, CNT_list2;
	//vec3D* ptr1 = NULL;
	//vec3D* ptr2 = NULL;
	//bool is_in = false;
	//for (int i = 0; i < p1_list.size(); i++)
	//{
	//	x1 = p1_list[i], x2 = p2_list[i];
	//	Clipping_straight_CNT(x1, x2, is_in);
	//	if (is_in)
	//	{
	//		CNT_list1.push_back(x1), CNT_list2.push_back(x2);
	//	}
	//}
	double size_l = _l_CNT / _n_divided;
	for (int i = 0; i < p1_list.size(); i++)
	{
		double l = (p1_list[i] - p2_list[i]).get_length();
		int N_segment = int(l / size_l) + 2;
		vector<vec3D> node_list;
		vector<int> index;
		for (int j = 0; j < N_segment + 1; j++)
		{
			vec3D pos = p1_list[i] + (p2_list[i] - p1_list[i])*j / N_segment;
			node_list.push_back(pos);
		}
		identify_straight_CNT_at_boundary(p1_list[i], p2_list[i], index, _lx, _ly, _lz);
		CNT_net.push_back(node_list);
		bc_Info.push_back(index);
	}
	return;
}
void Create_wavy_CNT_net(int num_CNT, double l_CNT, int n_divided, double ratio, double lx, double ly, double lz)
{
	double r = 0.335, d_vdw = 0.34, E = 860, mu = 0.17;
	vec3D x_min[3][3][3];
	//Create the minimal coordinate of the 27 boxes
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				x_min[i][j][k].value((i - 1)*lx, (j - 1)*ly, (k - 1)*lz);
			}
		}
	}
	int num_input = 0;
	vector<vector<vec3D>> curve_list;
	vec3D p_begin;
	srand((int)time(0));
	bool is_survive;
	double extension = l_CNT / n_divided;
	while (num_input < num_CNT)
	{
		p_begin.value(lx*rand() / RAND_MAX, ly*rand() / RAND_MAX, lz*rand() / RAND_MAX);
		vector<vec3D> new_curve = Generate_wavy_CNT(p_begin, ratio, n_divided, l_CNT);
		is_survive = true;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					if (is_survive)
					{
						vector<vec3D> trial_CNT = new_curve;
						move_curve(trial_CNT, x_min[i][j][k]);
						for (int l = 0; l < 27 * num_input; l++)
						{
							if (!is_curve_valid(trial_CNT, curve_list[l], 2 * r + d_vdw))
							{
								is_survive = false;
							}
						}
					}
				}
			}
		}
		if (is_survive)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						vector<vec3D> trial_CNT = new_curve;
						move_curve(trial_CNT, x_min[i][j][k]);
						curve_list.push_back(trial_CNT);
					}
				}
			}
			num_input = num_input + 1;
			cout << num_input << endl;
		}
	}
	cout << "Calculate effective curve" << endl;
	vec3D min(0, 0, 0), max(lx, ly, lz);
	vector<vector<vec3D>> effective_curve_list;
	for (int i = 0; i < curve_list.size(); i++)
	{
		vector<vector<vec3D>> local_effective_curve;
		effective_curve(curve_list[i], local_effective_curve, min, max);
		for (int i = 0; i < local_effective_curve.size(); i++)
		{
			effective_curve_list.push_back(local_effective_curve[i]);
		}
	}
	cout << "Distance check" << endl;
	cout << "The number of effective curve is " << effective_curve_list.size() << endl;
	for (int i = 0; i < effective_curve_list.size(); i++)
	{
		for (int j = i + 1; j < effective_curve_list.size(); j++)
		{
			double dis_mini = min_curve_dis(effective_curve_list[i], effective_curve_list[j]);
			if (dis_mini < 2 * r + d_vdw)
			{
				cout << "Minimal distance conflict!" << endl;
				cout << i << " " << j << " " << dis_mini << endl;
			}
		}
	}
	cout << "Output CNT net" << endl;
	ofstream output_tec, output_k;
	output_tec.open("CNT_grid.dat");
	output_k.open("CNT_grid.k");
	for (int i = 0; i < effective_curve_list.size(); i = i++)
	{
		TCNT_grid CNT_grid(effective_curve_list[i]);
		CNT_grid.generate_CNT_grid(r, r*0.6, 1, 2);
		CNT_grid.output_CNT_tecplot(output_tec);
		if (effective_curve_list[i].size() != 2)
		{
			CNT_grid.output_CNT_k_file(output_k, E, mu);
		}
	}
	output_k << "*END" << endl;
}
void TCNT_net_generator::Create_wavy_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info,int output_control)
{
	vec3D x_min[3][3][3];
	//Create the minimal coordinate of the 27 boxes
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				x_min[i][j][k].value((i - 1)*_lx, (j - 1)*_ly, (k - 1)*_lz);
			}
		}
	}
	int num_input = 0;
	double volume_CNT = 0.25*_d_CNT*_d_CNT*pi*_l_CNT;
	int num_target = int(_lx * _ly*_lz*_volume_fraction / volume_CNT);
	vector<vector<vec3D>> curve_list;
	vec3D p_begin;
	srand((int)time(0));
	bool is_survive;
	while (num_input < num_target)
	{
		p_begin.value(_lx*rand() / RAND_MAX, _ly*rand() / RAND_MAX, _lz*rand() / RAND_MAX);
		vector<vec3D> new_curve = Generate_wavy_CNT(p_begin, _curvature_ratio, _n_divided, _l_CNT);
		is_survive = true;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					if (is_survive)
					{
						vector<vec3D> trial_CNT = new_curve;
						move_curve(trial_CNT, x_min[i][j][k]);
						for (int l = 0; l < 27 * num_input; l++)
						{
							if (!is_curve_valid(trial_CNT, curve_list[l], _d_CNT + _d_VDW))
							{
								is_survive = false;
							}
						}
					}
				}
			}
		}
		if (is_survive)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						vector<vec3D> trial_CNT = new_curve;
						move_curve(trial_CNT, x_min[i][j][k]);
						curve_list.push_back(trial_CNT);
					}
				}
			}
			num_input = num_input + 1;
			//cout << num_input << endl;
			if (output_control==0)
			{
				if (num_input == 1)
				{
					cout << setw(4) << 100 * num_input / num_target << "%";
				}
				else
				{
					cout << '\b' << '\b' << '\b' << '\b' << '\b';
					cout << setw(4) << 100 * num_input / num_target << "%";
				}
			}
		}
	}
	if (output_control==0)
	{
		cout << endl;
		cout << "Calculate effective curve" << endl;
	}
	for (int i = 0; i < curve_list.size(); i++)
	{
		vector<vector<vec3D>> sub_CNT;
		vector<vector<int>> sub_bc_Info;
		Clipping_curved_CNT(curve_list[i], sub_CNT, sub_bc_Info);
		for (int j = 0; j < sub_CNT.size(); j++)
		{
			CNT_net.push_back(sub_CNT[j]);
			bc_Info.push_back(sub_bc_Info[j]);
		}
	}
	for (int i = 0; i < CNT_net.size(); i++)
	{
		if (CNT_net[i].size() == 2)
		{
			CNT_net[i].resize(3);
			vec3D pos = (CNT_net[i][0] + CNT_net[i][1])*0.5;
			CNT_net[i][2] = CNT_net[i][1];
			CNT_net[i][1] = pos;
		}
	}
}
void TCNT_net_generator::Create_wavy_CNT_net_advanced(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info, int output_control)
{
	vec3D x_min[3][3][3];
	//Create the minimal coordinate of the 27 boxes
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				x_min[i][j][k].value((i - 1)*_lx, (j - 1)*_ly, (k - 1)*_lz);
			}
		}
	}
	int num_input = 0;
	double volume_CNT = 0.25*_d_CNT*_d_CNT*pi*_l_CNT;
	int num_target = int(_lx * _ly*_lz*_volume_fraction / volume_CNT);
	//vector<vector<vec3D>> curve_list;
	vec3D p_begin;
	srand((int)time(0));
	bool is_survive;
	while (num_input < num_target)
	{
		p_begin.value(_lx*rand() / RAND_MAX, _ly*rand() / RAND_MAX, _lz*rand() / RAND_MAX);
		vector<vec3D> new_curve = Generate_wavy_CNT(p_begin, _curvature_ratio, _n_divided, _l_CNT);
		is_survive = true;
		vector<vector<vec3D>> sub_CNT_list;
		vector<vector<int>> sub_bc_info_list;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					if (is_survive)
					{
						vector<vec3D> trial_CNT = new_curve;
						move_curve(trial_CNT, x_min[i][j][k]);
						vector<vector<vec3D>> sub_CNT;
						vector<vector<int>> sub_bc_Info;
						Clipping_curved_CNT(trial_CNT, sub_CNT, sub_bc_Info);
						for (int l = 0; l < sub_CNT.size(); l++)
						{
							for (int  m = 0; m < CNT_net.size(); i++)
							{
								if (!is_curve_valid(sub_CNT[l], CNT_net[m], _d_CNT + _d_VDW))
								{
									is_survive = false;
								}
							}
						}
						if (is_survive)
						{
							for (int l = 0; l < sub_CNT.size(); l++)
							{
								sub_CNT_list.push_back(sub_CNT[l]);
								sub_bc_info_list.push_back(sub_bc_Info[l]);
							}
						}
					}
				}
			}
		}
		if (is_survive)
		{
			for (int i = 0; i < sub_CNT_list.size(); i++)
			{
				CNT_net.push_back(sub_CNT_list[i]);
				bc_Info.push_back(sub_bc_info_list[i]);
			}
			num_input = num_input + 1;
			if (output_control == 0)
			{
				if (num_input == 1)
				{
					cout << setw(4) << 100 * num_input / num_target << "%";
				}
				else
				{
					cout << '\b' << '\b' << '\b' << '\b' << '\b';
					cout << setw(4) << 100 * num_input / num_target << "%";
				}
			}
		}
	}
	if (output_control == 0)
	{
		cout << endl;
		cout << "Calculate effective curve" << endl;
	}
	//for (int i = 0; i < curve_list.size(); i++)
	//{
	//	vector<vector<vec3D>> sub_CNT;
	//	vector<vector<int>> sub_bc_Info;
	//	Clipping_curved_CNT(curve_list[i], sub_CNT, sub_bc_Info);
	//	for (int j = 0; j < sub_CNT.size(); j++)
	//	{
	//		CNT_net.push_back(sub_CNT[j]);
	//		bc_Info.push_back(sub_bc_Info[j]);
	//	}
	//}
	for (int i = 0; i < CNT_net.size(); i++)
	{
		if (CNT_net[i].size() == 2)
		{
			CNT_net[i].resize(3);
			vec3D pos = (CNT_net[i][0] + CNT_net[i][1])*0.5;
			CNT_net[i][2] = CNT_net[i][1];
			CNT_net[i][1] = pos;
		}
	}
	return;
}
void Generate_straight_CNT(vec3D p1, vec3D p2, int n_divided, double r, double l, int m, int n, ofstream& output_k, ofstream& output_tec, double E, double mu)
{
	vector<vec3D> node_list;
	for (int i = 0; i < n_divided + 1; i++)
	{
		vec3D pos = p1 + (p2 - p1)*i / n_divided;
		node_list.push_back(pos);
	}
	TCNT_grid CNT_grid(node_list);
	CNT_grid.generate_CNT_grid(r, l, m, n);
	CNT_grid.output_CNT_k_file(output_k, E, mu);
	CNT_grid.output_CNT_tecplot(output_tec);
	return;
}
void TCNT_net_generator::Generate_CNT_net_grid(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info, ofstream &output_k, ofstream &output_tec, ofstream &output_bc)
{
	output_bc << CNT_net.size() << " " << _m << " " << _n << endl;

	for (int i = 0; i < CNT_net.size(); i = i++)
	{
		TCNT_grid CNT_grid(CNT_net[i]);
		CNT_grid.generate_CNT_grid(0.5*_d_CNT, 0.3*_d_CNT, _m, _n);
		CNT_grid.output_CNT_tecplot(output_tec);
		CNT_grid.output_CNT_k_file(output_k, _E, _mu);
		for (int j = 0; j < 12; j++)
		{
			output_bc << bc_Info[i][j] << " ";
		}
		output_bc << endl;
	}
	output_k << "*END" << endl;
}
void TCNT_net_generator::input_RVE_info(ifstream &input)
{
	string a;
	while (true)
	{
		next_keyword(input);
		input >> a;
		if (a == "*RVE_SIZE")
		{
			next_data(input);
			input >> _lx >> _ly >> _lz;
		}
		else if (a == "*CNT_SIZE")
		{
			next_data(input);
			input >> _d_CNT >> _l_CNT;
		}
		else if (a == "*CNT_MATERIAL_PROPERTY")
		{
			next_data(input);
			input >> _E >> _mu;
		}
		else if (a == "*CNT_GRID_SIZE")
		{
			next_data(input);
			input >> _m >> _n >> _n_divided;
		}
		else if (a == "*CNT_TYPE")
		{
			next_data(input);
			input >> _type >> _curvature_ratio;
		}
		else if (a == "*CNT_VOLUME_FRACTION")
		{
			next_data(input);
			input >> _volume_fraction;
		}
		else if (a == "*VDW_DISTANCE")
		{
			next_data(input);
			input >> _d_VDW;
		}
		else if (a == "*END")
		{
			return;
		}
	}
}
void TCNT_net_generator::Create_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> &bc_Info,int output_control)
{
	if (_type==0)
	{
		Create_straight_CNT_net_advanced(CNT_net,bc_Info,output_control);
	}
	else if (_type==1)
	{
		Create_wavy_CNT_net(CNT_net, bc_Info,output_control);
	}
	else
	{
		cout << "Error: Invalid CNT network type!!" << endl;
	}
	return;
}
void TCNT_net_generator::Clipping_curved_CNT(vector<vec3D> &curved_CNT, vector<vector<vec3D>> &sub_CNT, vector<vector<int>> &bc_info)
{
	int num_segment = curved_CNT.size() - 1;
	vector<int> segment_state_list(num_segment);
	for (int i = 0; i < num_segment; i++)
	{
		segment_state_list[i] = segment_state(curved_CNT[i], curved_CNT[i + 1]);
	}
	vector<int> segment_start, segment_end;
	bool is_in = false;
	for (int i = 0; i < num_segment; i++)
	{
		if (segment_state_list[i]==-1)
		{
			if (!is_in)
			{
				is_in = true;
				segment_start.push_back(i);
			}
		}
		else
		{
			if (is_in)
			{
				is_in = false;
				segment_end.push_back(i-1);
			}
		}
		if (i==num_segment-1 && is_in)
		{
			segment_end.push_back(num_segment-1);
		}
	}
	if (segment_start.size() != segment_end.size())
	{
		cout << "Error in effective_curve: the size of index_start and index_end must be equal!" << endl;
		system("Pause");
		exit(0);
	}
	int num_sub_CNT = segment_start.size();
	sub_CNT.resize(num_sub_CNT);
	for (int i = 0; i < num_sub_CNT; i++)
	{
		vector<int> index(12);
		for (int j = segment_start[i]; j <= segment_end[i]+1; j++)
		{
			sub_CNT[i].push_back(curved_CNT[j]);
		}
		int state_start, state_end;
		if (segment_start[i]==0)
		{
			state_start = -1;
		}
		else
		{
			state_start = segment_state_list[segment_start[i]-1];
			index[state_start] = 1;
		}
		if (segment_end[i]==num_segment-1)
		{
			state_end = -1;
		}
		else
		{
			state_end = segment_state_list[segment_end[i] + 1];
			index[state_end + 6] = 1;
		}
		bc_info.push_back(index);
	}
	//vector<int> index_start, index_end;
	//bool is_in = false;
	//for (int i = 0; i < curve.size() - 1; i++)
	//{
	//	if (is_in_range(curve[i], curve[i + 1], x_min, x_max))
	//	{
	//		if (!is_in)
	//		{
	//			is_in = true;
	//			index_start.push_back(i);
	//		}
	//	}
	//	else
	//	{
	//		if (is_in)
	//		{
	//			is_in = false;
	//			index_end.push_back(i);
	//		}
	//	}
	//	if (i == curve.size() - 2 && is_in)
	//	{
	//		index_end.push_back(curve.size() - 1);
	//	}
	//}
	//if (index_start.size() != index_end.size())
	//{
	//	cout << "Error in effective_curve: the size of index_start and index_end must be equal!" << endl;
	//	system("Pause");
	//	exit(0);
	//}
	//int num_effective_CNT = index_start.size();
	//effective_curve_list.resize(num_effective_CNT);
	//for (int i = 0; i < num_effective_CNT; i++)
	//{
	//	for (int j = index_start[i]; j <= index_end[i]; j++)
	//	{
	//		effective_curve_list[i].push_back(curve[j]);
	//	}
	//}
	return;
}
int TCNT_net_generator::segment_state(vec3D p1, vec3D p2)
{
	bool is_p1_inside = is_inside_RVE(p1);
	bool is_p2_inside = is_inside_RVE(p2);
	if (is_p1_inside && is_p2_inside)
	{
		return -1;          //Segment is inside the RVE
	}
	else if (!is_p1_inside && !is_p2_inside)
	{
		return -2;          //Segment is outside the RVE
	}
	else if (is_p1_inside)
	{
		if (p2.x<=0) return 0;
		else if (p2.y <= 0) return 1;
		else if (p2.z <= 0) return 2;
		else if (p2.x >= _lx) return 3;
		else if (p2.y >= _ly) return 4;
		else if (p2.z >= _lz) return 5;
	}
	else if (is_p2_inside)
	{
		if (p1.x <= 0) return 0;
		else if (p1.y <= 0) return 1;
		else if (p1.z <= 0) return 2;
		else if (p1.x >= _lx) return 3;
		else if (p1.y >= _ly) return 4;
		else if (p1.z >= _lz) return 5;
	}
	cout << "Error in segment_state" << endl;
	system("Pause");
	exit(0);
	return 0;
}
bool TCNT_net_generator::is_inside_RVE(vec3D p)
{
	if (p.x>0 && p.y>0 && p.z>0 && p.x<_lx && p.y<_ly && p.z<_lz)
	{
		return true;
	}
	else
	{
		return false;
	}
}
void TCNT_net_generator::Clipping_straight_CNT(vec3D &x1, vec3D &x2, bool &is_in)
{
	is_in = true;
	//x direction cutting
	double ratio;
	if (x1.x < 0 && x2.x < 0)
	{
		is_in = false;
		return;
	}
	else if (x1.x < 0 && x2.x >= 0)
	{
		ratio = -x2.x / (x1.x - x2.x);
		x1 = x2 + (x1 - x2)*ratio;
	}
	else if (x1.x >= 0 && x2.x < 0)
	{
		ratio = -x2.x / (x1.x - x2.x);
		x2 = x2 + (x1 - x2)*ratio;
	}
	if (x1.x > _lx && x2.x > _lx)
	{
		is_in = false;
		return;
	}
	else if (x1.x <= _lx && x2.x > _lx)
	{
		ratio = (_lx - x2.x) / (x1.x - x2.x);
		x2 = x2 + (x1 - x2)*ratio;
	}
	else if (x1.x > _lx && x2.x <= _lx)
	{
		ratio = (_lx - x2.x) / (x1.x - x2.x);
		x1 = x2 + (x1 - x2)*ratio;
	}
	//y direction cutting
	if (x1.y < 0 && x2.y < 0)
	{
		is_in = false;
		return;
	}
	else if (x1.y < 0 && x2.y >= 0)
	{
		ratio = -x2.y / (x1.y - x2.y);
		x1 = x2 + (x1 - x2)*ratio;
	}
	else if (x1.y >= 0 && x2.y < 0)
	{
		ratio = -x2.y / (x1.y - x2.y);
		x2 = x2 + (x1 - x2)*ratio;
	}
	if (x1.y > _ly && x2.y > _ly)
	{
		is_in = false;
		return;
	}
	else if (x1.y <= _ly && x2.y > _ly)
	{
		ratio = (_ly - x2.y) / (x1.y - x2.y);
		x2 = x2 + (x1 - x2)*ratio;
	}
	else if (x1.y > _ly && x2.y <= _ly)
	{
		ratio = (_ly - x2.y) / (x1.y - x2.y);
		x1 = x2 + (x1 - x2)*ratio;
	}
	//z direction cutting
	if (x1.z < 0 && x2.z < 0)
	{
		is_in = false;
		return;
	}
	else if (x1.z < 0 && x2.z >= 0)
	{
		ratio = -x2.z / (x1.z - x2.z);
		x1 = x2 + (x1 - x2)*ratio;
	}
	else if (x1.z >= 0 && x2.z < 0)
	{
		ratio = -x2.z / (x1.z - x2.z);
		x2 = x2 + (x1 - x2)*ratio;
	}
	if (x1.z > _lz && x2.z > _lz)
	{
		is_in = false;
		return;
	}
	else if (x1.z <= _lz && x2.z > _lz)
	{
		ratio = (_lz - x2.z) / (x1.z - x2.z);
		x2 = x2 + (x1 - x2)*ratio;
	}
	else if (x1.z > _lz && x2.z <= _lz)
	{
		ratio = (_lz - x2.z) / (x1.z - x2.z);
		x1 = x2 + (x1 - x2)*ratio;
	}
	return;
}