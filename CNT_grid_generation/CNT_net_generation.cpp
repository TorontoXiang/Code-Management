#include "CNT_net_generation.h"
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
		//theta = 2*pi*rand() / RAND_MAX;
		dl.value(l_CNT*sin(theta)*cos(phy), l_CNT*sin(theta)*sin(phy), l_CNT*c_theta);
		//dl.value(l_CNT*cos(theta)*cos(phy), l_CNT*cos(theta)*sin(phy), l_CNT*sin(theta));
		//dx.value(1, 1, 1);
		//theta = pi*0.5;
		//phy = pi*0.5;
		//dl.value(l_CNT*cos(theta)*cos(phy), l_CNT*cos(theta)*sin(phy), l_CNT*sin(theta));
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
	ofstream output_k, output_tec;
	output_k.open("CNT_grid.k");
	output_tec.open("CNT_grid.dat");
	double size_l = 0.67 * 3;
	for (int i = 0; i < CNT_list1.size(); i++)
	{
		double l = (CNT_list1[i] - CNT_list2[i]).get_length();
		int n_divided = int(l / size_l) + 2;
		Generate_straight_CNT(CNT_list1[i], CNT_list2[i], n_divided, 0.335, 0.335*0.6, 1, 2, output_k, output_tec, 860, 0.17);
	}
	output_k << "*END" << endl;
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