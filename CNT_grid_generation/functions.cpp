#include "functions.h"
#include "TCNT_grid.h"
double min_dis(vec3D x1, vec3D x2, vec3D x3, vec3D x4)
{
	vec3D p = x1 - x3, n1 = x2 - x1, n3 = x4 - x3;
	double n1n1 = n1 * n1, n3n3 = n3 * n3, n1n3 = n1 * n3, n1p = n1 * p, n3p = n3 * p, pp = p * p;
	double delta = n1n1 * n3n3 - n1n3 * n1n3;
	double t, s;
	if (delta > 1e-15)
	{
		//n1 and n3 are not parallel
		t = (n1n3*n3p - n3n3 * n1p) / delta;
		s = (n3p + n1n3 * t) / n3n3;
		if (t >= 0 && t <= 1 && s >= 0 && s <= 1)
		{
			return sqrt(pp + 2 * n1p*t - 2 * n3p*s - 2 * n1n3*s*t + n1n1 * t*t + n3n3 * s*s);
		}
		else
		{
			double d1, d2, d3, d4, s;
			//minimal distance from x1 to x3x4
			s = (x1 - x3)*n3 / n3n3;
			if (s < 0) s = 0;
			else if (s > 1) s = 1;
			d1 = (x3 - x1 + (x4 - x3)*s).self_multuply();
			//minimal distance from x2 to x3x4
			s = (x2 - x3)*n3 / n3n3;
			if (s < 0) s = 0;
			else if (s > 1) s = 1;
			d2 = (x3 - x2 + (x4 - x3)*s).self_multuply();
			//minimal distance from x3 to x1x2
			s = (x3 - x1)*n1 / n1n1;
			if (s < 0) s = 0;
			else if (s > 1) s = 1;
			d3 = (x1 - x3 + (x2 - x1)*s).self_multuply();
			//minimal distance from x4 to x1x2
			s = (x4 - x1)*n1 / n1n1;
			if (s < 0) s = 0;
			else if (s > 1) s = 1;
			d4 = (x1 - x4 + (x2 - x1)*s).self_multuply();
			return sqrt(minval(minval(d1, d2), minval(d3, d4)));
		}
	}
	else
	{
		//n1 and n3 are parallel
		double v1 = -n3p, v2 = -n1n3 - n3p, v3 = n3n3 - n1n3 - n3p, v4 = n3n3 - n3p;
		if (v1*v2 <= 0) s = 0, t = -n3p / n1n3;
		else if (v2*v3 <= 0) t = 1, s = (n1n3 + n3p) / n3n3;
		else if (v3*v4 <= 0) s = 1, t = (n3n3 - n3p) / n1n3;
		else if (v4*v1 <= 0) t = 0, s = n3p / n3n3;
		else
		{
			double min_d = minval(minval(abs(v1), abs(v2)), minval(abs(v3), abs(v4)));
			if (min_d == abs(v1)) t = 0, s = 0;
			else if (min_d == abs(v2)) t = 1, s = 0;
			else if (min_d == abs(v3)) t = 1, s = 1;
			else t = 0, s = 1;
		}
	}
	return sqrt(pp + 2 * n1p*t - 2 * n3p*s - 2 * n1n3*s*t + n1n1 * t*t + n3n3 * s*s);;
}
void Create_straight_CNT_net(int num_CNT, double l_CNT,double lx, double ly, double lz)
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
	vec3D dx,dl;
	double phy,theta;
	double pi = 3.141592654;
	srand((int)time(0));
	bool is_survive;
	while (num_input<num_CNT)
	{
		dx.value(lx*rand() / RAND_MAX, ly*rand() / RAND_MAX, lz*rand() / RAND_MAX);
		phy = 2 * pi*rand() / RAND_MAX;
		theta = pi*rand() / RAND_MAX;
		dl.value(l_CNT*cos(theta)*cos(phy), l_CNT*cos(theta)*sin(phy), l_CNT*sin(theta));
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
						for (int l = 0; l < 27*num_input; l++)
						{
							if (min_dis(p1_trial,p2_trial,p1_list[l],p2_list[l])<=2*r+d_vdw)
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
	vector<vec3D> CNT_list1,CNT_list2;
	vec3D* ptr1=NULL;
	vec3D* ptr2=NULL;
	bool is_in = false;
	for (int i = 0; i < p1_list.size(); i++)
	{
		is_in = true;
		//x direction cutting
		x1 = p1_list[i], x2 = p2_list[i];
		double ratio;
		if (x1.x<0 && x2.x<0)
		{
			is_in = false;
		}
		else if (x1.x<0 && x2.x>=0)
		{
			ratio = -x2.x / (x1.x - x2.x);
			x1 = x2 + (x1 - x2)*ratio;
		}
		else if (x1.x>=0 && x2.x< 0)
		{
			ratio = -x2.x / (x1.x - x2.x);
			x2 = x2 + (x1 - x2)*ratio;
		}
		if (x1.x>lx && x2.x>lx)
		{
			is_in = false;
		}
		else if (x1.x<=lx && x2.x>lx)
		{
			ratio = (lx - x2.x) / (x1.x - x2.x);
			x2 = x2 + (x1 - x2)*ratio;
		}
		else if (x1.x>lx && x2.x<=lx)
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
	cout << "The total length is: "<<total_length << endl;

	//Output the CNT grid
	ofstream output_k,output_tec;
	output_k.open("CNT_grid.k");
	output_tec.open("CNT_grid.dat");
	double size_l = 0.67*3;
	for (int i = 0; i < CNT_list1.size(); i++)
	{
		double l = (CNT_list1[i] - CNT_list2[i]).get_length();
		int n_divided = int(l / size_l)+2;
		Generate_straight_CNT(CNT_list1[i], CNT_list2[i], n_divided, 0.335, 0.335*0.6, 1, 2, output_k,output_tec, 860, 0.17);
	}
	output_k << "*END" << endl;
}
void Create_wavy_CNT_net(int num_CNT, double l_CNT, int n_divided, double lx, double ly, double lz)
{
	double r = 0.335, d_vdw = 0.34, ratio = 3, E = 860, mu = 0.17;
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
							if (!is_curve_valid(trial_CNT,curve_list[l],2*r+d_vdw))
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
		for (int j = i+1; j < effective_curve_list.size(); j++)
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
	ofstream output_tec,output_k;
	output_tec.open("CNT_grid.dat");
	output_k.open("CNT_grid.k");
	for (int i = 0; i < effective_curve_list.size(); i=i++)
	{
		TCNT_grid CNT_grid(effective_curve_list[i]);
		CNT_grid.generate_CNT_grid(r, r*0.6, 1, 2);
		CNT_grid.output_CNT_tecplot(output_tec);
		if (effective_curve_list[i].size()!=2)
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
vector<vec3D> Generate_wavy_CNT(vec3D p_begin, double ratio, int n_divided, double l)
{
	double pi = 3.141592653;
	double curvature = ratio * 2 * pi / l;
	double dl = l / n_divided;
	double A[3][3], T[3][3],temp[3][3];
	A[0][0] = A[1][1] = A[2][2] = 1;
	A[0][1] = A[0][2] = A[1][0] = A[1][2] = A[2][0] = A[2][1] = 0;
	vector<vec3D> node_list;
	node_list.push_back(p_begin);
	double theta, phy;
	//srand((int)time(0));
	for (int i = 1; i < n_divided+1; i++)
	{
		if (i==1) theta = pi * rand() / RAND_MAX;          //First orientation is fully random
		else theta = curvature*dl;// *rand()*dl / RAND_MAX;   //The rotation of the rest segment is curvature

		phy = 2 * pi*rand() / RAND_MAX;
		vec3D z(sin(theta)*cos(phy), sin(theta)*sin(phy), cos(theta));
		vec3D dx=z*dl;
		dx = dx.multiply_by_matrix(A);
		node_list.push_back(node_list[i - 1] + dx);
		T_matrix(z, T);
		matrix_multiply(A, T, temp);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				A[i][j] = temp[i][j];
			}
		}
	}
	return node_list;
}
void T_matrix(vec3D z, double(&T)[3][3])
{
	double s_t, c_t, s_p, c_p;
	double eps = 1e-20;
	if (abs(z.z - 1) < 1e-15)
	{
		T[0][0] = 1; T[0][1] = T[0][2] = 0;
		T[1][1] = 1; T[1][0] = T[1][2] = 0;
		T[2][0] = T[2][1] = 0; T[2][2] = 1;
	}
	else
	{
		c_t = z.z; s_t = sqrt(1 - c_t * c_t);
		double l = sqrt(z.x*z.x + z.y*z.y);
		c_p = -z.y / l; s_p = -z.x / l;
		T[0][0] = c_p * c_p + s_p * s_p*c_t; T[0][1] = -c_p * s_p + s_p * c_p*c_t; T[0][2] = z.x;
		T[1][0] = -s_p * c_p + c_p * s_p*c_t; T[1][1] = s_p * s_p + c_p * c_p*c_t; T[1][2] = z.y;
		T[2][0] = s_t * s_p; T[2][1] = s_t * c_p; T[2][2] = z.z;
	}
	return;
}
void matrix_multiply(double(&A1)[3][3], double(&A2)[3][3], double(&A)[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			A[i][j] = 0;
			for (int k = 0; k < 3; k++)
			{
				A[i][j] = A[i][j] + A1[i][k] * A2[k][j];
			}
		}
	}
	return;
}
bool is_curve_valid(vector<vec3D> &curve1, vector<vec3D> &curve2, double dis_min)
{
	bool is_valid = true;
	vec3D x_min1, x_min2, x_max1, x_max2, x_min, x_max;
	curve_range(curve1, x_min1, x_max1); curve_range(curve2, x_min2, x_max2);
	//for (int i = 0; i < curve1.size()-1; i++)
	//{
	//	for (int j = 0; j < curve2.size()-1; j++)
	//	{
	//		if (min_dis(curve1[i],curve1[i+1],curve2[j],curve2[j+1])<dis_min)
	//		{
	//			is_valid=false;
	//		}
	//	}
	//}
	if (effective_range(x_min1,x_max1,x_min2,x_max2,x_min,x_max,2*dis_min))
	{
		for (int i = 0; i < curve1.size()-1; i++)
		{
			if (is_in_range(curve1[i],curve1[i+1],x_min,x_max))
			{
				for (int j = 0; j < curve2.size()-1; j++)
				{
					if (is_in_range(curve2[j], curve2[j + 1], x_min, x_max))
					{
						if (min_dis(curve1[i],curve1[i+1],curve2[j],curve2[j+1])<dis_min)
						{
							is_valid=false;
						}
					}
				}
			}
		}
	}
	return is_valid;
}
void curve_range(vector<vec3D> &curve, vec3D &x_min, vec3D &x_max)
{
	x_min.value(1e20, 1e20, 1e20), x_max.value(-1e20, -1e20, -1e20);
	for (int i = 0; i < curve.size(); i++)
	{
		x_min.x = minval(curve[i].x, x_min.x);
		x_min.y = minval(curve[i].y, x_min.y);
		x_min.z = minval(curve[i].z, x_min.z);
		x_max.x = maxval(curve[i].x, x_max.x);
		x_max.y = maxval(curve[i].y, x_max.y);
		x_max.z = maxval(curve[i].z, x_max.z);
	}
	return;
}
bool effective_range(vec3D &x_min1, vec3D &x_max1, vec3D& x_min2, vec3D x_max2, vec3D &x_min, vec3D &x_max, double min_dis)
{
	x_min.x = maxval(x_min1.x - min_dis, x_min2.x - min_dis);
	x_max.x = minval(x_max1.x + min_dis, x_max2.x + min_dis);
	x_min.y = maxval(x_min1.y - min_dis, x_min2.y - min_dis);
	x_max.y = minval(x_max1.y + min_dis, x_max2.y + min_dis);
	x_min.z = maxval(x_min1.z - min_dis, x_min2.z - min_dis);
	x_max.z = minval(x_max1.z + min_dis, x_max2.z + min_dis);
	if (x_min.x>x_max.x || x_min.y>x_max.y || x_min.z>x_max.z)
	{
		return false;
	}
	return true;
}
bool is_in_range(vec3D p1, vec3D p2,vec3D x_min, vec3D x_max)
{
	bool is_in = true;
	if (p1.x<x_min.x || p1.x>x_max.x || p1.y<x_min.y || p1.y>x_max.y || p1.z<x_min.z || p1.z>x_max.z)
	{
		is_in = false;
	}
	if (p2.x<x_min.x || p2.x>x_max.x || p2.y<x_min.y || p2.y>x_max.y || p2.z<x_min.z || p2.z>x_max.z)
	{
		is_in=false;
	}
	return is_in;
}
void move_curve(vector<vec3D> &curve, vec3D motion)
{
	for (int i = 0; i < curve.size(); i++)
	{
		curve[i] = curve[i] + motion;
	}
	return;
}
double min_curve_dis(vector<vec3D> &curve1, vector<vec3D> &curve2)
{
	double dis_min = 1e20;
	for (int i = 0; i < curve1.size()-1; i++)
	{
		for (int j = 0; j < curve2.size()-1; j++)
		{
			dis_min = minval(min_dis(curve1[i], curve1[i + 1], curve2[j], curve2[j + 1]), dis_min);
		}
	}
	return dis_min;
}
void effective_curve(vector<vec3D> &curve, vector<vector<vec3D>> &effective_curve_list, vec3D x_min, vec3D x_max)
{
	vector<int> index_start, index_end;
	bool is_in = false;
	for (int i = 0; i < curve.size()-1; i++)
	{
		if (is_in_range(curve[i],curve[i+1],x_min,x_max))
		{
			if (!is_in)
			{
				is_in = true;
				index_start.push_back(i);
			}
		}
		else
		{
			if (is_in)
			{
				is_in = false;
				index_end.push_back(i);
			}
		}
		if (i== curve.size() - 2 && is_in)
		{
			index_end.push_back(curve.size()-1);
		}
	}
	if (index_start.size() != index_end.size())
	{
		cout << "Error in effective_curve: the size of index_start and index_end must be equal!" << endl;
		system("Pause");
		exit(0);
	}
	int num_effective_CNT = index_start.size();
	effective_curve_list.resize(num_effective_CNT);
	for (int i = 0; i < num_effective_CNT; i++)
	{
		for (int j = index_start[i]; j <= index_end[i]; j++)
		{
			effective_curve_list[i].push_back(curve[j]);
		}
	}
	return;
}
