#include "functions.h"
#include <algorithm>
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
double min_dis(vec3D x1, vec3D x2, vec3D x3, vec3D x4, double &ratio1, double &ratio2)
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
			ratio1 = t; ratio2 = s;
			return sqrt(pp + 2 * n1p*t - 2 * n3p*s - 2 * n1n3*s*t + n1n1 * t*t + n3n3 * s*s);
		}
		else
		{
			vector<double> d;
			double t[4], s[4];
			d.resize(4);
			//minimal distance from x1 to x3x4
			t[0] = 0;
			s[0] = (x1 - x3)*n3 / n3n3;
			if (s[0] < 0) s[0] = 0;
			else if (s[0] > 1) s[0] = 1;
			d[0] = (x3 - x1 + (x4 - x3)*s[0]).self_multuply();
			//minimal distance from x2 to x3x4
			t[1] = 1;
			s[1] = (x2 - x3)*n3 / n3n3;
			if (s[1] < 0) s[1] = 0;
			else if (s[1] > 1) s[1] = 1;
			d[1] = (x3 - x2 + (x4 - x3)*s[1]).self_multuply();
			//minimal distance from x3 to x1x2
			s[2] = 0;
			t[2] = (x3 - x1)*n1 / n1n1;
			if (t[2] < 0) t[2] = 0;
			else if (t[2] > 1) t[2] = 1;
			d[2] = (x1 - x3 + (x2 - x1)*t[2]).self_multuply();
			//minimal distance from x4 to x1x2
			s[3] = 1;
			t[3] = (x4 - x1)*n1 / n1n1;
			if (t[3] < 0) t[3] = 0;
			else if (t[3] > 1) t[3] = 1;
			d[3] = (x1 - x4 + (x2 - x1)*t[3]).self_multuply();
			int min_pos = (int)(std::min_element(d.begin(), d.end()) - d.begin());
			ratio1 = t[min_pos]; ratio2 = s[min_pos];
			return sqrt(d[min_pos]);
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
	ratio1 = t; ratio2 = s;
	return sqrt(pp + 2 * n1p*t - 2 * n3p*s - 2 * n1n3*s*t + n1n1 * t*t + n3n3 * s*s);;
}
vector<vec3D> Generate_wavy_CNT(vec3D p_begin, double ratio, int n_divided, double l,int direction)
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
		if (i == 1)
		{
			if (direction==0)
			{
				//First orientation is fully random
				double c_theta = 1.0*rand() / RAND_MAX;  
				theta = acos(c_theta);
			}
			else if (direction==1)
			{
				//First orientation is along z direction
				theta = 0;
			}
			else if (direction==-1)
			{
				//First orientation is along -z direction
				theta = pi;
			}
			else
			{
				cout << "In valid direction in generating wavy CNT" << endl;
				system("Pause");
				exit(0);
			}
		}
		else
		{
			theta = curvature * dl;// *rand()*dl / RAND_MAX;   //The rotation of the rest segment is curvature
		}
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
	else if (abs(z.z + 1) < 1e-15)
	{
		T[0][0] = 1; T[0][1] = T[0][2] = 0;
		T[1][1] = -1; T[1][0] = T[1][2] = 0;
		T[2][0] = T[2][1] = 0; T[2][2] = -1;
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
int identify_segment_state(vec3D p1, vec3D p2, vec3D x_min, vec3D x_max)
{

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
double min_curve_dis(vector<vec3D> &curve1, vector<vec3D> &curve2, Selectrical_node &node1, Selectrical_node &node2)
{
	double dis_min = 1e20;
	double dis;
	double ratio1, ratio2;
	for (int i = 0; i < curve1.size() - 1; i++)
	{
		for (int j = 0; j < curve2.size() - 1; j++)
		{
			dis = min_dis(curve1[i], curve1[i + 1], curve2[j], curve2[j + 1], ratio1, ratio2);
			if (dis<dis_min)
			{
				dis_min = dis;
				node1.Segment_ID = i; node1.ratio = ratio1;
				node2.Segment_ID = j; node2.ratio = ratio2;
			}
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
void effective_curve(vector<vec3D> &curve, vector<vector<vec3D>> &effective_curve_list, vector<vector<int>> &bc_Info, vec3D x_min, vec3D x_max)
{
	vector<link_node> node_list;
	vector<link_node*> head_list;
	for (int i = 0; i < curve.size(); i++)
	{
		link_node new_node(curve[i]);
		node_list.push_back(new_node);
	}
	for (int i = 0; i < curve.size()-1; i++)
	{
		node_list[i].next_node = &node_list[i + 1];
	}
	link_node* new_head = &node_list[0];
	head_list.push_back(new_head);
}
void identify_straight_CNT_at_boundary(vec3D p1, vec3D p2, vector<int> &index, double lx, double ly, double lz)
{
	double eps = 1e-10;
	index.resize(12);
	for (int i = 0; i < 12; i++)
	{
		index[i] = 0;
	}
	if (abs(p1.x) < eps) index[0] = 1;
	if (abs(p1.x-lx) < eps) index[1] = 1;
	if (abs(p1.y) < eps) index[2] = 1;
	if (abs(p1.y - ly) < eps) index[3] = 1;
	if (abs(p1.z) < eps) index[4] = 1;
	if (abs(p1.z - lz) < eps) index[5] = 1;
	if (abs(p2.x) < eps) index[6] = 1;
	if (abs(p2.x - lx) < eps) index[7] = 1;
	if (abs(p2.y) < eps) index[8] = 1;
	if (abs(p2.y - ly) < eps) index[9] = 1;
	if (abs(p2.z) < eps) index[10] = 1;
	if (abs(p2.z - lz) < eps) index[11] = 1;
	return;
}
