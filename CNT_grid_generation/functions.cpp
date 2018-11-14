#include "functions.h"
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
		theta = pi * pi*rand() / RAND_MAX;
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
	cout << total_length << endl;
}