#include "TCNT_net.h"
#include <algorithm>
void new_line(ifstream& input)
{
	char a;
	do
	{
		input.get(a);
	} while (a != '\n');
}
vector<vec3D> SCNT::calculate_centroid_line()
{
	vector<vec3D> centroid_line;
	int cell_in_truncation = (n + 1)*(m + 1) * 4 + (n + 1)*(n + 1);
	int num_node = nume / cell_in_truncation + 1;
	centroid_line.resize(num_node);
	double area = 0, cell_area = 0;
	vec3D moment;
	moment.x = moment.y = moment.z = 0;
	int cell_id;
	vec3D p1, p2, p3, p4, cell_moment;
	for (int i = 0; i < num_node; i++)
	{
		for (int j = 0; j < cell_in_truncation; j++)
		{
			if (i==num_node-1)
			{
				cell_id = (i-1) * cell_in_truncation + j;
				p1 = node_list[IEN_list[cell_id][4]-1]; p2 = node_list[IEN_list[cell_id][5]-1];
				p3 = node_list[IEN_list[cell_id][6]-1]; p4 = node_list[IEN_list[cell_id][7]-1];
			}
			else
			{
				cell_id = i * cell_in_truncation + j;
				p1 = node_list[IEN_list[cell_id][0]-1]; p2 = node_list[IEN_list[cell_id][1]-1];
				p3 = node_list[IEN_list[cell_id][2]-1]; p4 = node_list[IEN_list[cell_id][3]-1];
			}
			calculate_moment(p1, p2, p3, p4, cell_moment, cell_area);
			moment = moment + cell_moment;
			area = area + cell_area;
		}
		centroid_line[i] = moment / area;
		moment.value(0, 0, 0); area = 0;
	}
	return centroid_line;
}
void TCNT_net::input_CNT_net()
{
	ifstream input_CNT_info, input_bc_info;
	input_CNT_info.open("CNT_grid.dat");
	input_bc_info.open("CNT_bc.dat");
	int m, n;
	input_bc_info >> _num_CNT >> m >> n;
	_bc_info.resize(_num_CNT);
	_CNT_net.resize(_num_CNT);
	for (int i = 0; i < _num_CNT; i++)
	{
		_bc_info[i].resize(12);
		for (int j = 0; j < 12; j++)
		{
			input_bc_info >> _bc_info[i][j];
		}
	}
	string a, b, c;
	int num_input = 0;
	while (num_input<_num_CNT)
	{
		num_input = num_input + 1;
		while (true)
		{
			input_CNT_info >> a;
			if (a=="ZONE")
			{
				SCNT CNT(m,n);
				input_CNT_info >> b >> CNT.nump >> c >> CNT.nume;
				CNT.node_list.resize(CNT.nump);
				CNT.IEN_list.resize(CNT.nume);
				new_line(input_CNT_info);
				for (int i = 0; i < CNT.nump; i++)
				{
					input_CNT_info >> CNT.node_list[i].x >> CNT.node_list[i].y >> CNT.node_list[i].z;
					new_line(input_CNT_info);
				}
				for (int i = 0; i < CNT.nume; i++)
				{
					CNT.IEN_list[i].resize(8);
					for (int j = 0; j < 8; j++)
					{
						input_CNT_info >> CNT.IEN_list[i][j];
					}
				}
				_CNT_net[num_input - 1] = CNT.calculate_centroid_line();
				break;
			}
		}
	}
	cout << "Finish" << endl;
	return;
}
void calculate_moment(vec3D &p1, vec3D &p2, vec3D &p3, vec3D &p4, vec3D &moment, double &area)
{
	double s1 = triagnle_area(p1, p2, p3);
	double s2 = triagnle_area(p1, p3, p4);
	vec3D c1 = (p1 + p2 + p3) / 3.0;
	vec3D c2 = (p1 + p3 + p4) / 3.0;
	moment = c1 * s1 + c2 * s2;
	area = s1 + s2;
	return;
}
double triagnle_area(vec3D &p1, vec3D &p2, vec3D &p3)
{
	double nx = (p2.y - p1.y)*(p3.z - p1.z) - (p2.z - p1.z)*(p3.y - p1.y);
	double ny = (p2.z - p1.z)*(p3.x - p1.x) - (p2.x - p1.x)*(p3.z - p1.z);
	double nz = (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
	double s = sqrt(nx*nx + ny * ny + nz * nz);
	return 0.5*s;
}
void TCNT_net::Create_tunneling_resistance()
{
	_num_tunneling_R = 0;
	for (int i = 0; i < _num_CNT-1; i++)
	{
		for (int j = i+1; j < _num_CNT; j++)
		{
			Selectrical_node node1, node2;
			double dis_min = min_curve_dis(_CNT_net[i], _CNT_net[j], node1, node2);
			if (dis_min <_d_tunneling)
			{
				_num_tunneling_R = _num_tunneling_R + 1;
				node1.CNT_ID = i; node2.CNT_ID = j;
				node1.Electric_node_ID = 2 * _num_tunneling_R - 2;
				node2.Electric_node_ID = 2 * _num_tunneling_R - 1;
				double R = 100;
				Sresistance_T new_tunneling_resistance(node1, node2, R);
				Sresistance new_resistance(2 * _num_tunneling_R - 2, 2 * _num_tunneling_R - 1, R);
				_RT_list.push_back(new_tunneling_resistance);
				_resistance_net.push_back(new_resistance);
			}
		}
	}
	return;
}
void TCNT_net::Create_intrinsic_resistance()
{
	for (int i = 0; i < _num_CNT; i++)
	{
		Intrinsic_resistance_on_CNT(i);
	}
	return;
}
void TCNT_net::Intrinsic_resistance_on_CNT(int CNT_ID)
{
	Selectrical_node node_begin, node_end;
	node_begin.CNT_ID = CNT_ID; node_begin.Electric_node_ID = 2 * _num_tunneling_R + 2 * CNT_ID; 
	node_begin.Segment_ID = 0; node_begin.ratio = 0;
	node_end.CNT_ID = CNT_ID; node_end.Electric_node_ID = 2 * _num_tunneling_R + 2 * CNT_ID+1; 
	node_end.Segment_ID = _CNT_net[CNT_ID].size() - 2; node_end.ratio = 1;
	vector<Selectrical_node> Electric_node_list;
	Electric_node_list.push_back(node_begin);
	Electric_node_list.push_back(node_end);
	for (int i = 0; i < _num_tunneling_R; i++)
	{
		if (_RT_list[i].node1.CNT_ID==CNT_ID)
		{
			Electric_node_list.push_back(_RT_list[i].node1);
		}
		if (_RT_list[i].node2.CNT_ID==CNT_ID)
		{
			Electric_node_list.push_back(_RT_list[i].node2);
		}
	}
	vector<double> CNT_segment_length;
	for (int i = 0; i < _CNT_net[CNT_ID].size()-1; i++)
	{
		CNT_segment_length.push_back((_CNT_net[CNT_ID][i + 1] - _CNT_net[CNT_ID][i]).get_length());
	}
	sort(Electric_node_list.begin(), Electric_node_list.end());
	for (int i = 0; i < Electric_node_list.size()-1; i++)
	{
		double l_total;
		if (Electric_node_list[i].Segment_ID== Electric_node_list[i+1].Segment_ID)
		{
			l_total = (Electric_node_list[i + 1].ratio - Electric_node_list[i].ratio)*CNT_segment_length[Electric_node_list[i + 1].Segment_ID];
		}
		else
		{
			double l_begin = (1 - Electric_node_list[i].ratio)*CNT_segment_length[Electric_node_list[i].Segment_ID];
			double l_end= Electric_node_list[i+1].ratio*CNT_segment_length[Electric_node_list[i+1].Segment_ID];
			l_total = l_begin+l_end;
			for (int j = Electric_node_list[i].Segment_ID+1; j <= Electric_node_list[i + 1].Segment_ID-1; j++)
			{
				l_total = l_total + CNT_segment_length[j];
			}
		}
		Sresistance new_resistance(Electric_node_list[i].Electric_node_ID, Electric_node_list[i + 1].Electric_node_ID, l_total);
		_resistance_net.push_back(new_resistance);
	}
	return;
}