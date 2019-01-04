#include "TCNT_net.h"
#include <algorithm>
#include <cmath>
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
	cout << "Input the CNT net..." << endl;
	ifstream input_CNT_info, input_bc_info;
	input_CNT_info.open("CNT_grid.dat");
	input_bc_info.open("bc_Info.k");
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
	return;
}
void TCNT_net::input_CNT_net(vector<vector<vec3D>> &CNT_net, vector<vector<int>> bc_Info)
{
	_num_CNT = CNT_net.size();
	_CNT_net.resize(_num_CNT);
	_bc_info.resize(_num_CNT);
	for (int i = 0; i < _num_CNT; i++)
	{
		_CNT_net[i] = CNT_net[i];
		_bc_info[i] = bc_Info[i];
	}
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
int identify_RVE_face(string face_name)
{
	if (face_name == "x_min") return 0;
	else if (face_name == "y_min") return 1;
	else if (face_name == "z_min") return 2;
	else if (face_name == "x_max") return 3;
	else if (face_name == "y_max") return 4;
	else if (face_name == "z_max") return 5;
	else
	{
		cout << "Invalid face name!" << endl;
		system("Pause");
		exit(0);
	}
}
void TCNT_net::Create_tunneling_resistance()
{
	//Constants:
	double pi = 3.141592653;
	double d_cnt = 0.67e-9;      //The diameter of the CNT
	double d_cutoff = 1.4e-9;    //The maximum tunneling thickness (m)
	double Me = 9.10938e-31;     //The mass of a electron (kg)
	double h = 6.62607004e-34;   //The Planck constant (J*s)
	double hh = h/(2*pi);        //Dirac constant h/2pi (J*s)
	double d_Vdw = 0.34e-9;      //van der wall separaion distance (nm)

	double dtunnel = hh / sqrt(8 * Me*_dE);  //Bao's paper
	double Tmax = exp(-d_Vdw / dtunnel); //in case of contact with no insulation between them(d = D + dvdw)
	double R_tunnel_min = 12905.4 / (Tmax*_M); // where 12905.4 is h / 2e ^ 2
	_num_tunneling_R = 0;
	for (int i = 0; i < _num_CNT-1; i++)
	{
		for (int j = i+1; j < _num_CNT; j++)
		{
			Selectrical_node node1, node2;
			double dis_min = min_curve_dis(_CNT_net[i], _CNT_net[j], node1, node2);
			dis_min = dis_min * 1e-9;       //Change the unit to meter
			if (dis_min<=d_cnt+d_cutoff)
			{
				_num_tunneling_R = _num_tunneling_R + 1;
				node1.CNT_ID = i; node2.CNT_ID = j;
				node1.Electric_node_ID = 2 * _num_tunneling_R - 2 + 2 * _num_CNT;
				node2.Electric_node_ID = 2 * _num_tunneling_R - 1 + 2 * _num_CNT;
				double R_tunnel;
				if (dis_min <d_cnt+d_Vdw)
				{
					R_tunnel = R_tunnel_min;

				}
				else
				{
					double T = exp(-(dis_min - d_cnt) / dtunnel);
					R_tunnel = 12905.4 / (T*_M);
				}
				Sresistance_T new_tunneling_resistance(node1, node2, R_tunnel);
				Sresistance new_resistance(node1.Electric_node_ID, node2.Electric_node_ID, R_tunnel);
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
	double pi = 3.141592653;
	double d_cnt = 0.67e-9;      //The diameter of the CNT
	Selectrical_node node_begin, node_end;
	node_begin.CNT_ID = CNT_ID; node_begin.Electric_node_ID = 2 * CNT_ID; 
	node_begin.Segment_ID = 0; node_begin.ratio = 0;
	node_end.CNT_ID = CNT_ID; node_end.Electric_node_ID = 2 * CNT_ID + 1;
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
		double R_intrinsic = 4 * l_total*1e-9 / (pi*_sigma_CNT*d_cnt*d_cnt);
		if (R_intrinsic<1e-20)
		{
			R_intrinsic = 1e-5;
		}
		Sresistance new_resistance(Electric_node_list[i].Electric_node_ID, Electric_node_list[i + 1].Electric_node_ID, R_intrinsic);
		_resistance_net.push_back(new_resistance);
	}
	return;
}
double TCNT_net::Calculate_effective_resistance(string positive_face, string negative_face, int output_control)
{
	if (output_control==0)
	{
		cout << "Create the effective resistance net" << endl;
	}
	Create_tunneling_resistance();
	Create_intrinsic_resistance();
	int positive_face_id = identify_RVE_face(positive_face);
	int negative_face_id = identify_RVE_face(negative_face);
	int num_positive_node=0, num_negative_node=0;
	int num_electrical_node = 2 * _num_CNT + 2 * _num_tunneling_R;
	vector<int> Identify(num_electrical_node);
	vector<int> ID(num_electrical_node);
	for (int i = 0; i < num_electrical_node; i++)
	{
		Identify[i] = 0;
	}
	for (int i = 0; i < _num_CNT; i++)
	{
		if (_bc_info[i][positive_face_id]==1)
		{
			num_positive_node = num_positive_node + 1;
			Identify[2*i]=1;
		}
		if (_bc_info[i][negative_face_id] == 1)
		{
			num_negative_node = num_negative_node + 1;
			Identify[2 * i] = -1;
		}
		if (_bc_info[i][6+positive_face_id] == 1)
		{
			num_positive_node = num_positive_node + 1;
			Identify[2 * i + 1] = 1;
		}
		if (_bc_info[i][6 + negative_face_id] == 1)
		{
			num_negative_node = num_negative_node + 1;
			Identify[2 * i + 1] = -1;
		}
	}
	if (num_positive_node==0 || num_negative_node==0)
	{
		return 1e20;
	}
	int num_DOF = 1;
	for (int i = 0; i < num_electrical_node; i++)
	{
		if (Identify[i]==1)
		{
			ID[i] = 1;
		}
		else if (Identify[i] == -1)
		{
			ID[i] = 0;
		}
		else
		{
			num_DOF = num_DOF + 1;
			ID[i] = num_DOF;
		}
	}
	for (int n = 0; n < _resistance_net.size(); n++)
	{
		_resistance_net[n].i = ID[_resistance_net[n].i];
		_resistance_net[n].j = ID[_resistance_net[n].j];
	}
	if (output_control==0)
	{
		cout << "Calcuate the effective resistance..." << endl;
	}
	TResisitance_net Modified_Node_Analysis;
	Modified_Node_Analysis.input_resistance_net(_resistance_net);
	Modified_Node_Analysis.Generate_Sparse_Matrix();
	double effective_R = Modified_Node_Analysis.Calculate_effective_resistance();
	if (output_control==0)
	{
		cout << "The effective resistance of the CNT net is: " << abs(effective_R) << endl;
	}
	return effective_R;
}