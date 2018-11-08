#include <vector>
#include "data_structure.h"
#include <cmath>
#include <fstream>
using namespace std;
void CNT_grid_generation(vector<int> &node_list);
void generate_local_nodes(double r,double l,int n,int m,int &nume,int &nump,int**IEN,vec3D* node_list);
int calculate_index(int i, int j, int i_max2, int i_max1,int j_max1,int nump1);
void T_matrix(vec3D z, double(&T)[3][3]);
//Calculate the Transformation matrix between two coordinate system 1 and 2
//z: the coordinate of z-axis of system 2 in system 1
int main()
{
	int nume, nump;
	int** IEN=NULL;
	vec3D* node_list=NULL;
	generate_local_nodes(1, 0.8, 10,10,nume,nump,IEN,node_list);
	return 0;
}
void CNT_grid_generation(vector<int> &node_list)
{

}
void generate_local_nodes(double r, double l,int n,int m, int &nume, int &nump, int**IEN, vec3D* node_list)
{
	//Generate grid on two curve portion
	double pi = 3.141592654;
	int i_max1 = 4 * n + 4, j_max1 = m + 1;
	int i_max2 = n + 1, j_max2 = n + 1;
	int nump1 = i_max1 * (j_max1+1), nume1 = i_max1*j_max1;
	int nump2 = (i_max2 - 1)*(j_max2 - 1);
	int nume2 = i_max2 * j_max2;
	nume = nume1 + nume2; nump = nump1 + nump2;
	node_list = new vec3D[nump1 + nump2];
	IEN = new int*[nume1 + nume2];
	for (int i = 0; i < nume1+nume2; i++)
	{
		IEN[i] = new int[4];
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
		double theta = -3*pi/4+2 * pi * i / i_max1;
		pos[i][0].value(r*cos(theta), r*sin(theta), 0);
		if (i<=n+1)
		{
			pos[i][j_max1].value(-0.5*l + l * i / (n + 1), -0.5*l, 0);
		}
		else if (i>n+1 && i<=2*(n+1))
		{
			pos[i][j_max1].value(0.5*l, -0.5*l + l * (i - n - 1) / (n + 1), 0);
		}
		else if (i>2*(n+1) && i<=3*(n+1))
		{
			pos[i][j_max1].value(0.5*l - l * (i - 2*n - 2) / (n + 1), 0.5*l, 0);
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
		for (int j = 0; j < j_max1+1; j++)
		{
			int index = (i_max1)*j + i;
			node_list[index] = pos[i][j];
		}
	}
	for (int i = 0; i < i_max1; i++)
	{
		for (int j = 0; j < j_max1; j++)
		{
			int index = i_max1 * j + i;
			int n1 = i_max1*j + i, n2 = i_max1*j + (i + 1)%i_max1, n3 = i_max1*(j + 1) + (i + 1)%i_max1, n4 = i_max1*(j + 1) + i;
			IEN[index][0] = n1, IEN[index][1] = n2, IEN[index][2] = n3, IEN[index][3] = n4;
		}
	}
	//Generate the inner grid
	vec3D o(-0.5*l, -0.5*l, 0);
	for (int i = 0; i < i_max2-1; i++)
	{
		for (int j = 0; j < j_max2-1; j++)
		{
			vec3D temp((i + 1)*l / i_max2, (j + 1)*l / j_max2, 0);
			int index = (i_max2-1)*j + i + nump1;
			node_list[index] = o + temp;
		}
	}
	for (int i = 0; i < i_max2; i++)
	{
		for (int j = 0; j < j_max2; j++)
		{
			int index = i_max2 * j + i + nume1;
			int n1 = calculate_index(i, j, i_max2, i_max1, j_max1, nump1);
			int n2 = calculate_index(i+1, j, i_max2, i_max1, j_max1, nump1);
			int n3 = calculate_index(i+1, j+1, i_max2, i_max1, j_max1, nump1);
			int n4 = calculate_index(i, j+1, i_max2, i_max1, j_max1, nump1);
			IEN[index][0] = n1, IEN[index][1] = n2, IEN[index][2] = n3, IEN[index][3] = n4;
		}
	}
	ofstream output;
	output.open("temp.dat");
	output << "TITLE = \"Tecplot Grid\"" << endl;
	output << "VARIABLES = \"X\",\"Y\",\"Z\"" << endl;
	output << "ZONE F=FEPOINT,N=" << nump1+nump2<< "," << "E=" << nume1+nume2 << "," << "ET=QUADRILATERAL" << endl;
	for (int i = 0; i < nump1+nump2; i++)
	{
		output << node_list[i].x << " "<<node_list[i].y <<" "<<node_list[i].z << endl;
	}
	for (int i = 0; i < nume1+nume2; i++)
	{
		output << IEN[i][0]+1 << " " << IEN[i][1]+1 << " " << IEN[i][2]+1 << " " << IEN[i][3]+1 << endl;
	}
	return;
}
int calculate_index(int i, int j, int i_max2, int i_max1,int j_max1,int nump1)
{
	int index;
	int i1, j1 = j_max1;
	if (i==0)
	{
		i1 = (4 * i_max2 - j) % (4 * i_max2); index = i_max1 * j1 + i1;
	}
	else if (j==0)
	{
		i1 = i; index = i_max1 * j1 + i1;
	}
	else if (i==i_max2)
	{
		i1 = i_max2 + j; index = i_max1 * j1 + i1;
	}
	else if (j==i_max2)
	{
		i1 = 3 * i_max2 - i; index = i_max1 * j1 + i1;
	}
	else
	{
		i = i - 1; j = j - 1; index = (i_max2 - 1)*j + i + nump1;
	}
	return index;
}
void T_matrix(vec3D z, double(&T)[3][3])
{
	double s_t, c_t, s_p, c_p;
	double eps = 1e-20;
	if (z.z==1)
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
		T[2][0] = s_t * s_p; T[2][1] = s_t * c_p; T[3][3] = z.z;
	}
	return;
}