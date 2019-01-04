#pragma once
#include "data_structure.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
using namespace std;
struct Selectrical_node
{
	int CNT_ID;
	int Segment_ID;
	double ratio;
	int Electric_node_ID;

	bool operator ==(Selectrical_node const &other) { return (Segment_ID + ratio == other.Segment_ID + other.ratio); };
	bool operator <(Selectrical_node const &other) { return (Segment_ID + ratio < other.Segment_ID + other.ratio); };
	bool operator >(Selectrical_node const &other) { return (Segment_ID + ratio > other.Segment_ID + other.ratio); };
};

double min_dis(vec3D x1, vec3D x2, vec3D x3, vec3D x4);
//Retrun the minimal distance of segment x1x2 and x3x4
double min_dis(vec3D x1, vec3D x2, vec3D x3, vec3D x4, double &ratio1, double &ratio2);
//Retrun the minimal distance of segment x1x2 and x3x4 and the information of the two related points
vector<vec3D> Generate_wavy_CNT(vec3D p_begin, double ratio, int n_divided, double l,int direction=0);
//p_begin - The start point of the CNT
//ratio - The parameter for controlling the curvature
//n_divided - The number of the division along the length
//l - The radius and length of the CNT
//direction - The direction of the first segment:0 - random;1 - along +z;-1 - along -z
void T_matrix(vec3D z, double(&T)[3][3]);
//Calculate the Transformation matrix between two coordinate system 1 and 2 such that x1=Tx2
//z - z is the z-axis in system 2 and its coordinate in system 1
void matrix_multiply(double(&A1)[3][3], double(&A2)[3][3], double(&A)[3][3]);
//Calculate A=A1*A2
bool is_curve_valid(vector<vec3D> &curve1, vector<vec3D> &curve2,double min_dis);
//Determine whether the minimal distance between two curves is larger than min_dis
void curve_range(vector<vec3D> &curve, vec3D &x_min, vec3D &x_max);
//Calculate the range of a curve x_min and x_max
bool effective_range(vec3D &x_min1, vec3D &x_max1, vec3D& x_min2, vec3D x_max2, vec3D &x_min, vec3D &x_max, double min_dis);
//Calculate the overlapping region of two curve
bool is_in_range(vec3D p1, vec3D p2, vec3D x_min, vec3D x_max);
//Determine whether any part of a segment is inside a range
int identify_segment_state(vec3D p1, vec3D p2, vec3D x_min, vec3D x_max);
//Determine the state of a segment in curved CNT
//Return -1 - inside RVE
//Return -2 - outside RVE
//Return 0 - cross x_min;Return 1 - cross y_min;Return 2 - cross z_min
//Return 3 - cross x_max;Return 4 - cross y_max;Return 5 - cross z_max
void move_curve(vector<vec3D> &curve, vec3D motion);
//Move a curve by motion
double min_curve_dis(vector<vec3D> &curve1, vector<vec3D> &curve2);
//Calculate the minimal distance of two curve
double min_curve_dis(vector<vec3D> &curve1, vector<vec3D> &curve2, Selectrical_node &node1, Selectrical_node &node2);
//Calcualte the minimal distance of two curve and return the information of the two related points
void effective_curve(vector<vec3D> &curve, vector<vector<vec3D>> &effective_curve_list, vec3D x_min, vec3D x_max);
//Get the effective curve in a box

void identify_straight_CNT_at_boundary(vec3D p1, vec3D p2, vector<int> &index,double lx,double ly,double lz);
//Identify whether the endpoint of a straight CNT is located at the boundary