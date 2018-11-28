#pragma once
#include "data_structure.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
using namespace std;
double min_dis(vec3D x1, vec3D x2, vec3D x3, vec3D x4);
//Retrun the minimal distance of segment x1x2 and x3x4
vector<vec3D> Generate_wavy_CNT(vec3D p_begin, double ratio, int n_divided, double l,int direction=0);
//p_begin - The start point of the CNT
//ratio - The parameter for controlling the curvature
//n_divided - The number of the division along the length
//l - The radius and length of the CNT
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
void move_curve(vector<vec3D> &curve, vec3D motion);
//Move a curve by motion
double min_curve_dis(vector<vec3D> &curve1, vector<vec3D> &curve2);
//Calculate the minimal distance of two curve
void effective_curve(vector<vec3D> &curve, vector<vector<vec3D>> &effective_curve_list, vec3D x_min, vec3D x_max);
//Get the effective curve in a box