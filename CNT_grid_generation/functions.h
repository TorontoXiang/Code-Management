#pragma once
#include "data_structure.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
using namespace std;
double min_dis(vec3D x1, vec3D x2, vec3D x3, vec3D x4);
//Retrun the minimal distance of segment x1x2 and x3x4
void Create_straight_CNT_net(int num_CNT,double l_CNT,double lx,double ly,double lz);
//Create a random straight CNT net
//lx,ly,lz - The length of the box
//num_CNT - The number of the CNT put in the box
//l_CNT - The length of the CNT
void Create_wavy_CNT_net(int num_CNT, double l_CNT, int n_divided, double ratio, double lx, double ly, double lz, double E, double mu);
//Create a random wavy CNT net
//lx,ly,lz - The length of the box
//num_CNT - The number of the CNT put in the box
//l_CNT - The length of the CNT
//n_divided - The division of each CNT
//ratio - control the curvature of the CNT
//E,mu - The material property of the CNT
void Generate_straight_CNT(vec3D p1, vec3D p2, int n_divided, double r, double l, int m, int n, ofstream& output_k, ofstream& output_tec, double E, double mu);
//Generate a straight CNT grid and input file if the two endpoints are given
//p1,p2 - The two endpoints
//n_divided - The number of the division along the length
//r,l - The radius and length of the CNT
//m,n - The parameters for generating the truncation mesh
//output_k,output_tec - The output of the CNT grid
//E,mu - The material property of the CNT
vector<vec3D> Generate_wavy_CNT(vec3D p_begin, double ratio, int n_divided, double l);
//p_begin - The start point of the CNT
//ratio - The parameter for controlling the curvature
//n_divided - The number of the division along the length
//r,l - The radius and length of the CNT
//m,n - The parameters for generating the truncation mesh
//output_k,output_tec - The output of the CNT grid
//E,mu - The material property of the CNT
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