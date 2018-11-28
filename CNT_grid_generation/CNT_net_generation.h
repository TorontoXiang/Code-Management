#pragma once
#include "TCNT_grid.h"
void Create_straight_CNT_net(int num_CNT, double l_CNT, double lx, double ly, double lz);
//Create a random straight CNT net
//lx,ly,lz - The length of the box
//num_CNT - The number of the CNT put in the box
//l_CNT - The length of the CNT
void Create_wavy_CNT_net(int num_CNT, double l_CNT, int n_divided, double ratio, double lx, double ly, double lz);
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