#pragma once
#include "data_structure.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
using namespace std;
double min_dis(vec3D x1, vec3D x2, vec3D x3, vec3D x4);
//Retrun the minimal distance of segment x1x2 and x3x4
void Create_straight_CNT_net(int num_CNT,double l_CNT,double lx,double ly,double lz);
//Create a random CNT net
//lx,ly,lz - The length of the box
//num_CNT - The number of the CNT put in the box
//l_CNT - The length of the CNT