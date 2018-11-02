#ifndef SHAPE_FUNCTION
#define SHAPE_FUNCTION
#include "data_structure.h"

//Define the shape functions for 8 node brick element
double N_brick_8(int i,double xi,double eta,double zeta);
vec3D GN_brick_8(int i,double xi,double eta,double zeta);

#endif