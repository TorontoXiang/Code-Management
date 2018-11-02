#include "Shape_function.h"
const int c_brick_8[8][3]=
{
	{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
	{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}

};
double N_brick_8(int i,double xi,double eta,double zeta)
{
	return 0.125*(1+c_brick_8[i][0]*xi)*(1+c_brick_8[i][1]*eta)*(1+c_brick_8[i][2]*zeta);
}
vec3D GN_brick_8(int i,double xi,double eta,double zeta)
{
	vec3D result;
	result.x=0.125*c_brick_8[i][0]*(1+c_brick_8[i][1]*eta)*(1+c_brick_8[i][2]*zeta);
	result.y=0.125*c_brick_8[i][1]*(1+c_brick_8[i][0]*xi)*(1+c_brick_8[i][2]*zeta);
	result.z=0.125*c_brick_8[i][2]*(1+c_brick_8[i][0]*xi)*(1+c_brick_8[i][1]*eta);
	return result;
}