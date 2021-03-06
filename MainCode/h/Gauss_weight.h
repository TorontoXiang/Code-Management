#ifndef GAUSS_WEIGHT
#define GAUSS_WEIGHT
const double Gauss_weight[5][5]=
{
	{              2.0,               0.0,               0.0,               0.0,               0.0},
	{              1.0,               1.0,               0.0,               0.0,               0.0},
	{          5.0/9.0,           8.0/9.0,           5.0/9.0,               0.0,               0.0},
	{0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454,               0.0},
	{0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189}
};
const double Gauss_coordinate[5][5]=
{
	{               0.0,                0.0,               0.0,               0.0,               0.0},
	{-0.577350269189626,  0.577350269189626,               0.0,               0.0,               0.0},
	{-0.774596669241483,                0.0, 0.774596669241483,               0.0,               0.0},
	{ -0.86113631159405,  -0.33998104358486,  0.33998104358486,  0.86113631159405,               0.0},
	{-0.906179845938664, -0.538469310105683,               0.0, 0.538469310105683, 0.906179845938664}
};
#endif // !GAUSS_WEIGHT
