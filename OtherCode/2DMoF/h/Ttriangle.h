#ifndef TRIANGLE
#define TRIANGLE
#include "public.h"
#include "Tpolygon.h"
struct triangle
{
	triangle(vec2D v1,vec2D v2,vec2D v3):p1(v1),p2(v2),p3(v3){};
	vec2D p1,p2,p3;
	
	void intersecting_with_polygon(Tpolygon* polygon,double& area,vec2D& moment);
	void triangle_segment_calculation(vec2D q1,vec2D q2,double& area,vec2D& moment);
};

#endif