#ifndef  TTETRAHEDRON
#define  TTETRAHEDRON

#include "data_structure.h"
#include "Ttriangle.h"
#include "Tpolyhedron.h"
class Ttetrahedron
{
public:
	Ttetrahedron(){};
	Ttetrahedron(vec3D p1,vec3D p2,vec3D p3,vec3D p4);
	Ttetrahedron(vec3D p1,vec3D p2,vec3D p3,vec3D p4,int temp);  //Constructor for calculating distortion


	void intersect_with_polyhedron(Tpolyhedron* mp,double &ivolume,double &ix,double &iy,double &iz);
	//Calculate the volume and moment of the intersecting portion between a tetrahedron and a polyhedron
	void calculate_vertex_range(vec3D &coor_min,vec3D &coor_max);
	//Calculate the coordinate range of this tetrahedron

	//Access functions
	double G_volume(){return _volume;};
	vec3D G_centroid(){return (_v1+_v2+_v3+_v4)*0.25;};
	double G_distort_ratio(){return _distort_ratio;};
	vec3D G_vertex(int i);
private:
	vec3D _v1,_v2,_v3,_v4;                //vertexes of the tetrahedron
	double _transfer_matrix[3][3];        //transfer matrix of the tetrahedron
	double _volume;                       //volume of the tetrahedron
	double _distort_ratio;                //The distort ratio of the tetrahedron

	void transfer_vertex(vec3D &p);
	//Transfer a point into the standard space
	vec3D transfer_result(vec3D &p,double &ivolume);
	//Transfer the volume integration on x,y,z in standard space into physical space
};
#endif