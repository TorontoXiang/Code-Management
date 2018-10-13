#ifndef TPOLYGON
#define TPOLYGON
//Define the class for the polygon
#include "public.h"
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
struct offset
{
	offset():_offset(0),_span(-1),_id(-1){};
	double _offset;          //The offset of a vertex
	double _span;            //The span of this vertex
	int _id;                 //The corresponding vertex id

	bool operator <(offset const &other) const {return _offset>other._offset;};
	bool operator >(offset const &other) const {return _offset<other._offset;};
	bool operator ==(offset const &other) const {return _offset==other._offset;};
};
class Tpolygon
{
public:
	Tpolygon(int num_edge):_num_vertex(num_edge){_vertex_list.resize(num_edge);_offset_list.resize(num_edge);};
	Tpolygon(){};
	void Add_vertex(vec2D& new_vertex,int id);
	//Add new vertex to the polygon
	void calculate_offset_and_span(vec2D& n);
	//Calculate the offset of each vertex with respected to n
	void calculate_properities();
	//Calculate the area and the centroid of a polygon
	void flood_algorithm(double fraction,double& offset_solution,double& length);
	//Solve the volume equation by the flood algorithm
	//n & fraction: the input normal and given volume fraction
	//offset_solution: the offset of the target solution
	//length: the length of the surface
	void rotation_algorithm(double fraction,vec2D& centroid_ref,int id,discontinuous_point (&dp)[2]);
	//Calculate the discontinuity theta for vertex i with the given fraction
	//frcation & id: the given frcation and the target vertex
	//dp: the information of the two discontinuous points
	//centroid_below:the centroid below the surface
	void calculate_geometry_for_derivatives(vec2D& n,double offset_solution,double& area_below,vec2D& centroid_below);
	//Calculate the variables for calculating the derivatives
	//centroid_below: the centroid below the surface
	void divide_polygon(vec2D& n,double fraction,Tpolygon* &below,Tpolygon* &above);
	//Divided the polygon for the given n and fraction
	void plot_polygon(ofstream& ouput,double value);
	//Output the polygon with a value
	int trial_id_base(vec2D& n,double fraction,int& order);
	//Calculate the initial trial id_base

	void calculate_range(vec2D& x_range,vec2D& y_range);

	void S_num_edge(int num){_vertex_list.resize(num);_offset_list.resize(num);_num_vertex=num;};
	void S_vertex(vec2D p,int i){_vertex_list[i]=p;};
	double G_area(){return _area;};
	int G_num_vertex(){return _num_vertex;};
	vec2D G_centroid(){return _centroid;};
	vec2D G_vertex(int i){return _vertex_list[i];};

private:
	int _num_vertex;                  //The edge number of the polygon
	vector<vec2D> _vertex_list;    //The coordinate list of the vertex in clock order
	vector<offset> _offset_list;   //The offset list of the vertex in incresing order
	double _area;
	vec2D _centroid;

	double calculate_span(int id);
	//Calcualte the span of vertex id

	friend class TMoF2D;
};


#endif
