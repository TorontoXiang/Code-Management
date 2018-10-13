#ifndef TPOLYHEDRON
#define TPOLYHEDRON
//Define the class for arbitrary polyhedron
//In order to represent the non-convex polyhedrons, the polyhedron is made up with vertexes and several triangle pieces
#include "data_structure.h"
#include <vector>
#include <string>
#include "public_function.h"
#include <fstream>
using namespace std;
//Define the structure of the cutting ring (the polygon generated from the approximate surface) in Tpolyhedron
struct Scutting_ring
{
	vector<vec3D> _point_list;   //Store the cutting point in a vector which will be changed in every iteration step
	                             //These points are in the order of outward normal of the approxiamte surface
	int _length;                 //The vertexes in _point_list
	int _num_ring;               //The number of cutting ring in this structure
	                             //One approxiamte may generate more than one rings if the polyhedron is non-convex
	int _strat_point_id[10];     //The start position of every ring in _point_list
	                             //Assume that the maximal number of the rings is 10 by one approximate surface

	Scutting_ring():_length(0),_num_ring(0){};
	void add_point(vec3D &new_point);   
	//Add a new cutting point to ring
	//If length=size,use push_back
	void final_one_ring();
	//Final one ring if more than one rings exist
	void calcualte_volume_integral(vec3D &n,double &s,double &int_1,vec3D &int_x,string type);
	//Calculate the integral on the ring(s) for volume(type="volume") and centroid(type="centroid")
	//n:the outward normal of the cutting plane
	//s:the cutting ring(s) area
	//int_1:integral on z for volume
	//int_x:intrgral on zx,zy,zz for centroid
	vec3D calculate_surface_integral(string type);
	//Calculate the integral on the cutting ring
	//type="centroid":calculate the ring centroid
	//type="derivative":calculate the needed variables for analytical derivative
	void determine_ring_order(vec3D &n,int &order,vec3D_int &ring_centroid);
	//Determine the order of this ring (when calling this function, the structure only contains one ring)
	//order=0:along the direction of n;order=1:opposite to the direction of n
	//ring_centroid:the centroid of this ring
	void calculate_local_coordinate_system(double theta,double phy,vec3D &l1,vec3D &l2,vec3D &l3);
	//Calculate the local coordinate system in analytical derivative calculation
	void clear_ring();
	//Clear this ring
	void print_points();
	//Print the points of the ring onto the screen

};
class Tpolyhedron
{
public:
	Tpolyhedron():_numv(0),_nump(0),_volume(0){};

	void add_vertex(vec3D_double &new_vertex);  
	void add_vertex(vec3D &new_coordinate);
	//Add a vertex into the vertex_list
	void add_piece(int n1,int n2,int n3);   
	//Add a triangle piece into the piece_list
	void clear_the_polyhedron();
	//Clear the vertex and pieces of the polyhedron but remain the memory
	double calculate_polyhedron_volume();
	//Calculate the volume of the polyhedron
	vec3D calculate_polyhedron_centroid();
	//Calculate the volume and centroid of the polyhedron
	vec3D_double calculate_polyhedron_center();
	//Calculate the center of the polyhedron
	void calculate_vertex_range(vec3D &coor_min,vec3D &coor_max);
	//Calculate the range of the vertexes in this polyhedron
	void plot_this_polyehdron(ofstream &output,double value=0);
	//Plot the polyhedron
	//value: The value on each vertex (for ploting)
	void move_polyhedron(vec3D &n,double d);
	//Move the polyehdron along n

	//--------------------------------------------------------------------
	//Operations in surface reconstrcution phase
	//The equation of the plane is nx-d=0
	//---------------------------------------------------------------------
	void calculate_piece_topology();                    
	//Calculate adjacent of every piece 
	void calculate_sign_distance(vec3D &n,double d);
	//Set _variable(the sign distance to the surface n*x+d=0)for every vertex
	vec3D calculate_fraction_below(vec3D &n,double d);    
	//Calculate volume of the polyhedron below the cutting plane (first component of the result)
	//Calculate the derivative with respect to d (the area of the cutting ring) ( second component of the result)
	vec3D calcualte_centroid_below(vec3D &n,double d);      
	//Calculate center of the polyhedron below the cutting plane
	void cut_polyhedron_by_plane(vec3D &n,double d,Tpolyhedron* poly_below,Tpolyhedron* poly_above);   
	//Calculate two sub_polyhedron by the optimal cutting plane
	void calculate_altitude_range(vec3D &n,double &d_min,double &d_max);
	//Calculate the altitude range with respected to a given normal
	//-----------------------------------------------------------------------

	//Access functions
	vec3D G_vertex(int i){return _vertex_list[i]._coor;};
	vec3D_double* G_vertex_info(int i){return &_vertex_list[i];};
	int G_numv(){return _numv;};
	int G_nump(){return _nump;};
	int* G_IEN(int i){return _piece_list[i]._IEN;};
	vec3D G_vertex_on_piece(int piece_id,int IENid){return _vertex_list[_piece_list[piece_id]._IEN[IENid]]._coor;}
	double G_volume(){return _volume;};
	vec3D G_centroid(){return _centroid;};
	Scutting_ring* G_cutting_ring(){return &_cutting_ring;};
	int G_num_surface_piece(){return _surface_piece_list.size();};
	int G_surface_piece_id(int i){return _surface_piece_list[i];}

	//Value functions
	void S_vertex(int i,vec3D &coor){_vertex_list[i]._coor=coor;};
private:
	struct Spiece                     //Triangle piece of the polyhedron
	{
		int _IEN[3];                  //The adjacent triangle
		int _adjacent[3];             //The adjacent piece (edge vertex_id 0:0-1;1:1-2;2:2-0)
		int _type[2];                 //Piece type 
		                              //first component:0-cross;1-above;2-below 
		                              //p.s:If there are two points in a piece are locating on the cutting plane, this piece is across
		                              //second component:0-general;i-[i-1]point degenerate(i=1,2,3);i-[i-4]edge degenerate (i=4,5,6)
		int _is_visit;                //whether the piece is visited 0-not;1-yes;
		Spiece(){};
		Spiece(int n1,int n2,int n3);
	};
	vector<vec3D_double> _vertex_list;   //The vertexes of the polyhedron
	vector<Spiece> _piece_list;          //The triangle pieces of the polyehedron
	int _numv;                           //The number of vertexes
	int _nump;                           //The number of triangle pieces
	double _volume;                      //The polyhedron volume
	vec3D _centroid;                     //The polyhedron centroid
	Scutting_ring _cutting_ring;         //The cutting by an approximate surface
	vector<int> _surface_piece_list;     //The pieces of material surface (The pieces from cutting plane)


	//-------------------------------------------------------------
	//Private functions in surface reconstruction phase
	//-------------------------------------------------------------
	void classify_pieces(vec3D &n,double d); 
	//Classify the type of piece according to its vertexes' sign distance
	vec3D calculate_integral_below(vec3D &n,double d,string type);
	//Calculate integral of the polyhedron below the cutting plane
	//type="volume":below volume is stored in the first component of the result and the cutting plane area is stored in the second component
	//typr="centroid":the result is the centroid below the cutting plane
	void calculate_initial_cutting_point(int i,vec3D_int &initial_point,int &next_piece_id);   
	//Calculate the initial cutting point and the piece id next across piece
	//i:the initial cross piece
	//initial_point:the initial point of this cutting ring
	//next_piece:the next across piece id
	void calculate_next_cutting_point(int current_piece_id,vec3D_int &pre_point,int &next_piece_id,vec3D_int &next_point);
	//Calculate the next cutting point and the next cross piece id
	//current_piece_id: the piece is being operated
	//pre_point: the pre-cutting point
	//next_piece_id: the next cross piece id
	//next_point: the next cutting point
	void calculate_integral_on_subpiece(vec3D_int &pre_point,vec3D_int &next_point,int &current_piece_id,double &int_1,vec3D &int_x,string type);
	//Calculate the integration of the below sub-piece when the pre-point and next point are obtained
	void calculate_subpiece_point_below(vec3D_int &pre_point,vec3D_int &next_point,int &current_piece_id,vec3D (&point_list)[6],int &num_point);
	//Calculate the subpiece below the cutting plane
	//point_list:the points in outward order less than 6 points in the list
	//num_point: the number of point in point list
	//num_point=3:the subpiece is triangle
	//num_point=4:the subpiece is quadrangle
	//num_point=6:current piece is two point degenerate and the bilateral pieces are all below the cutting plane (almost do not occur)  
	void allocate_subpieces(vec3D_int &pre_point,vec3D_int &next_point,int current_piece_id,int id_pre,int id_next,Tpolyhedron* poly_below,Tpolyhedron* poly_above);
	//Add the sub-pieces information into the below and above polyhedron
	void classify_new_point(vector<int> &id_in_ring,vector<vec3D> &new_point_list,vec3D_int &new_point_in_ring);
	//Classify the type of a newly generated point on cutting ring 
	//Determine whether its a new vertex and its global vertex id
	void allocate_pieces_in_ring(int ring_start_position,vector<int> &id_in_ring,Tpolyhedron* poly_below,Tpolyhedron* poly_above,int order);
	//Allocate the pieces in cutting ring into poly_below and poly_above
	void compress_vertex_id(vector<int> &id_in_ring,vector<vec3D> &new_point_list,Tpolyhedron* ploy_ptr);
	//Compress the vertex id in the generated ploy_ptr from cutting process
	int calculate_common_edge_id(int &piece_id1,int &piece_id2,int &j);   
	//If _piece_list[piece_id2]._adjacent[j]=piece_id1,
	//calculate * such that _piece_list[piece_id1]._adjacent[*]=piece_id2
	//Used to determine the cutting point type
	int find_certain_piece_id(int &piece_id1,int &piece_id2,int &_vertex_id);
	//Find a piece id such that id!=piece_id1 && id!=piece_id2 && _type[0]=0
	//&& one of a vertex id in this piece is _vertex_id
	//Used to find the next piece of a degenerated cutting point
	int find_vertex_local_id(int &piece_id,int &vertex_id);                                      
	//Find local id * such that _piece_list[piece_id]._IEN[*]=vertex_id
};
#endif