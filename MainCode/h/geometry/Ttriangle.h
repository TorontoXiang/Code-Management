#ifndef TTRIANGLE
#define TTRIANGLE

#include "data_structure.h"
#include "public_function.h"
#include <iostream>
class Ttriangle
{
public:
	Ttriangle(vec3D v1,vec3D v2,vec3D v3);
	Ttriangle(){};

	void column_integration(double& volume,double& ix,double& iy,double& iz);
	//Calculate the integration of the intersected column of the given triangle piece
	void initialize_face(vec3D v1,vec3D v2,vec3D v3);
	//Initialize a new face
private:
	vec3D _v1;vec3D _v2;vec3D _v3;//The coordinate of the triangle vertexes
	struct Snode
	{
		vec3D _v;            //Coordinate of the node
		double _variable;    //The value of the node
		int _discard;        //Whether the node is required to be discarded
		Snode* _next;         //Point to the next node
		Snode():_variable(0),_discard(0),_next(NULL){};
	};
	int _used;              //Used RAM of the ring which is less than 13
	Snode _ring[20];         //Node vector:In order to avoid new and delete operation,we give enough space to store the node
	Snode* _head;            //Point to a node of the ring

	void value_ring(int i);
	//Value every node according to the cutting plane
	//p=1 denotes the x=0 plane
	//p=2 denotes the y=0 plane
	//p=3 denotes the z=0 plane
	//p=4 denotes the x+y=1 plane
	//p=5 denotes the x+y+z=1 plane
	void get_head();
	//The head of the ring must be inside of this cutting plane, namely _variable>0
	void insert_new_node();
	//Insert new node to the ring if the values of the adjacent nodes are in different sign and mark the nodes which will be discarded
	void cut_outside_node();
	//Cut the nodes with discard=1
	void project_onto_plane();
	//Project nodes into plane z=1-x-y which are higher than this plane
	void final_integration(double& volume,double&ix,double& iy,double& iz);
	//Calculate the integrations on 1,x,y,z of the final ring
	
};
#endif
