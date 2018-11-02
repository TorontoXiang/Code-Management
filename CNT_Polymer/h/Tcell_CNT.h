#ifndef TCELL_CNT
#define TCELL_CNT
#include "Tnode.h"
#include "data_structure.h"
#include <vector>
#include "TMAT_base.h"
#include <iostream>
using namespace std;
//Define cell shell
class Tcell_CNT
{
public:
	Tcell_CNT() {};
	Tcell_CNT(int cell_id, Tnode_CNT* (&node_ptr)[2], double area, TMAT_base* mat);


	void calculate_corner_force();
	//Calculate the corner force
	void assemble_corner_force();
	//Assemble the corner force into nodal force and clear the corner force
protected:
	int _id;                    //Cell id
	TMAT_base* _mat_ptr;        //Cell material
	Tnode_CNT* _node_ptr[2];    //Cell node pointer
	vec3D _corner_force[2];     //The nodal force from this cell
	double _area;               //The area of the bar
	double _l0;                 //The initial length of the bar
	int _location_id;           //The location of the CNT in the polymer
	double _ix, _iy, _iz;       //The isoparametic coordinate of the CNT cell in the polymer

	//Friend Class List
	friend class Tbody;

};
#endif
#pragma once
