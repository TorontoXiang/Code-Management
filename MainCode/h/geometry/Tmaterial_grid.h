#pragma once
#ifndef TMATERIAL_GRID
#define TMATERIAL_GRID
//Define the material grid
#include "Ttetrahedron.h"
#include "TMAT_base.h"
#include <iostream>
#include <fstream>
#include "readin.h"
using namespace std;
class Tmaterial_grid
{
public:
	struct SIEN
	{
		int IEN[4];
	};
	Tmaterial_grid(int mid,bool is_pseudo):_nume(0),_total_volume(0),_material_id(mid),_is_peusdo(is_pseudo){_total_centroid.value(0,0,0);};
	void input_material_grid(ifstream &input);
	//Input the material grid
	void calculate_geometry_property();
	//Calculate the total volume and centroid of this material grid
	void plot_material_grid(ofstream &output);
	//Plot the material grid

	int _material_id;
	int _nume;   
	int _nump;
	double _total_volume;
	bool _is_peusdo;           //Whether the material surface is a peusdo surface
	vec3D _total_centroid;
	vector<SIEN> _IEN_list;
	vector<vec3D> _node_list;
	vector<Ttetrahedron> _material_tet_list;
	Smaterial _material;
	SEOS _EOS;
};
#endif