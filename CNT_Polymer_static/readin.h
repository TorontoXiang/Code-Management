#ifndef READIN
#define READIN
//Define the class for readin input file
#include <iostream>
#include <fstream>
#include <vector>
#include "data_structure.h"
using namespace std;
class Tbody_base;
//Tmeporary structure for brick
struct Scell_8
{
	int cell_id;
	int IEN[8];
	int mid;
};
//Temporary structure for node
struct Snode_temp
{
	int id;            //The node id
	double x,y,z;      //The node coordinate
};
struct Snode_group
{
	int group_id;                //Node group id
	vector<int> node_id;         //Node id in this group
};
struct Spart
{
	int part_id;         //Part id
	int section_id;      //Section id
	int material_id;     //Material id
	int EOS_id;          //EOS id
};
struct Smaterial
{
	int material_id;     
	double Youngs;
	double Possion;
};
struct Sregular_grid
{
	double x_min[3];
	double x_max[3];
	int nx,ny,nz;
};
struct Sregular_grid_bc
{
	string FaceName;       //FaceName=x_min,x_max,y_min,y_max,z_min,z_max
	int direction;         //Constriant direction
	double value;          //The value of the constraint
};
//---------------------------------------------------------------------
//Structure for reading keyword file
//---------------------------------------------------------------------
struct Skeyword
{
	Skeyword() {};
	vector<Snode_temp> node_list;
	vector<Scell_8> cell_8_list;
	vector<Snode_group> node_group_list;
	vector<Spart> part_list;
	vector<Smaterial> material_list;
	Sregular_grid regular_grid;
	vector<Sregular_grid_bc> regular_bc_list;

	void clear_keyword();
};
//------------------------------------------------------------------------------
//Functions in reading keyword file
//------------------------------------------------------------------------------
void read_in_keyword_file(ifstream& input,Skeyword& keyword);   //Read in the keyword file and store in keyword structure
void new_line(ifstream& input);     //Move to the beginning of the next line
bool exam_keyword(ifstream& input); //Whether the beginning of the next functional line is keyword
void next_keyword(ifstream& input); //Come to the next keyword from current position
void next_data(ifstream& input);    //Finish this line and scape the note to the next functional line

#endif