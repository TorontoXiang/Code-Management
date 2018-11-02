#ifndef READIN
#define READIN
//Define the class for readin input file
#include <iostream>
#include <fstream>
#include <vector>
#include "TMAT_base.h"
#include "data_structure.h"
using namespace std;
class Tbody_base;
//Tmeporary structure for brick
struct Scell_8
{
	int cell_id;
	int IEN[8];
	int part_id;
};
struct Scell_2
{
	int cell_id;
	int IEN[2];
	int part_id;
};
//Temporary structure for node
struct Snode_temp
{
	int id;            //The node id
	double x,y,z;      //The node coordinate
};
//Temporary structures for control cards
struct Sboundary_type
{
	int id;            //Boundary condition will be applied on temp_node_list[id]
	int type;          //0-fixed boundary 1-velocity boundary
	int pos[3];        //Constraint in position
	int rot[3];        //Constraint in rotation
	int dofvel[3];     //Constraint in velocity
	double vel[3];     //Velocity boundary condition
	double x,y,z,v;    //Velocity boundary for cylinder and sphere
	
};
struct Snode_group
{
	int group_id;                //Node group id
	vector<int> node_id;         //Node id in this group
};
struct Sload_on_node
{
	int id;          //Load will be applied on temp_node_list[id]
	int dof;         //1,2,3: Fx,Fy,Fz;4,5,6: Mx,My,Mz 
	int curve_id;    //The Load curve
};
struct Spart
{
	int part_id;         //Part id
	int section_id;      //Section id
	int material_id;     //Material id
	int EOS_id;          //EOS id
};
struct Ssection
{
	int scetion_id;
	double area;          //The area of the section
};
struct Smaterial
{
	int material_id;         //Material id
	string material_type;    //The type of the material
	double density;
	double Youngs;
	double Possion;
	double Sigma_y;
	double ET;
	double A,B,n,C,epso;    //Parameters for Johnson Cook model
};
struct Stime
{
	double endtime;        //The terminal time
	double CFL;            //The CFL number
	double motiontime;     //The time for applying the displacement
};
struct Sarti_vis
{
	int av_type;
	double k1,k2;
};
struct Sregular_grid
{
	vec3D x_min;
	vec3D x_max;
	int nx,ny,nz;
};
struct Sinitial_velocity_temp
{
	int part_id;
	vec3D velocity;
};
struct Scurve
{
	int curve_id;
	vector<double> v1;      //The first column
	vector<double> v2;      //The second column
};
struct Sregular_grid_bc
{
	string FaceName;       //FaceName=x_min,x_max,y_min,y_max,z_min,z_max
	int bc_type;           //0-fix position,1-fix velocity
	int direction;         //Constriant direction
	double velocity;       //The given velocity (Only for velocity boundary condition)
};
//---------------------------------------------------------------------
//Structure for reading keyword file
//---------------------------------------------------------------------
struct Skeyword
{
	Skeyword() { is_initial_vel = false; is_regular_grid = false; };
	vector<Snode_temp> node_list;
	vector<Scell_8> cell_8_list;
	vector<Scell_2> cell_2_list;
	vector<Sboundary_type> boundary_list;
	vector<Snode_group> node_group_list;
	vector<Scurve> curve_list;
	vector<Sload_on_node> load_list;
	vector<Spart> part_list;
	vector<Smaterial> material_list;
	vector<Sinitial_velocity_temp> initial_velocity_list;
	vector<Ssection> section_list;

	Stime time_control;
	bool is_initial_vel;
	vec3D initial_vel;
	bool is_regular_grid;

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

TMAT_base* generate_material(Smaterial mat,int type);
//Generate a material pointer
//mat: the information of the material
//type: 0:for 3D 1:for shell
#endif