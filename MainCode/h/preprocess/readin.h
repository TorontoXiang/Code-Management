#ifndef READIN
#define READIN
//Define the class for readin input file
#include <iostream>
#include <fstream>
#include <vector>
#include "TMAT_base.h"
#include "TEOS.h"
#include "Tviscosity_base.h"
#include "Tinteraction_base.h"
using namespace std;
class Tbody_base;
//Temporary structure for shell
struct Scell_4
{
	int cell_id;         //The id of the shell
	int IEN[4];          //The connected node for the cell
	int part_id;         //The part id of the cell
};
//Tmeporary structure for brick
struct Scell_8
{
	int cell_id;
	int IEN[8];
	int part_id;
};
//Brick cell with variable information
struct Scell_8_with_variable : public Scell_8
{
	double density, internal_energy;
};
//Temporary structure for node
struct Snode_temp
{
	int id;            //The node id
	double x,y,z;      //The node coordinate
};
//Node with velocity information
struct Snode_with_vel : Snode_temp
{
	double vx, vy, vz;
};
//Temproary structure for particle
struct Sparticle_temp
{
	int group_id;
	int part_id;
	vec3D position;
	vec3D velocity;
	double volume;
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
struct Ssection_shell
{
	int section_id;      //Section id
	int nGauss;          //The number of Gauss point in the shell
	double thickness;    //The initial thickness of the shell
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
struct SEOS
{
	int EOS_id;
	string EOS_type;
	double c0,c1,c2,c3,c4,c5,c6;
	double internal_energy_per_volume;
	double p0,gama,b;
	double A,B,R1,R2,w;
	double gama0,s;
};
struct Scurve
{
	int curve_id;
	vector<double> v1;      //The first column
	vector<double> v2;      //The second column
};
struct Stime
{
	double endtime;        //The terminal time
	double CFL;            //The CFL number
};
struct Sarti_vis
{
	int av_type;
	double k1,k2;
};
struct Sinteraction
{
	string interaction_name;
	Tbody_base* body1;
	Tbody_base* body2;

	string calculate_interaction_type_id();
	//Return the interaction type id
};
struct Sremapping_scheme
{
	string name;
	int num_remapping;
	double detonation_time;
	double tolerance_angle;
	double tolerance_length;
	double DtReductionRatio;
	double initialDt;
	vector<double> remapping_time_point;


	void exam_remapping_scheme();
	//Exam whether the remapping scheme is valid
};
struct SNodePerturbation
{
	double x,y,z;
	double ratio;
	int id;
	int dimension;
};
struct SImplosion
{
	double d0;
	double N;
	double r_in,r_out;
};
//Temporary structure of MPM background
struct Sbackground_temp
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
struct SMPM_boundary_condrition
{
	int strat_subscript[3];
	int terminate_subscript[3];
	int constraint[3];
};
struct SMPM_background_load
{
	int strat_subscript[3];
	int terminate_subscript[3];
	vec3D load;
};
struct SMPM_particle_load
{
	int part_id;
	int particle_id;
	vec3D load;
};
struct Sbackground_property
{
	int part_id;
	int nx_begin,ny_begin,nz_begin;
	int nx_end,ny_end,nz_end;
};
//---------------------------------------------------------------------
//Structure for reading keyword file
//---------------------------------------------------------------------
struct Skeyword
{
	Skeyword() { num_particle_part = 0; is_implosion = false; hourglass_option = true; is_initial_vel = false; };
	vector<Snode_temp> node_list;
	vector<Snode_with_vel> vel_node_list;
	vector<Sparticle_temp> particle_list;
	vector<Scell_4> cell_4_list;
	vector<Scell_8> cell_8_list;
	vector<Scell_8_with_variable> variable_cell_8_list;
	vector<Sboundary_type> boundary_list;
	vector<Snode_group> node_group_list;

	vector<Sload_on_node> load_list;
	vector<Spart> part_list;
	vector<Ssection_shell> section_shell_list;
	vector<Smaterial> material_list;
	vector<SEOS> EOS_list;
	vector<Scurve> curve_list;
	Stime time_control;
	Sbackground_temp background;
	vector<Sinitial_velocity_temp> initial_velocity_list;
	vector<SMPM_boundary_condrition> MPM_boundary_condition_list;
	vector<SMPM_background_load> MPM_background_load_list;
	vector<SMPM_particle_load> MPM_particle_load_list;
	vector<Sbackground_property> ALEMPM_background_property;
	int num_particle_part;
	bool is_implosion;
	bool hourglass_option;
	SImplosion im;
	SNodePerturbation im_node_perturbation;
	bool is_initial_vel;
	vec3D initial_vel;

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

TEOS_base* generate_EOS(SEOS EOS,Smaterial mat);
//Generate an EOS pointer

TMAT_base* generate_material(Smaterial mat,SEOS EOS,int type);
//Generate a material pointer
//mat: the information of the material
//type: 0:for 3D 1:for shell

Tviscosity_base* generate_AV(Sarti_vis arti_vis);
//Generate a artificial viscosity type

Tinteraction_base* generate_interaction(Sinteraction interaction);
//Generate a interaction type

template<class T1,class T2> Tinteraction_base* generate_interaction_with_type(Sinteraction interaction);
//Generate a interaction type with specific type
#endif