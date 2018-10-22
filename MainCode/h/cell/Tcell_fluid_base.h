#ifndef TCELL_FLUID_BASE
#define TCELL_FLUID_BASE
//Define the base class of fluid cell
//Containing the public portions of pure cell and mix cell
#include "Tcell_base.h"
#include "Tnode_fluid.h"
#include "TMAT_base.h"
#include "Tpolyhedron.h"
#include "Tviscosity_base.h"
#include "TMPM_particle.h"
class Tnode_fluid;
class Tviscosity_base;
class Tcell_fluid_base : public Tcell_base
{
public:
	//Define the fluid cell internal connection matrix
	static int cell_face[6][4];       //The local index of 6 cell faces
	static int cell_edge_v[12][2];    //The two vertexes of 12 edges
	static int cell_edge_f[12][2];    //The two faces vertical to 12 edges
	static int cell_edge_f1[12][2];   //The two faces include 12 edges
	static int cell_corner_v[8][3];   //The three vertexes connected to 8 vertexes
	static int cell_corner_f[8][3];   //The three faces connected to 8 vertexes
public:
	Tcell_fluid_base(){};
	Tcell_fluid_base(int cell_id,int num_material,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat);
	//Generate a fluid cell
	//num_material:the maximal number of material in cell;1-Tcell_pure_fluid;>1-Tcell_mix_fluid
	//mat:the material pointer of the first material (background material) and other materials will
	//be read in later
	Tcell_fluid_base(int cell_id,int num_material,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat[]);
	//mat[]:All the material will be read in this function
//=======================================================================
//Virtual functions
//=======================================================================
	virtual void output_cell_monitor(ofstream& output,double t){};
//-----------------------------------------------------------------------
//Functions for calculating corner force
//-----------------------------------------------------------------------
protected:
	virtual bool calculate_corner_stress_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3],double dt){return true;};
	//Calculate the corner force due to pressure
	//Closure model will be involved in Tcell_mix_fluid
	//Return whether it is failed in calculating corner_force from stress
	virtual void calculate_corner_hourglass_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3]){};
	//Calculate the corner force due to hourglass
	//Multipul materials will be considered in Tcell_min_fluid
	virtual void set_pressure(double t){};

//--------------------------------------------------------------------
//Functions for physical
//--------------------------------------------------------------------
public:
	virtual void calculate_corner_mass(){};
	//Calculate the corner mass by corner volume and density
	//Tcell_pure_fluid only use one material
	//Tcell_mix_fluid use multipul materials
	virtual bool update_state(double dt){return true;};
	//Update the material state at corrector step
	virtual void reset_state_after_remapping(){};
	//Reset the cell state after remapping phase
	virtual void calculate_average_variable(double  &av_density,double &av_pressure,double &av_soundspeed){};
	//Calculate the average density,pressure and soundspeed
	virtual void assemble_corner_inertance(){};
	//Assemble the corner inertance into node
	virtual double calculate_cell_internal_energy(){return 0;};
	//Calculate the internal energy in the cell
	virtual void reset_material_initial_state(int mid){};
	//Set the initial state of material mid after overlapping
	virtual double calculate_time_step_particle(){return 100000;};
	//Calculate the time step for particles
//------------------------------------------------------------------------
//Functions in ALE process
//------------------------------------------------------------------------
public:
	virtual void calculate_gradient_in_cell(){};
	//Calculate the gradient of variables
	//Tcell_pure_fluid use the information of cell polyhedron
	//Tcell_mix_fluid use the information of material polyhedron
	virtual void limit_gradient_in_cell(){};
	//Limit the gradient by BJ limiter
	//Tcell_pure_fluid use the information of cell vertexes
	//Tcell_mix_fluid use the information of material polyhedron vertexes
	virtual void clear_cell_variable();
	//Clear the variable on the old cell after remapping
	//More information will be cleared for Tcell_mix_fluid
	virtual void accumulate_variables(int mid,double volume,vec3D int_x,double mass,double internal_energy){};
	//Accumulate the variables of material from the intersection portion
	virtual void reset_geometry_initial_condition(){};
	//Reset the initial material fraction and material centroid of mixed cell after overlapping
	virtual void reset_conservative_variable_initial_condition(){};
	//Reset the initial material mass and material internal energy of mixed cell before time integration

//------------------------------------------------------------------------
//Functions for multi-material
//------------------------------------------------------------------------
	virtual void update_material_centroid(){};
	//Update the material centroid before surface reconstruction
	virtual void delete_material_polyhedron(){};
	//Delete the material polyhedron after remapping
	virtual int is_surface_cell(int mid){return 0;};
	//Determine whether it is a surface cell of material mid
	//If it is,return the number of surface piece of this cell
	virtual void add_surface_cell_surface(int mid,int** piece_list,int vertex_id_start,int &piece_id){};
	//Add the surface piece of a surface cell into piece_list
	//mid:the material id
	//piece_list:the current piece_list
	//vertex_id_start:the start global vertex id of this cell
	//piece_id:input the start global piece id of this cell and output terminate global piece id of this cell
//-------------------------------------------------------------------------
//Functions for ALEMPM
//-------------------------------------------------------------------------
public:
	virtual void distribute_particle_force(TMPM_particle* particle_ptr){};
	//Distribute the particle force to the cell node
	virtual void distribute_particle_variables(TMPM_particle* particle_ptr,vec3D interval){}
	//Distribute the particle mass to the cell node
	virtual void update_particle_DOF(TMPM_particle* particle_ptr,double dt){};
	//Update the DOF of particles
	virtual void update_particle_stress(TMPM_particle* particle_ptr,double dt){};
	//Update the state of particles
protected:
	void calculate_DN_matrix(vec3D (&DN)[8],double &detJ,double xi,double eta,double zeta);
	//Calculate the dervative of shape functions with respect to x,y,z at particle
	void calcualte_Jacobi(double (&Jacobi)[3][3],vec3D (&GN)[8],double xi,double eta,double zeta);
	//Calculate the GN matrix and Jacobi matrix at particle
	virtual vec3D calculate_standard_coordinate(vec3D& particle_coordinate,vec3D interval){vec3D temp;return temp;};
	//Calculate the standard coordinate of a particle in this cell
//------------------------------------------------------------------------
//Access functions
//------------------------------------------------------------------------
public:
	virtual Tpolyhedron* G_material_polyhedron(int i){return NULL;};
	//Return the material polyhedron of material i
	//Return _cell_polyhedron for Tcell_pure_fluid
	virtual vec3D G_material_centroid(int i){vec3D temp;return temp;};
	//Return the material centroid of material i
	//Return cell centroid for Tcell_pure_fluid
	virtual vec3D G_g_ro(int i){vec3D temp;return temp;};
	//Return the density gradient of material i
	virtual vec3D G_g_roe(int i){vec3D temp;return temp;};
	//Return the density*internal of material i
	virtual double G_material_fraction(int i){return 0;};
	//Return the material fraction of material i
	virtual double G_max_spd_particle(){return 0;};
	//Return the maximal soundspeed of particles in cell
	virtual double G_average_pressure(){return 0;};
//------------------------------------------------------------------------
//Value functions
//------------------------------------------------------------------------
	virtual void S_material_polyhedron(int i,Tpolyhedron* poly_ptr){};
	//Set the material polyhedron for material i
	virtual void S_mass(int i,double mass){};
	//Set the mass for material i
	virtual void S_max_spd_particle(double max_spd){};
	//Set the maximal soundspeed of particles in cell
//========================================================================
//Common functions for Tcell_pure_fluid and Tcell_mix_fluid
//========================================================================
//------------------------------------------------------------------------
//Functions for physical
//------------------------------------------------------------------------
public:
	void calculate_corner_force(Tviscosity_base* Arti_Vis,double dt);
	void calculate_corner_force_without_hourglass(Tviscosity_base* Arti_Vis,double dt);
	//Calculate the corner in the cell
	void calculate_viscosity_force();
	//Calculate the fluid viscosity force by FEM approach
	virtual void assemble_corner_force();
	//Assemble the corner force into nodal force
	void calculate_velocity_gradient(vec3D &dvx,vec3D &dvy,vec3D &dvz);
	//Calculate the graident of velocity
	void assemble_nodal_mass_moment();
	//Assemble mass and moment to node after remapping phase
	double calculate_time_step();
	//Calculate the time step
	void reset_TR_conditions(string type);
	//Reset the initial conditions for 2D and 3D Taylor-Rayleigh problem
	void reset_dam_breaking_conditions(double dz,int nz,int mat_id,int num_material);
	//Reset the initial conditions for dam breaking problem
//------------------------------------------------------------------------
//Functions for geometry
//------------------------------------------------------------------------
public:
	void create_cell_polyhedron();
	//Create the cell polyhedron
	void update_cell_polyhedron();
	//Update the cell polyhedron before remapping phase	
	void calculate_coordinate_range();
	//Calculate _coor_min and _coor_max of the cell 
	double calculate_maximal_edge();
	//Calculate the maximal edge of the cell
	double calculate_cell_volume(string type,bool &negative);
	//Calculate the cell volume and corner volume
	//type="Full":calcualte volume at time step
	//type="Temporary":calcualte volume at predictor step
	void calculate_distorted_ratio(double &ratio_angle,double &ratio_length);
	//Calculate the angle and length distorted ratio
protected:
	void calculate_auxiliary_area(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3]);
	//Calculate the assistant and surfacial direction area
	double calculate_characteristic_length();
	//Calculate the characteristic length
	void calculate_segment_center(vec3D &c_cell,vec3D (&c_face)[6],vec3D (&c_edge)[12],string type);
	//Calculate the center of cell,face and edge

//--------------------------------------------------------------------------
//Functions for shape function interpolation
//--------------------------------------------------------------------------
public:
	vec3D calculate_std_coordinate(vec3D physical_coordiante);
	//Calculate the coordinate in the standard space of the cell
	vec3D shape_function_interpolate(vec3D std_coordinate,string name);
	//Interpolate the position or velocity by shape function
protected:
	void iteration_for_std_coordinate(vec3D &iterator_coordinate,vec3D &physical_coordinate);
	//The iteration for calculating the standard coordiante
	//Iterator_coordiante will be convergent to the result
	double get_iteration_matrix(int i,int j,vec3D &iterator_coordinate);
	//Calculate the iteration matrix
	
//-----------------------------------------------------------------------
//Functions for cell topology
//-----------------------------------------------------------------------
public:
	int vertex_reverse(int nid);
	//Calculate * s.t _node_ptr[*]->_id=nid
	int cell_corner_v_reverse(int p0,int p1);
	//Calculate * s.t. cell_corner_v[p0][*]=p1
	Tnode_fluid* access_extended_node(int local_face_id,int local_node_id);
	//Calculate the extended node through _node_ptr[n] and face f 
	void set_cell_connection();
	//Set _adjacent and _around
	void add_to_node();
	//Add this cell to its nodes

//------------------------------------------------------------------------------
//Access functions
//------------------------------------------------------------------------------
public:
	Tpolyhedron* G_cell_polyhedron(){return &_cell_polyhedron;};
	virtual void S_material_fraction(double fraction,int mat_id){};
	virtual void S_material_centroid(vec3D centroid,int mat_id){};
	virtual void S_material_base(TMAT_base* mat){};
	virtual void S_material_base(TMAT_base* mat[]){};
//------------------------------------------------------------------------------
//Member variables
//-------------------------------------------------------------------------------
protected:
	Tnode_fluid* _node_ptr[8];                   //Point to node formation
	Tcell_fluid_base* _adjacent[6];              //The adjacent cells
	vector<Tcell_fluid_base*> _around;           //Cells which have at least one common node
	double _cell_volume;                         //Cell volumes
	double _cell_length;                         //Cell characteristic length
	Tpolyhedron _cell_polyhedron;                //The cell polyhedron
	double _corner_mass[8];                      //The 8 corner mass
	double _corner_volume[8];                    //The 8 corner volume
	vec3D _corner_moment[8];                     //The 8 corner moment
	vec3D _corner_hourglass_force[8];            //The 8 corner force by hourglass
	vec3D _corner_viscosity_force[8];            //The 8 corner viscosity force by artificial

	vec3D _coor_min;
	vec3D _coor_max;
	double _volume_intersection;                 //Cell volume calculated by the accumulation of intersected portion
	int _flag;

	friend class Tcell_intersect;
	friend class Tcell_pure_fluid;
	friend class Tedge_viscosity;
	friend class TFEM_fluid_viscosity;
	friend class TFEM_viscosity;
	friend class Tbody_Lagrangian;
	friend class Tnode_fluid;
	friend class Tbody_ALE;
	friend class Tgrid_smooth_continum;
	friend class Tcell_intersect;
	friend class Tbody_explosion;
	friend class Tbody_MPM;
	friend class Tbody_ALEMPM;
	friend class Tcell_base;
	friend class Tbody_ALE_with_explosion;
};
#endif