#ifndef TCELL_MIXED_FLUID
#define TCELL_MIXED_FLUID
//Define the class of mixed fluid
#include "Tcell_fluid_base.h"
class Tcell_mixed_fluid : public Tcell_fluid_base
{
public:
	Tcell_mixed_fluid(){};
	Tcell_mixed_fluid(int cell_id,int num_material,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat);
	Tcell_mixed_fluid(int cell_id,int num_material,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat[]);

//=======================================================================
//Virtual functions from Tcell_fluid_base
//=======================================================================
//-----------------------------------------------------------------------
//Functions for calculating corner force
//-----------------------------------------------------------------------
protected:
	virtual bool calculate_corner_stress_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3],double dt);
	//Calculate the corner force due to pressure
	//Closure model will be involved in Tcell_mix_fluid
	//Return whether it is failed in calculating corner_force from stress
	virtual void calculate_corner_hourglass_force(vec3D (&surfacial)[8][3],vec3D (&assistant)[8][3]);
	//Calculate the corner force due to hourglass
	//Multipul materials will be considered in Tcell_min_fluid
	virtual void set_pressure(double t);

//--------------------------------------------------------------------
//Functions for physical
//--------------------------------------------------------------------
public:
	virtual void calculate_corner_mass();
	//Calculate the corner mass by corner volume and density
	//Tcell_pure_fluid only use one material
	//Tcell_mix_fluid use multipul materials
	virtual bool update_state(double dt);
	//Update the material state at corrector step
	//Return whether error occurs in update process
	virtual void reset_state_after_remapping();
	virtual void reset_state_after_overlapping_exp();
	//Reset the cell state after remapping phase
	virtual void calculate_average_variable(double  &av_density,double &av_pressure,double &av_soundspeed);
	//Calculate the average density,pressure and soundspeed
	virtual void assemble_corner_inertance();
	//Assemble the corner inertance into node
	virtual double calculate_cell_internal_energy();
	//Calculate the internal energy in the cell
	virtual void reset_material_initial_state(int mid);
	//Set the initial state of material mid after remapping
	virtual void output_cell_monitor(ofstream& output,double t);
	//Output the cell information
//------------------------------------------------------------------------
//Functions in ALE process
//------------------------------------------------------------------------
public:
	virtual void calculate_gradient_in_cell();
	//Calculate the gradient of variables
	//Tcell_pure_fluid use the information of cell polyhedron
	//Tcell_mix_fluid use the information of material polyhedron
	virtual void limit_gradient_in_cell();
	//Limit the gradient by BJ limiter
	//Tcell_pure_fluid use the information of cell vertexes
	//Tcell_mix_fluid use the information of material polyhedron vertexes
	virtual void clear_cell_variable();
	//Clear the variable on the old cell after remapping
	virtual void accumulate_variables(int mid,double volume,vec3D int_x,double mass,double internal_energy);
	//Accumulate the variables of material from the intersection portion
	virtual void reset_geometry_initial_condition();
	//Reset the initial condition of material fraction and material centroid in mixed cell after overlapping
	virtual void reset_conservative_variable_initial_condition();
	//Reset the initial material mass and material internal energy of mixed cell before time integration
//------------------------------------------------------------------------
//Functions for multi-material
//------------------------------------------------------------------------
	virtual void update_material_centroid();
	//Update the material centroid before surface reconstruction
	virtual void delete_material_polyhedron();
	//Delete the material polyhedron after remapping
	virtual int is_surface_cell(int mid);
	//Determine whether it is a surface cell of material mid
	//If it is,return the number of surface piece of this cell
	virtual void add_surface_cell_surface(int mid,int** piece_list,int vertex_id_start,int &piece_id);
	//Add the surface piece of a surface cell into piece_list
	//mid:the material id
	//piece_list:the current piece_list
	//vertex_id_start:the start global vertex id of this cell
	//piece_id_start:the start global piece id of this cell
//-----------------------------------------------------------------------------
//Access functions
//-----------------------------------------------------------------------------
	virtual Tpolyhedron* G_material_polyhedron(int i){return _material_polyhedron[i];};
	//Return _cell_polyhedron for Tcell_pure_fluid
	virtual vec3D G_material_centroid(int i){return _material_centroid[i];};
	//Return cell centroid for Tcell_pure_fluid
	virtual vec3D G_g_ro(int i){return _g_ro[i];};
	//Return the density gradient of material i
	virtual vec3D G_g_roe(int i){return _g_roe[i];};
	//Return the density*internal of material i
	virtual double G_material_fraction(int i){return _material_fraction[i];};
	//Return the material fraction of material i
	virtual double G_average_pressure(){return _average_pressure;};
//------------------------------------------------------------------------
//Value functions
//------------------------------------------------------------------------
	virtual void S_material_polyhedron(int i,Tpolyhedron* poly_ptr){_material_polyhedron[i]=poly_ptr;};
	//Set the material polyhedron for material i
	virtual void S_material_fraction(double fraction,int mat_id){_material_fraction[mat_id]=fraction;};
	virtual void S_material_centroid(vec3D centroid,int mat_id){_material_centroid[mat_id]=centroid;};
	virtual void S_material_base(TMAT_base* mat[]);
//-----------------------------------------------------------------------------
//Member variables
//------------------------------------------------------------------------------
protected:
	vector<double> _material_mass;               //The materials' mass
	vector<double> _material_internal_energy;    //The materials' internal energy
	vector<vec3D> _g_ro;                         //The gradients' of material density
	vector<vec3D> _g_roe;                        //The gradients' of material density*internal energy
	vector<double> _material_fraction;           //The materials' volume fraction
	vector<double> _fraction_increment;          //The increment of material fraction in predict step
	vector<vec3D> _material_centroid;            //The materials' material centroid
	vector<vec3D> _material_centroid_std;        //The materials' material centroid in standard space
	vector<Tpolyhedron*> _material_polyhedron;   //The materials' polyhedron
	double _average_pressure;                    //The pressure from closure model

	int classify_mixed_cell(int (&id)[10]);
	//Determine the materials in the cell
	//Retrun the number of different materials in this cell
	//id[i]: id[i] is the material id of ith material in this cell
	void filter_fraction();
	//Modify the volume fraction if it is close to zero or one
};
#endif