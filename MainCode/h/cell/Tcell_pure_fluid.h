#ifndef TCELL_PURE_FLUID
#define TCELL_PURE_FLUID
//Define the class for pure fluid cell
#include "Tcell_fluid_base.h"
class Tcell_pure_fluid : public Tcell_fluid_base
{
public:
	Tcell_pure_fluid(){};
	Tcell_pure_fluid(int cell_id,Tnode_fluid* (&node_ptr)[8],TMAT_base* mat);

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
	//Reset the cell state after remapping phase
	virtual void calculate_average_variable(double  &av_density,double &av_pressure,double &av_soundspeed);
	//Calculate the average density,pressure and soundspeed
	virtual void assemble_corner_inertance();
	//Assemble the corner inertance into node
	virtual double calculate_cell_internal_energy();
	//Calculate the internal energy in the cell
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
	virtual void reset_conservative_variable_initial_condition();
//-----------------------------------------------------------------------------
//Access functions
//-----------------------------------------------------------------------------
	virtual Tpolyhedron* G_material_polyhedron(int i){return &_cell_polyhedron;};
	//Return _cell_polyhedron for Tcell_pure_fluid
	virtual vec3D G_material_centroid(int i){return _cell_polyhedron.G_centroid();};
	//Return cell centroid for Tcell_pure_fluid
	virtual vec3D G_g_ro(int i){return _g_ro;};
	//Return the density gradient of material i
	virtual vec3D G_g_roe(int i){return _g_roe;};
	//Return the density*internal of material i
	virtual void S_mass(int i,double mass){_material_mass=mass;};
	//Set the material mass
	virtual void S_material_base(TMAT_base* mat){_gausspoint[0]=mat;};
//-----------------------------------------------------------------------------
//Member variables
//------------------------------------------------------------------------------
protected:
	double _material_mass;               //The material mass
	double _material_internal_energy;    //The material internal energy
	vec3D _g_ro;                         //The gradient of material density
	vec3D _g_roe;                        //The gradient of material density*internal energy

};
#endif