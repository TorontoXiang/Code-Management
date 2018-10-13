#ifndef TVOLUME_EQUATION_SOLVER
#define TVOLUME_EQUATION_SOLVER
//Define the solver for the local volume enforcement
#include "Tpolyhedron.h"
#include "data_structure.h"
class Tvolume_equation_solver
{
public:
	Tvolume_equation_solver():_poly(NULL),_num_Newton_fails(0),_num_solving(0){};
	double calculate_plane_constant(vec3D &n,double fraction,double predict=1e10);
	double calculate_plane_constant_advanced(vec3D &n,double fraction,double predict=1e10);
	//Calcualte the d in the plane nx+d=0 with a reference volume fraction
	//predict: the predicted constant which will be a lagre number if it is not given
	void set_target_polyhedron(Tpolyhedron* poly);
	//Set the target polyehdron
	vec3D calculate_fraction_below(vec3D &n,double d);    
	//Calculate volume of the polyhedron below the cutting plane (first component of the result)
	//Calculate the derivative with respect to d (the area of the cutting ring) ( second component of the result)
	//void clear_sub_volume();
	////Cleat the sub_volume after this polyhedron has operated
	int G_num_Newton_fails(){return _num_Newton_fails;};
	int G_num_solving(){return _num_solving;};
private:
	Tpolyhedron* _poly;            //The polyhedron to be operated
	int _num_Newton_fails;         //The number of failure in solving the equation
	double _num_solving;           //The number of solving the equation

	bool terminate(double d_min,double d_max);
	//Determine whether d_min and d_max are in a same interval of _altitude_list
	int find_altitude_position(double altitude);
	//Find the position of altitude in _altitude_list
	
};
#endif
