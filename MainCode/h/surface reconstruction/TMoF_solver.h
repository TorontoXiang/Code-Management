#ifndef TMOF_SOLVER
#define TMOF_SPLVER
//Define the solver for MoF surface reconstrcution
#include "Tpolyhedron.h"
#include "Tvolume_equation_solver.h"
#include "public_function.h"
#include <vector>
using namespace std;

class TMoF_solver
{
	struct SMoF_input                //The needed information of a MoF process
	{
		Tpolyhedron* _poly;          //The target polyhedron
		double _fraction_reference;  //The fraction of the material
		vec3D _centroid_reference;   //The reference centroid of the material

		SMoF_input(Tpolyhedron* poly);
	};
//	struct SMoF2D_input
//	{

//	};
	struct Siterator_point          //The information of an iterator in line searching
	{
		vec2D _point;               //The (theta,phy) of this iterator
		double _d;                  //The plane constant
		double _value;              //The objective function
		vec2D _gradient;            //The objective function gradient of this iterator
		vec2D _d_derivative;        //The plane constant gradient of this iterator
	};
public:
	TMoF_solver(int maximal_num_material,int type,int dimension);
	void set_initial_condition(Tpolyhedron* cell_poly,double fraction_list[],vec3D centroid_list[],int num_material);
	//Set the initial cell polyhedron and material fractions and their centroids
	void MoF_for_pairs(Tpolyhedron* poly_in,double (&fraction)[2],vec3D (&centroid)[2],Tpolyhedron* (&sub_poly)[2]);
	//MoF reconstruction for two material pairs
	double MoF_surface_reconstrcution();
	//Perform the MoF surface reconstruction process
	//Return the discrepance of approximate centroid and reference centroid
	void Finalize();
	//Reset variables after a reconstruction

	//Access functions
	Tvolume_equation_solver* G_volume_solver(){return &_volume_equation_solver;};
	Tpolyhedron* G_result(int i){return _result[i];};
	int G_nun_non_convergence(){return _num_non_convergence;};
	int G_num_material(){return _num_material;};

	//----------------------------------------------------------------------------
	//                              Test functions
	//--------------------------------------------------------------------------
	vec2D derivative_comparison(Tpolyhedron* poly_in,double fraction,vec3D centroid);
	//Compare the partical derivatives of numerical and analytical
	vec2D patch_test(Tpolyhedron* poly_in,double fraction);
	//Test the MoF solver with given normal n and fraction
	double three_material_patch_test(Tpolyhedron* poly_in,double (&fraction)[2],vec2D (&angle)[2]);
	//Test the MoF solver for three materials plane
private:
	Tvolume_equation_solver _volume_equation_solver;  //The solver for volume equation
	int _maximal_num_material;                        //The maximal number of material in the computation
	int _num_material;                                //The number of material to be reconstructed
	vector<double> _fraction_list;                    //The volume fraction of each material
	vector<vec3D> _centroid_list;                     //The material centroid of each material
	vector<int> _material_id;                         //The id of candidate materials
	Tpolyhedron* _cell_poly;                          //The input cell polyhedron
	vector<Tpolyhedron*> _result;                     //The output pure material polyhedron 
	int _derivative_type;                             //The type of derivative calculation
	int _dimension;                                   //The dimension of the problem
	int _num_non_convergence;                         //The number of non-convergence cases

	vec2D calculate_initial_guess(SMoF_input &MoF_input);
	//Calculate the initial guess of the optimization
	vec3D calculate_iterator_info(SMoF_input &MoF_input,double theta,double phy,double d,vec2D &d_derivative);
	//Calculate the information at iterator point
	//d: plane constant at this iterator point
	//ouput._coor.x: value of objective function
	//ouput._coor.y: objective function derivative with respected to theta
	//ouput._coor.z: objective function derivative with respected to phy
	//d_derivative: plane constant derivative respected to theta and phy
	vec2D BFGS_iteration(SMoF_input &MoF_input,vec2D initial_geuss,bool &is_failed);
	//BFGS optimization process
	//Input the initial guess and output the optimal result
	//discrepance: The discrepance of centroids of below and above parts
	TMoF_solver::Siterator_point line_searching(SMoF_input &MoF_input,Siterator_point &iterator,vec2D direction,bool &is_failed);
	//Line searching in BFGS method
	//Input the former point and output the next point
	vec2D limit_angle(vec2D direction);
	//limit the angle into [0,pi],[0,2pi]
	void classify_mix_cell();
	//Classify the number of materials in the cell
	double three_material_trial(int i,Tpolyhedron* (&sub_poly)[3]);
	//The trial solution for 3 material
	//In the first division,material i is individual and materials (i+1)%3,(i+2)%3 are together
	//and they will be distinguished in the second division
	
};
#endif