#ifndef TFE_GRID
#define TFE_GRID
//Define the class of the FEM grid
#include "data_structure.h"
#include "Tcell.h"
#include <vector>
#include "readin.h"
#include "TMKL_solver.h"
#include <string>
#include "readin.h"
using namespace std;
class Tgrid
{
public:
	Sstress** Calculate_cell_stress();
	//Calcaulte the stress at the 8 Gauss points of each cell
	void Output_Tecplot(string file_name,int mid=-1);
	//Calculate the stress and output;
protected:
	vector<Snode> _node_list;
	vector<Tcell> _cell_list;
	int _nume, _nump;

	//Variables for solving equilibrium equation
	int _num_freedom_degree;                //The number of freedom degree
	int _num_fixed;                         //The number of fiexed freedom degree
	vector<int> _ID;                        //The freedom degree id of every node at 3 direction
											//>0:  The freedom degree id
											//<0:  The displacement on this freedom degree is fixed
	int* _iK;
	int* _jK;
	double* _K;                             //Sparse stiffness matrix         
	double* _p;	                            //_p can be used as:
											//1. The p in CG iteration in Tgrid_Polymer
											//2. The RHS of the equilibrium equation in Tgrid_CNT
											//3. Temporary store the RHS from boundary condition in Tgrid_Polymer when initilizing the CG iteration                             
	double* _reacting_force;                //The reacting force from the boundary condition
	                                        //_num_fixed dimension for Tgrid_Polymer
	                                        //3*_nump dimension for Tgrid_CNT
	double* _dis;                           //The displacement at each freedom of degree



	//The material property
	vector<Smat> _mat_list;

	void calculate_ID();
	//Calculate the ID vector and the freedom degree of the system
	void calculate_cell_displacement(int cell_id,double(&d_cell)[8][3]);
	//Calculate the displacement of each cell
	void detect_boundary_cells();
	//Find the cells connect to node with boundary condition
};
class Tgrid_CNT;
class Tgrid_Polymer : public Tgrid
{
public:
	Tgrid_Polymer() {};
	void Calculate_stiffness_matrix();
	//Only calculate the sparse stiffness matrix but does not perform the numerical factorization
	//Allocate the memory for the variables in CG iteration
	void Input_Polymer();
	//Input the regular polymer grid
	vec3D calculate_external_load(string face_name);
	//Calculate the total external load on the polymer boundary
	void calculate_reacting_force();
	//Calculate the reacting force of each constraint
	void calculate_load_from_constraint();
	//Calculate _load_constraint
	void assemble_reacting_froce_from_CNT(Tgrid_CNT* grid_CNT);
	//Assemble the reacting force from CNT

	//Functions for CG iteration
	void assemble_Fp(Tgrid_CNT* grid_CNT);
	//Assemble the Fp from a CNT grid
	void calculate_cell_p(int cell_id, double(&p_cell)[8][3]);
	//Calculate cell p in CG iteration
	void initialize_Fp();
	//Set _Fp to 0
	void Set_F0();
	//Set _F0=_Fp in the initialization of CG iteration
	void Calculate_Kp();
	//Calculate Kp
	void CG_initialization();
	//Initialization the CG iteration
	void CG_update();
	//Update the variables in CG iteration
	double calculate_diff();
	//Return the maximal value of _r

protected:
	double _x_min[3], _x_max[3];
	double _interval[3];
	int _nx, _ny, _nz;

	//Variables for CG iteration
	double* _r;
	double* _Kp;
	double* _Fp;
	double* _F0;

	void apply_regular_bc(Sregular_grid_bc& bc);
	void calculate_iK_jK();
	//Calculate iK and jK of the stifness matrix
	void assemble_equations();
	//Assemble the equations: calculate K matrix
	friend class Tgrid_CNT;
};
class Tgrid_CNT : public Tgrid
{
public:
	Tgrid_CNT() {};
	void Create_MKL_solver();
	//Create the MKL solver and do the numerical factorization
	void Input_CNT(Skeyword& keyword);
	//Generate this CNT from keyword
	void calculate_CNT_location(Tgrid_Polymer* grid_polymer);
	//Calculate the location of the surface node the CNT grid
	void calculate_CNT_boundary_displacement(Tgrid_Polymer* grid_polymer);
	//Calculate the boundary displacement of a CNT
	void calculate_CNT_boundary_p(Tgrid_Polymer* grid_polymer);
	//Calculate the boundary p a CNT in CG_iteration
	void Solving_equilibrium_equation();
	//Solving the equilibrium equation (The RHS is only the load from constraint)
	void calculate_reacting_force();
	//Calculate the reacting force of each constraint
	void calculate_load_from_constraint();
	//Calculate _load_constraint
	void output_tecplot(ofstream& output);
	//Output the results of CNT
protected:
	void detect_boundary_nodes();
	//Detect the boundary nodes of the CNT

	TMKL_solver MKL_solver;                 //The MKL solver for this grid
	double* _dis_whole;                      //The displacement of all

	int* _iKB;
	int* _jKB;
	double* _KB;                           //The sparse matrix for calculating the reacting force

	void calculate_iK_jK();
	//Calculate iK and jK of the stifness matrix
	void assemble_equations();
	//Assemble the equations: calculate K and KB matrix
	friend class Tgrid_Polymer;
};
#endif
