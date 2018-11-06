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
	void Create_MKL_solver();
	//Create the MKL solver and do the numerical factorization
	Sstress** Calculate_cell_stress();
	//Calcaulte the stress at the 8 Gauss points of each cell
	void Output_Tecplot(string file_name);
	//Calculate the stress and output;
	void calculate_reacting_force();
	//Calculate the reacting force of each constraint
	void calculate_load_from_constraint();
	//Calculate _load_constraint
protected:
	vector<Snode> _node_list;
	vector<Tcell> _cell_list;
	int _nume, _nump;

	//Variables for solving equilibrium equation
	int _num_freedom_degree;                //The number of freedom degree
	int _num_non_zero;                      //The number of non-zero elements in the stifness matrix
	int _num_fixed;                         //The number of fiexed freedom degree
	vector<int> _ID;                        //The freedom degree id of every node at 3 direction
											//>0:  The freedom degree id
											//<0:  The displacement on this freedom degree is fixed
	int* _iK;
	int* _jK;
	double* _K;                             //Sparse stiffness matrix                       
	double* _load_constraint;               //The right hand side of from the displacement constraint
	double* _reacting_force;                //The reacting force from the boundary condition
	double* _dis;                           //The displacement at each freedom of degree
	TMKL_solver MKL_solver;                 //The MKL solver for this grid



	//The material property
	vector<Smat> _mat_list;

	void calculate_ID();
	//Calculate the ID vector and the freedom degree of the system
	void calculate_iK_jK();
	//Calculate iK and jK of the stifness matrix
	void assemble_equations();
	//Assemble the equations: calculate K and _load
	void calculate_cell_displacement(int cell_id,double(&d_cell)[8][3]);
	//Calculate the displacement of each cell
	void detect_boundary_cells();
	//Find the cells connect to node with boundary condition
};
class Tgrid_Polymer : public Tgrid
{
public:
	Tgrid_Polymer() {};
	void Create_MKL_solver();
	//Tgrid::Create_MKL_solver + allocate the memory for _load_total
	void Input_Polymer();
	//Input the regular polymer grid
	vec3D calculate_external_load(string face_name);
	//Calculate the total external load on the polymer boundary
	void Solving_equilibrium_equation();
	//Solving the equilibrium equation (The RHS is the sum of load from constraint and the CNT)
	void calculate_load_from_CNT(Tgrid_CNT* grid_CNT,int phase);
	//Calculate the force from a CNT grid
	//phase=0:load equals to _load
	//phase=1:load equals to _F0
	//phase=2:load equals to _Fp
	void calculate_cell_p(int cell_id, double(&p_cell)[8][3]);
	//Calculate cell p in CG iteration
	void initialize_load_total(int phase);
	//Set _load_total to 0
	void Calculate_b();
	//Calculate _b vector
	void Calculate_Kp();
	//Calculate Kp
	void Calculate_Kd();

	void Calculate_Ap();
	//Calculate (K+T)p

	void Calculate_Ad();

	void CG_initialize_p();
	void CG_initialize_r();
	//Initialization the CG iteration
	void CG_update();
	//Update the variables in CG iteration
	double calculate_diff();
	//Return the maximal value of _r
	void output_dis(ofstream& output);

	void output_p()
	{
		for (int i = 0; i < _num_freedom_degree; i++)
		{
			cout << i << " " << _p[i] << endl;
		};
	}
	void calculate_rd()
	{
		for (int i = 0; i < _num_freedom_degree; i++)
		{
			_rd[i] = _b[i] - _Ad[i];
			cout << i << " " << _r[i] << " " << _rd[i] << endl;
		};
	}
protected:
	double _x_min[3], _x_max[3];
	double _interval[3];
	int _nx, _ny, _nz;

	double* _load_total;                //Load from constriant and the CNT

	//Variables for CG iteration
	double* _r;
	double* _p;
	double* _Kp;
	double* _Fp;
	double* _F0;
	double* _b;
	double* _Ap;

	double* _Kd;
	double* _Fd;
	double* _Ad;
	double* _rd;

	void apply_regular_bc(Sregular_grid_bc& bc);

	friend class Tgrid_CNT;
};
class Tgrid_CNT : public Tgrid
{
public:
	Tgrid_CNT() {};
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
protected:
	void detect_boundary_nodes();
	//Detect the boundary nodes of the CNT

	friend class Tgrid_Polymer;
};
#endif
