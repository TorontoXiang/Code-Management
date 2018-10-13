#ifndef TGRID_SMOOTH_CONTINUM
#define TGRID_SMOOTH_CONTINUM
//Define the class of continum analogy for grid smooth
#include "Tgrid_smooth.h"
#include <vector>
#include "Tcell_8H.h"
using namespace std;
class Tgrid_smooth_continum : public Tgrid_smooth_base
{
public:
	Tgrid_smooth_continum(Tbody_ALE* the_body);

	//Initilize the object
	//void initialize();
	//Generate a new grid in this_body
	virtual void generate_new_grid();

protected:
	vector<Tcell_8H> _pseudo_cell_list;     //The pseudo brick cells
	int _nump;
	int _nume;
	int _num_freedom_degree;                //The number of freedom degree
	int _num_non_zero;                     //The number of non-zero elements in the stifness matrix
	int _num_fixed;                         //The number of freedom degree which has non-zero given displacement
	vector<int> _ID;                        //The freedom degree id of every node at 3 direction
	                                        //>0:  The freedom degree id
	                                        //-1:  The displacement on this freedom degree is zero
	                                        //-2:  The displacement on this freedom degree is given (for the conode)
	int* _iK;
	int* _jK;
	double* _K;                             //Sparse stiffness matrix                       
	double* _load;                          //The right hand side of the equilibrium equation

	void create_pseudo_cell_list();
	//Create the pseduo cell list
	void calculate_ID();
	//Calculate the ID vector and the freedom degree of the system
	void calculate_iK_jK();
	//Calculate iK and jK of the stifness matrix
	void assemble_equations();
	//Assemble the equations: calculate K and _load
	void solve_equations(double* dis);
	//Solve the equations by MKL PARDISO
	void set_initial_and_final_position();
	//Store the initial and final position
	void calculate_current_position_at_conode(int n,int num_move);
	//Calculate the position of co-node of current iteration
	//n:this is the n th iteration;num_move:number of total iterations
	void solidify_node();
	//Set the current position as the initial position in the next iteration
};
#endif