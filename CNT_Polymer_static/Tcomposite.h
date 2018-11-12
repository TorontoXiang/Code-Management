#ifndef TCOMPOSITE
#define TCOMPOSITE
#include "TFE_grid.h"
class Tcomposite
{
public:
	Tcomposite() {};
	void Input();
	//Input the polymer grid and CNT grids
	void CG_iteration();
	//The CG_iteration to solve the problem
	void output_result();
	
private:
	Tgrid_Polymer _grid_Polymer;
	int _num_CNT;
	vector<Tgrid_CNT> _CNT_list;

	void CG_iteration_initialization();
	//Use _dis=0 as the initial value of the CG iteration to set _F0 and
	//Initialize the variables in CG iteration, namely calculating _r and _p in _grid_Polymer
	void Calculate_Fp();
	//Calculate _Fp in _grid_Polymer
	void Calculate_reacting_force();
	//Calculate the reacting froce
	void input_CNT_list();
	//Input the CNT grid
};
#endif

