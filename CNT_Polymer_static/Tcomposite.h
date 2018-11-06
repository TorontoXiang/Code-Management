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
	Tgrid_CNT _grid_CNT;

	void CG_iteration_initialization();
	//Use _dis=0 as the initial value of the CG iteration to set _F0 and
	//Initialize the variables in CG iteration, namely calculating _r and _p in _grid_Polymer

	void Calculate_F0_b();
	//Calculate _F0 and _b in _grid_Polymer
	void Calculate_Fp();
	//Calculate _Fp in _grid_Polymer

	void Calculate_Fd();
};
#endif

