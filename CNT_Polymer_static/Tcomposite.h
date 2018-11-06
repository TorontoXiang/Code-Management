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

	void Calculate_F0_b();
	//Calculate _F0 and _b in _grid_Polymer
	void Calculate_Fp();
	//Calculate _Fp in _grid_Polymer

	void Calculate_Fd();
};
#endif

