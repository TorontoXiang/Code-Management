#ifndef TCOMPOSITE
#define TCOMPOSITE
#include "TFE_grid.h"
class Tcomposite
{
public:
	Tcomposite() {};
	void Input();
	//Input the polymer grid and CNT grids
private:
	Tgrid_Polymer _grid_Polymer;
	vector<Tgrid_CNT> _grid_CNT_list;
	
};
#endif

