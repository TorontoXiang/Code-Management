#include "Tcell_base.h"
Tcell_base::Tcell_base(int cell_id,int nGauss)
{
	_id=cell_id;
	_gausspoint.resize(nGauss);
	return;
}