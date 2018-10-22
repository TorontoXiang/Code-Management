#include "Tnode.h"
#include "Tcell_shell.h"
#include "Tcell_base.h"
#include "TMATshell_elastic.h"
#include <iostream>
#include "Tregion.h"
#include <Windows.h>
#include "readin.h"
#include <omp.h>
//#include "Tgrid_smooth_continum.h"
#include"Tvolume_equation_solver.h"
#include "Tpolyhedron.h"
#include "data_structure.h"
#include <cmath>
#include "public_function.h"
#include "Tcell_fluid_base.h"
//#include "TMoF_solver.h"
using namespace std;
int main()
{
	omp_set_num_threads(16);
	Tregion region;
	if (!region.input_region())
	{
		cout<<"Error in openning input date!!"<<endl;
		system("Pause");
		exit(0);
	}
	//system("Pause");
	double t_begin=GetTickCount();
	//region.test();
	region.time_integration();
	double t_end=GetTickCount();
	int total_seconds=int((t_end-t_begin)/1000);
	int hour=total_seconds/3600;
	int min=(total_seconds-3600*hour)/60;
	int second=total_seconds-3600*hour-60*min;
	cout<<"Time consuming is "<<hour<<" hours "<<min<<" mins "<<second<<" seconds"<<endl;
	system("Pause");
	return 0;
}